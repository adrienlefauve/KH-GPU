function KH_from_params(Re, Ri, Pr, Lz, Nx, Nz, enoise, T, lin_stab_out_file, dns_out_file)

    # Load initial condition which includes flow parameters
    @load lin_stab_out_file * ".jld2" sigma sigma_max K k_max u_hat w_hat b_hat U B z Ri Re Pr iBC1 iBCN Lz

    Lx = 2 * pi / k_max # one billow
    stretch = 1

    @info "Grid setup..."
    grid_start_time = time_ns()
    if stretch == 1
        # Grid Generation with Hyperbolic Stretching
        stretching = 5.0  # Controls refinement intensity at center (the higher the denser around center)

        # Normalized height coordinate (-1 to 1)
        h(k) = (2 * k - Nz - 1) / (Nz - 1)

        # Center-focused stretching function using tanh
        # This creates a hyperbolic tangent distribution that clusters points near z=0
        d = 0.3
        Σ(k) = (1 - d) .* sinh(stretching * h(k)) / sinh(stretching) + d * h(k)

        # Generating function (scaled to domain size)
        # Maps the stretched coordinates to physical space [-Lz/2, Lz/2]
        z_faces(k) = (Lz / 2) * Σ(k)

        # Generate grid points
        zz = [z_faces(k) for k in 0:Nz]

        # Compute grid spacing (always positive)
        Δz = diff(zz)

        # Display grid size information
        # @info "   Minimum grid spacing: $(round(minimum(Δz), sigdigits=3))"
        # @info "   Maximum grid spacing: $(round(maximum(Δz), sigdigits=3))"
        @info "   Vertical grid spacing ratio (max/min): $(round(maximum(Δz)/minimum(Δz), sigdigits=3))"
    else
        zz = (-Lz / 2, Lz / 2)
    end
    grid = RectilinearGrid(GPU(), size=(Nx, Nz), x=(-Lx / 2, Lx / 2), z=zz,
        topology=(Oceananigans.Periodic, Oceananigans.Flat, Oceananigans.Bounded))
    grid_time = (time_ns() - grid_start_time) * 1e-9
    @info "   completed in $(prettytime(grid_time))"
    # No boundary conditions explicitly set - the BCs will default to free-slip and no-flux in z

    # Now, define a 'model' where we specify the grid, advection scheme, bcs, and other settings
    @info "Model setup..."
    model_start_time = time_ns()
    buoyancy = BuoyancyForce(BuoyancyTracer(), gravity_unit_vector=[0, 0, -1])
    model = NonhydrostaticModel(; grid,
        advection=WENO(),  # WENO advection scheme (more accurate but slower than e.g. UpwindBiasedFifthOrder())
        timestepper=:RungeKutta3, # 3rd order Runge-Kutta time-stepping
        tracers=(:b),  # name(s) of any tracers; here b is buoyancy
        buoyancy=buoyancy, # buoyancy acts on momentum. Use BuoyancyTracer() (rather than Buoyancy(model=BuoyancyTracer())) for compatibility with current Oceananigans on GPU.
        closure=(ScalarDiffusivity(ν=1 / Re, κ=1 / (Re * Pr))),  # sets kinematic viscosity and diffusivty, here just 1/Re since we are solving the non-dimensional equations 
        coriolis=nothing # this line tells the mdoel not to include system rotation (no Coriolis acceleration)
    )

    # Set initial conditions with either FGM from the linear stability shifted so that inflection point of braid is at x=0
    ϕ_b = angle.(b_hat)
    z_idx = searchsortedlast(z, 0)
    if z[z_idx] != 0
        t = -z[z_idx] / (z[z_idx+1] - z[z_idx])      # fractional distance between x1 and x2
        shift = ϕ_b[z_idx] + t * (ϕ_b[z_idx+1] - ϕ_b[z_idx])
    else
        shift = ϕ_b[z_idx]
    end

    u_hat_interp = linear_interpolation(z, u_hat, extrapolation_bc=Interpolations.Flat())
    w_hat_interp = linear_interpolation(z, w_hat, extrapolation_bc=Interpolations.Flat())
    b_hat_interp = linear_interpolation(z, b_hat, extrapolation_bc=Interpolations.Flat())

    u_init(x, z) = real(u_hat_interp(z) * exp.(im * k_max .* x .- im * shift .+ im * pi / 2))
    w_init(x, z) = real(w_hat_interp(z) * exp.(im * k_max .* x .- im * shift .+ im * pi / 2))
    b_init(x, z) = real(b_hat_interp(z) * exp.(im * k_max .* x .- im * shift .+ im * pi / 2))

    # Here, we start with a tanh function for buoyancy and add an FGM perturbation + a random perturbations
    uᵢ(x, z) = tanh.(z) + u_init(x, z) + sqrt(enoise) * randn()
    vᵢ(x, z) = 0
    wᵢ(x, z) = w_init(x, z) + sqrt(enoise) * randn()
    bᵢ(x, z) = Ri * tanh.(z) + b_init(x, z) + Ri * sqrt(enoise) * randn()

    # Send the initial conditions to the model to initialize the variables
    set!(model, u=uᵢ, v=vᵢ, w=wᵢ, b=bᵢ)

    # Decide on the initial timestep based on the grid size and maximum velocity
    max_Δt = 1e-2 # maximum allowable timestep
    Δx = Lx / Nx
    cfl_initial_target = 0.3
    U_max = 1.0  # this doesn't include the initial perturbations
    initial_Δt = cfl_initial_target * Δx / U_max
    if initial_Δt > max_Δt
        initial_Δt = max_Δt
    end

    # Now, we create a 'simulation' to run the model for a specified length of time
    simulation = Simulation(model, Δt=initial_Δt, stop_time=T)

    # TimeStepWizard is a callback that adjusts the timestep during the simulation
    wizard = TimeStepWizard(cfl=0.8, diffusive_cfl=0.8, max_change=1.2, max_Δt=max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # Progress messenger
    start_time = time_ns()
    progress(sim) = @printf("i: % 6d, sim time: % 7s, wall time: % 7s, Δt: %.2e, CFL: %.2e, DCFL: %.2e, remaining: % 7s\n",
        sim.model.clock.iteration,
        sim.model.clock.time,
        prettytime(1e-9 * (time_ns() - start_time)),
        sim.Δt,
        AdvectiveCFL(sim.Δt)(sim.model),
        DiffusiveCFL(sim.Δt)(sim.model),
        prettytime((T - sim.model.clock.time) * (1e-9 * (time_ns() - start_time) / sim.model.clock.time)))
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))

    model_time = (time_ns() - model_start_time) * 1e-9
    @info "   completed in $(prettytime(model_time))"

    # Outputs
    @info "Output setup..."
    u, v, w = model.velocities # unpack velocity `Field`s
    b = model.tracers.b # extract the buoyancy
    # Now, calculate secondary quantities
    # Oceananigans has functions to calculate derivatives on the model grid˚
    outputs_start_time = time_ns()
    # ω = ∂z(u) - ∂x(w) # The spanwise vorticity (uncomment to also save)
    χ = (1 / (Re * Pr)) * (∂x(b)^2 + ∂z(b)^2) # The dissipation rate of buoyancy variance
    ϵ = (1 / Re) * (2 * ∂x(u)^2 + 2 * ∂z(w)^2 + ∂x(w)^2 + ∂z(u)^2 + 2 * ∂z(u) * ∂x(w)) # The dissipation rate of kinetic energy, using the full strain rate tensor
    outputs_time = (time_ns() - outputs_start_time) * 1e-9
    @info "   completed in $(prettytime(outputs_time))"

    # Set up the output file
    simulation.output_writers[:xz_slices] =
        JLD2Writer(model, (; u, w, b, χ, ϵ),
            filename=dns_out_file * ".jld2",
            indices=(:, 1, :),
            schedule=TimeInterval(0.5),  # save every 0.5 time unit
            with_halos=false,
            overwrite_existing=true)

    # Now, run the simulation
    simulation_start_time = time_ns()
    run!(simulation)
    simulation_time = (time_ns() - simulation_start_time) * 1e-9
    n_iter = simulation.model.clock.iteration
    time_per_iter = simulation_time / n_iter
    @info "   simulation completed in $(prettytime(simulation_time))"
    @info "   total iterations: $n_iter"
    @info "   average time per iteration: $(round(time_per_iter, digits=4)) seconds"
    @info "   results saved in $dns_out_file"
end