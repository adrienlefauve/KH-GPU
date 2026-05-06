# Reads in the JLD2 output written by KH_from_params.jl and saves an animation of buoyancy and log10(chi).
function plot_KH(filename)

    @info "Generating plots for buoyancy and mixing..."
    plot_start_time = time_ns()

    # Read in the first iteration.  We do this to load the grid
    # filename * ".jld2" concatenates the extension to the end of the filename
    u_ic = FieldTimeSeries(filename * ".jld2", "u", iterations = 0)
    # v_ic = FieldTimeSeries(filename * ".jld2", "v", iterations = 0)
    w_ic = FieldTimeSeries(filename * ".jld2", "w", iterations = 0)
    b_ic = FieldTimeSeries(filename * ".jld2", "b", iterations = 0)
    # ω_ic = FieldTimeSeries(filename * ".jld2", "ω", iterations = 0)
    χ_ic = FieldTimeSeries(filename * ".jld2", "χ", iterations = 0)
    # ϵ_ic = FieldTimeSeries(filename * ".jld2", "ϵ", iterations = 0)
    # KE = FieldTimeSeries(filename * ".jld2", "E", iterations = 0)

    # Load in coordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xu, yu, zu = nodes(u_ic)
    # xv, yv, zv = nodes(v_ic)
    xw, yw, zw = nodes(w_ic)
    xb, yb, zb = nodes(b_ic)
    # xω, yω, zω = nodes(ω_ic)
    xχ, yχ, zχ = nodes(χ_ic)
    # xϵ, yϵ, zϵ = nodes(ϵ_ic)
    xlims = (minimum(xb), maximum(xb))
    ylims = (minimum(zb), maximum(zb))


    # Now, open the file with our data
    file_xz = jldopen(filename * ".jld2")

    # Extract a vector of iterations
    iterations = parse.(Int, keys(file_xz["timeseries/t"]))

    @info "Making an animation from saved data..."

    # First pass: find global min/max for consistent color limits
    b_global_min, b_global_max = Inf, -Inf
    χ_global_min, χ_global_max = Inf, -Inf

    for iter in iterations
        b_data = file_xz["timeseries/b/$iter"][:, 1, :]
        χ_data = file_xz["timeseries/χ/$iter"][:, 1, :]
        χ_log_data = log10.(abs.(χ_data) .+ 1e-10)
        
        b_global_min = min(b_global_min, minimum(b_data))
        b_global_max = max(b_global_max, maximum(b_data))
        χ_global_min = min(χ_global_min, minimum(χ_log_data))
        χ_global_max = max(χ_global_max, maximum(χ_log_data))
    end

    times = Float64[]

    # Here, we loop over all iterations
    anim = @animate for (i, iter) in enumerate(iterations)

        # @info "Drawing frame $i from iteration $iter..."

        # Load data for this iteration
        b_xz = file_xz["timeseries/b/$iter"][:, 1, :];
        χ_xz = file_xz["timeseries/χ/$iter"][:, 1, :];
        t = file_xz["timeseries/t/$iter"];

        # Store values for time tracking
        push!(times, t)

        # Create titles
        b_title = @sprintf("negative b, t = %s", round(t));
        χ_title = @sprintf("log₁₀(χ)");
        
        # Use global limits for consistent visualization across all frames
        b_lims = (b_global_min, b_global_max)
        χ_lims = (χ_global_min, χ_global_max)
        x_lims = (minimum(xb), maximum(xb))
        z_lims = (minimum(zb), maximum(zb))
        
        # Create plots with Plots.jl syntax
        b_xz_plot = heatmap(xb, zb, -b_xz', color = :viridis, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = x_lims, ylims = z_lims, clims = b_lims, title=b_title); 
        χ_xz_plot = heatmap(xχ, zχ, log10.(abs.(χ_xz') .+ 1e-10), color = :hot, xlabel = "x", ylabel = "z", aspect_ratio = :equal,  xlims = x_lims, ylims = z_lims, clims = χ_lims, title=χ_title);
        
        # Combine the sub-plots into a single figure (vertical layout)
        plot(b_xz_plot, χ_xz_plot,
            layout = (2, 1),
            size = (1200, 1200))
    end

    close(file_xz)

    # Save the animation to a file
    mp4(anim, filename*".mp4", fps = 10) 

    # Display plot time
    plot_time = (time_ns() - plot_start_time) * 1e-9
    @info "   Plots completed in $(prettytime(plot_time))"
end  
