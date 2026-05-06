# Reads in the JLD2 output written by KH_from_params.jl and saves a buoyancy-only animation
# (lighter alternative to plot_KH.jl which also plots log10(chi)).
function plot_KH_b_only(filename)

    @info "Generating plots..."
    plot_start_time = time_ns()

    # Read in the first iteration.  We do this to load the grid
    # filename * ".jld2" concatenates the extension to the end of the filename
    b_ic = FieldTimeSeries(filename * ".jld2", "b", iterations = 0)

    # Load in coordinate arrays
    # We do this separately for each variable since Oceananigans uses a staggered grid
    xb, yb, zb = nodes(b_ic)
    xlims = (minimum(xb), maximum(xb))
    ylims = (minimum(zb), maximum(zb))


    # Now, open the file with our data
    file_xz = jldopen(filename * ".jld2")

    # Extract a vector of iterations
    iterations = parse.(Int, keys(file_xz["timeseries/t"]))

    @info "Making an animation from saved data..."

    # First pass: find global min/max for consistent color limits
    b_global_min, b_global_max = Inf, -Inf


    for iter in iterations
        b_data = file_xz["timeseries/b/$iter"][:, 1, :]
        
        b_global_min = min(b_global_min, minimum(b_data))
        b_global_max = max(b_global_max, maximum(b_data))
    end

    times = Float64[]

    # Here, we loop over all iterations
    anim = @animate for (i, iter) in enumerate(iterations)

        # @info "Drawing frame $i from iteration $iter..."

        # Load data for this iteration
        b_xz = file_xz["timeseries/b/$iter"][:, 1, :];
        t = file_xz["timeseries/t/$iter"];

        # Store values for time tracking
        push!(times, t)

        # Create titles
        b_title = @sprintf("b, t = %.2f", t)

        
        # Use global limits for consistent visualization across all frames
        b_lims = (b_global_min, b_global_max)
        x_lims = (minimum(xb), maximum(xb))
        z_lims = (minimum(zb), maximum(zb))
        
        # Create plots with Plots.jl syntax
        b_xz_plot = heatmap(xb, zb, b_xz', color = :viridis, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlims = x_lims, ylims = z_lims, clims = b_lims, title=b_title); 
  
        # Combine the sub-plots into a single figure (vertical layout)
        plot(b_xz_plot,
            size = (1200, 1200))
    end

    close(file_xz)

    # Save the animation to a file
    mp4(anim, filename*".mp4", fps = 10) 

    # Display plot time
    plot_time = (time_ns() - plot_start_time) * 1e-9
    @info "   Plots completed in $(prettytime(plot_time))"
end  
