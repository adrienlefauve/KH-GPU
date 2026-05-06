using Oceananigans
using NCDatasets
using FilePathsBase: basename
using Printf

# Function to export snapshots from a JLD2 file containing u,w,chi,epsilon at various times 
# to a directory of NetCDF files at specified times. This function has been optimised to handle large datasets
# without running out of memory by processing one field at a time and freeing memory explicitly.
# It also logs memory usage at key points to help diagnose any issues.

function export_jld2_to_netcdf_snapshots(filename::String, times)

    function log_memory_usage(tag="")
        try
            path = "/proc/self/status"
            open(path, "r") do io
                for line in eachline(io)
                    if startswith(line, "VmRSS:")
                        parts = split(strip(line))
                        rss_kb = parse(Int, parts[2])
                        rss_gb = rss_kb / 1024 / 1024
                        total_gb = Sys.total_memory() / 1024^3
                        @info "MEM $tag: $(round(rss_gb, digits=1)) GiB ($(round(100 * rss_gb / total_gb, digits=1))% of node)"
                        return
                    end
                end
                @warn "MEM $tag: VmRSS line not found in /proc/self/status"
            end
        catch e
            @warn "MEM $tag: error reading /proc/self/status ($e)"
        end
    end
    log_memory_usage("start of function")

    times_to_export = Float64.(collect(times))

    # Prepare output folder
    base = splitext(basename(filename))[1]
    outdir = joinpath(dirname(filename), base * "_netcdf")
    isdir(outdir) || mkpath(outdir)

    # Probe a single field just to get grid + time info
    b_fts = FieldTimeSeries(filename, "b")
    all_times = b_fts.times
    grid = b_fts.grid
    b_fts = nothing  # saves memory
    GC.gc()  # force garbage collection to free memory

    x_nodes, y_nodes, z_nodes = nodes(grid, Center(), Center(), Center(); reshape=false)
    x = vec(x_nodes[:, 1, 1])
    z = vec(z_nodes[:, 1, 1])
    y = [0.0]  # 2D assumption

    println("Exporting NETCDF snapshots to: $outdir")

    for t_req in times_to_export
        idx = findfirst(isequal(t_req), all_times)
        if isnothing(idx)
            println("   !!! requested time $t_req not found. Skipping.")
            continue
        end

        i = idx::Int
        t = all_times[i]
        log_memory_usage("start of t=$t")
        outnc = joinpath(outdir, @sprintf("snapshot_t%07.2f.nc", t))
        isfile(outnc) && rm(outnc; force=true)

        ds = NCDataset(outnc, "c", format=:netcdf4_classic)

        # Define dims
        defDim(ds, "x", length(x))
        defDim(ds, "y", length(y))
        defDim(ds, "z", length(z))

        defVar(ds, "x", Float32, ("x",); deflatelevel=0)
        defVar(ds, "y", Float32, ("y",); deflatelevel=0)
        defVar(ds, "z", Float32, ("z",); deflatelevel=0)
        ds["x"][:] = Float32.(x)
        ds["y"][:] = Float32.(y)
        ds["z"][:] = Float32.(z)

        defVar(ds, "t", Float32, (); deflatelevel=0)
        ds["t"][] = Float32(t)

        # Process fields one at a time — no cache buildup
        for fname in ["b", "u", "w", "ϵ", "χ"]
            fts = FieldTimeSeries(filename, fname)
            data = Array(fts[i])
            fts = nothing

            # Handle staggered grid if necessary
            if fname == "u" && size(data, 1) == length(x) + 1
                data = data[1:end-1, :, :]
            elseif fname == "w" && size(data, 3) == length(z) + 1
                data = data[:, :, 1:end-1]
            end

            data_nc = permutedims(Float32.(data), (3, 2, 1))
            data = nothing  # drop before writing to NetCDF

            safe_name = fname == "ϵ" ? "eps" : fname == "χ" ? "chi" : fname

            defVar(ds, safe_name, Float32, ("z", "y", "x"); deflatelevel=0)
            ds[safe_name][:, :, :] = data_nc
            data_nc = nothing

            GC.gc()
            sleep(0.2)
            GC.gc()
            sleep(0.2)
            log_memory_usage("   after $fname at t=$t")
        end

        close(ds)
        println("Saved snapshot t=$t")
        GC.gc()
    end

    println("Done")
end