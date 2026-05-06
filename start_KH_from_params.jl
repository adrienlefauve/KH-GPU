using Base.Threads
using Dates
using CUDA
using JLD2
ENV["GKSwstype"] = "100"
using Plots
gr()
using Printf
using Oceananigans
using Oceananigans.OutputWriters
using Interpolations
using Statistics
using LinearAlgebra
using SpecialFunctions


@info "-----------------------------------------------------------------------------------------------"
@info "-----------------------------------------------------------------------------------------------"

# Log GPU information
CUDA.NVML.nvmlInit()
handle = CUDA.NVML.device()
props = CUDA.device()
@info "GPU name: $(CUDA.name(props))"
sm_count = CUDA.attribute(props, CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)
@info "Streaming Multiprocessors: $sm_count"

# Assumes Julia was launched from the repo root (the directory containing this script)
include("./code/SSF.jl")
include("./code/ddz.jl")
include("./code/ddz2.jl")
include("./code/ddz4.jl")
include("./lin_stab_from_params.jl")
include("./KH_from_params.jl")
include("./plot_KH.jl")
include("./plot_KH_b_only.jl")  # alternative lighter plotting (b only)

# Read parameter arrays
include("./params.jl")
params = Dict(
    :Re => Re_values,
    :Ri => Ri_values,
    :Pr => Pr_values,
    :Lz => Lz_values,
    :Nx => Nx_values,
    :Nz => Nz_values,
    :efgm => efgm_values,
    :enoise => enoise_values,
    :T => T_values,
)

# Find the maximum length across all parameter arrays
nruns = maximum(length.(values(params)))

# Validate: each array must have length == 1 or == nruns
for (name, arr) in params
    if !(length(arr) == 1 || length(arr) == nruns)
        error("Parameter $name has length $(length(arr)), which is not 1 or $nruns")
    end
end

# Helper to get i-th value or repeat a single value
getval(a, i) = length(a) == 1 ? a[1] : a[i]

# Main loop
runs_start_time = time_ns()
for irun in 1:nruns
    # Extract parameters for this run 
    @info "STARTING RUN $irun of $nruns ..."
    Re = getval(params[:Re], irun)
    Ri = getval(params[:Ri], irun)
    Pr = getval(params[:Pr], irun)
    Lz = getval(params[:Lz], irun)
    Nx = getval(params[:Nx], irun)
    Nz = getval(params[:Nz], irun)
    efgm = getval(params[:efgm], irun)
    enoise = getval(params[:enoise], irun)
    T = getval(params[:T], irun)

    # Set output paths
    today_str = Dates.format(today(), "yyyy-mm-dd")
    today_str = replace(today_str, "-" => "") # to match the required format
    base_path = pwd() # should be the KH_GPU_sweep directory
    out_path = joinpath(base_path, "out")
    lin_stab_out_name = @sprintf("KH_lin_stab_Re=%.1e_Ri=%0.2f_Pr=%d_Lz=%d_efgm=%.1e",
        Re, Ri, Pr, Lz, efgm)
    lin_stab_out_file = joinpath(out_path, lin_stab_out_name)
    dns_out_name = @sprintf("KH_%s_Re=%.1e_Ri=%0.2f_Pr=%d_Lz=%d_Nx=%d_Nz=%d_efgm=%.1e_enoise=%.1e_T=%d",
        today_str, Re, Ri, Pr, Lz, Nx, Nz, efgm, enoise, T)
    dns_out_file = joinpath(out_path, dns_out_name)

    # If needed, run linear stability analysis and save results in lin_stab_out_file
    if isfile(lin_stab_out_file * ".jld2")
        @info "Linear stability results already exist, skipping to next step..."
    else
        lin_stab_from_params(Re, Ri, Pr, Lz, efgm, lin_stab_out_file)
    end

    # Run DNS simulation using lin stab results and save results in dns_out_file
    KH_from_params(Re, Ri, Pr, Lz, Nx, Nz, enoise, T, lin_stab_out_file, dns_out_file)

    # Generate plots after simulation completion
    plot_KH(dns_out_file)
    @info "-----------------------------------------------"
end

runs_time = (time_ns() - runs_start_time) * 1e-9

@info "-----------------------------------------------------------------------------------------------"
@info "ALL RUNS COMPLETED IN $(prettytime(runs_time))"
@info "-----------------------------------------------------------------------------------------------"