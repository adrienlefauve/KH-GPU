function lin_stab_from_params(Re, Ri, Pr, Lz, fgm, lin_stab_out_file)
    Nz = 600
    z = LinRange(-Lz/2, Lz/2, Nz)
    U = tanh.(z)
    B = Ri * tanh.(z)

    @info "Linear stability computation..."
    lin_stab_start_time = time_ns()

    K = LinRange(0.35, 0.50, 16) # wavenumber range to sweep through

    iBC1=[0, 0] # changed from [0,1] check [0,0]
    iBCN=[0, 0]

    # initialise outputs
    sigma = complex(zeros(length(K)))
    lambda_w = complex(zeros(Nz, length(K)))
    lambda_b = complex(zeros(Nz, length(K)))

    Threads.@threads for k in eachindex(K)
        k_val = K[k]
        (sigma[k], lambda_w[:, k], lambda_b[:, k]) = SSF(z, U, B, k_val, 0, 1/Re, 1/(Pr*Re), iBC1, iBCN, 1)
    end

    k_index_max = argmax(real.(sigma)) # index of the k with the maximum real growth rate
    sigma_max = real.(sigma[k_index_max])
    k_max = K[k_index_max]
    w_hat = lambda_w[:, k_index_max]
    b_hat = lambda_b[:, k_index_max]
    # Compute uhat by divergence free condition
    D1 = ddz(z)  # 1st derivative matrix with 1-sided boundary terms
    u_hat = im / k_max * D1 * w_hat 

    nz1 = floor(Int,Nz/10)+1
    nz2 = floor(Int,Nz*9/10)

    # Rescale to a target energy
    background_energy = mean(1/2 *( U[nz1:nz2].^2+B[nz1:nz2].^2))
    mean_perturbation_energy = mean(1/2 *( abs2.(u_hat[nz1:nz2]) + abs2.(w_hat[nz1:nz2])+ abs2.(b_hat[nz1:nz2])))
    target_energy =fgm*background_energy
    scaling_factor = sqrt(target_energy / mean_perturbation_energy)

    w_hat*=scaling_factor
    u_hat*=scaling_factor
    b_hat*=scaling_factor

    # write output file
    Nz_lin_stab = Nz
    @save lin_stab_out_file*".jld2" sigma sigma_max K k_max u_hat w_hat b_hat U B z Ri Re Pr iBC1 iBCN Lz Nz_lin_stab
    lin_stab_time = (time_ns() - lin_stab_start_time) * 1e-9
    @info "   completed in $(prettytime(lin_stab_time))"
    @info "   results saved in $lin_stab_out_file"

end