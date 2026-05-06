function ddz_nonuni(z)
    # First derivative matrix for independent variable z.
    # 2nd order centered differences.
    # Use one-sided derivatives at boundaries.
    # adapted to take z not evenly spaced

    N = length(z)

    d = zeros(N, N)
    for n in 2:N-1
        del1 = z[n]-z[n-1];
        del2 = z[n+1]-z[n];
        d[n, n-1] = -del2/(2*del1*del2)
        d[n,n] = (del1-del2)/(2*del1*del2)
        d[n, n+1] = del1/(2*del1*del2)
    end
    del2 = z[2]-z[1];
    del3 = z[3]-z[2];
    d[1, 1] = -(2*del2*del3+del3^2)/(del2*del3*(del2+del3))
    d[1, 2] = (del2+del3)^2/(del2*del3*(del2+del3))
    d[1, 3] = -del2^2/(del2*del3*(del2+del3))

    delN = z[N]-z[N-1];
    delm1 = z[N-1]-z[N-2];
    d[N, N] = (2*delN*delm1+delm1^2)/(delN*delm1*(delN+delm1))
    d[N, N-1] = -(delN+delm1)^2/(delN*delm1*(delN+delm1))
    d[N, N-2] = delN^2/(delN*delm1*(delN+delm1))

    return d
end
