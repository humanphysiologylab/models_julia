function calc_iks!(du, u, p, t, a)

    @unpack gKs = p
    @unpack EK = a
    V, n = u[24], u[16]
    
    iks_n_inf = 1.0 / (1.0 + exp((V - 19.9) / (-12.7)))
    iks_n_tau = 0.4 * exp((-(((V - 20.0) / 20.0) ^ 2.0))) + 0.7
    du[16] = d_n = (iks_n_inf - n) / iks_n_tau
    IKs = gKs * n * (V - EK)
    
    @pack! a = iks_n_inf, iks_n_tau, IKs
    nothing
end
