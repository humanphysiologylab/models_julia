function calc_it!(du, u, p, t, a)

    @unpack gt = p
    @unpack EK = a
    V, it_r, it_s = u[24], u[22], u[23]
    
    it_r_inf = 1.0 / (1.0 + exp((V - 1.0) / (-11.0)))
    it_r_tau = 0.0035 * exp((-(((V + 0.0) / 30.0) ^ 2.0))) + 0.0015
    it_s_inf = 1.0 / (1.0 + exp((V + 40.5) / 11.5))
    it_s_tau = 0.025635 * exp((-(((V + 52.45) / 15.8827) ^ 2.0))) + 0.01414
    du[22] = d_it_r = (it_r_inf - it_r) / it_r_tau
    du[23] = d_it_s = (it_s_inf - it_s) / it_s_tau
    It = gt * it_r * it_s * (V - EK)
    
    @pack! a = it_r_inf, it_r_tau, it_s_inf, it_s_tau, It
    nothing
end
