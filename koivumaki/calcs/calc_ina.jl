function calc_ina!(du, u, p, t, a)

    @unpack Nao, F, FRT, PNa = p
    @unpack ENa = a
    V, h1, h2, m = u[24], u[19], u[20], u[21]
    
    h_inf = 1.0 / (1.0 + exp((V + 63.6) / 5.3))
    ina_h1_tau = 0.03 / (1.0 + exp((V + 35.1) / 3.2)) + 0.0003
    ina_h2_tau = 0.12 / (1.0 + exp((V + 35.1) / 3.2)) + 0.003
    ina_m_inf = 1.0 / (1.0 + exp((V + 27.12) / (-8.21)))
    ina_m_tau = 4.2e-05 * exp((-(((V + 25.57) / 28.8) ^ 2.0))) + 2.4e-05
    du[19] = d_h1 = (h_inf - h1) / ina_h1_tau
    du[20] = d_h2 = (h_inf - h2) / ina_h2_tau
    du[21] = d_m = (ina_m_inf - m) / ina_m_tau
    INa = PNa * ((m) ^ 3.0) * (0.9 * h1 + 0.1 * h2) * Nao * V * F * FRT * (exp((V - ENa) * FRT) - 1.0) / (exp(V * FRT) - 1.0)
    
    @pack! a = h_inf, ina_h1_tau, ina_h2_tau, ina_m_inf, ina_m_tau, INa
    nothing
end
