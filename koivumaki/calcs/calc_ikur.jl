function calc_ikur!(du, u, p, t, a)

    @unpack gKur = p
    @unpack EK = a
    V, ikur_r, ikur_s = u[24], u[17], u[18]
    
    ikur_r_inf = 1.0 / (1.0 + exp((V + 6.0) / (-8.6)))
    ikur_r_tau = 0.009 / (1.0 + exp((V + 5.0) / 12.0)) + 0.0005
    ikur_s_inf = 1.0 / (1.0 + exp((V + 7.5) / 10.0))
    ikur_s_tau = 0.59 / (1.0 + exp((V + 60.0) / 10.0)) + 3.05
    du[17] = d_ikur_r = (ikur_r_inf - ikur_r) / ikur_r_tau
    du[18] = d_ikur_s = (ikur_s_inf - ikur_s) / ikur_s_tau
    IKur = gKur * ikur_r * ikur_s * (V - EK)
    
    @pack! a = ikur_r_inf, ikur_r_tau, ikur_s_inf, ikur_s_tau, IKur
    nothing
end
