function calc_ikr!(du, u, p, t, a)

    @unpack gKr = p
    @unpack EK = a
    V, pa = u[24], u[15]
    
    ikr_pa_inf = 1.0 / (1.0 + exp((V + 15.0) / (-6.0)))
    ikr_pa_tau = 0.21718 * exp((-(((V + 20.1376) / 22.1996) ^ 2.0))) + 0.03118
    ikr_pi = 1.0 / (1.0 + exp((V + 55.0) / 24.0))
    du[15] = d_pa = (ikr_pa_inf - pa) / ikr_pa_tau
    IKr = gKr * pa * ikr_pi * (V - EK)
    
    @pack! a = ikr_pa_inf, ikr_pa_tau, ikr_pi, IKr
    nothing
end
