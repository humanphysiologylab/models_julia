function calc_ical!(du, u, p, t, a)

    @unpack Cao, RTF, ECa_app, gCaL, kCa, kCan, ical_fca_tau = p
    V, d, f1, f2, Cass, fca = u[24], u[10], u[11], u[12], u[9], u[13]
    
    f_inf = 1.0 / (1.0 + exp((V + 27.4) / 7.1))
    ical_d_inf = 1.0 / (1.0 + exp((V + 9.0) / (-5.8)))
    ical_d_tau = 0.0027 * exp((-(((V + 35.0) / 30.0) ^ 2.0))) + 0.002
    ical_f1_tau = 0.98698 * exp((-(((V + 30.16047) / 7.09396) ^ 2.0))) + 0.04275 / (1.0 + exp((V - 51.61555) / (-80.61331))) + 0.03576 / (1.0 + exp((V + 29.57272) / 13.21758)) - 0.00821
    ical_f2_tau = 1.3323 * exp((-(((V + 40.0) / 14.2) ^ 2.0))) + 0.0626
    du[10] = d_d = (ical_d_inf - d) / ical_d_tau
    du[11] = d_f1 = (f_inf - f1) / ical_f1_tau
    du[12] = d_f2 = (f_inf - f2) / ical_f2_tau

    E_Ca = RTF / 2 * log(Cao / Cass)  # or ECa_app
    ICaL = gCaL * d * fca * f1 * f2 * (V - E_Ca)

    ical_fca_inf = 1.0 - 1.0 / (1.0 + (kCa / Cass) ^ kCan)
    du[13] = d_fca = (ical_fca_inf - fca) / ical_fca_tau
    
    @pack! a = f_inf, ical_d_inf, ical_d_tau, ical_f1_tau, ical_f2_tau, ICaL, ical_fca_inf
    nothing
end
