function calc_ryr!(du, u, p, t, a)

    @unpack tau_adapt, tau_inact, tau_inactss, tau_act, tau_actss = p
    @unpack nu1, nu2, nu3, nuss = p
    @unpack Jrel_multiplier = p
    
    CaSR1, CaSR2, CaSR3, CaSR4 = u[1], u[2], u[3], u[4]
    Cai1, Cai2, Cai3, Cass = u[5], u[6], u[7], u[9]
    ryr_a1, ryr_a2, ryr_a3, ryr_ass = u[26], u[27], u[28], u[29]
    c1, c2, c3, css = u[30], u[31], u[32], u[33]
    o1, o2, o3, oss = u[34], u[35], u[36], u[37]
    
    SRCa1 = 1.0 - 1.0 / (1.0 + exp((CaSR1 - 0.3) / 0.1))
    SRCa2 = 1.0 - 1.0 / (1.0 + exp((CaSR2 - 0.3) / 0.1))
    SRCa3 = 1.0 - 1.0 / (1.0 + exp((CaSR3 - 0.3) / 0.1))
    SRCass = 1.0 - 1.0 / (1.0 + exp((CaSR4 - 0.3) / 0.1))
    ainf1 = 0.505 - 0.427 / (1.0 + exp((Cai1 * 1000.0 - 0.29) / 0.082))
    ainf2 = 0.505 - 0.427 / (1.0 + exp((Cai2 * 1000.0 - 0.29) / 0.082))
    ainf3 = 0.505 - 0.427 / (1.0 + exp((Cai3 * 1000.0 - 0.29) / 0.082))
    ainfss = 0.505 - 0.427 / (1.0 + exp((Cass * 1000.0 - 0.29) / 0.082))
    cinf1 = 1.0 / (1.0 + exp((Cai1 * 1000.0 - (ryr_a1 + 0.02)) / 0.01))
    cinf2 = 1.0 / (1.0 + exp((Cai2 * 1000.0 - (ryr_a2 + 0.02)) / 0.01))
    cinf3 = 1.0 / (1.0 + exp((Cai3 * 1000.0 - (ryr_a3 + 0.02)) / 0.01))
    cinfss = 1.0 / (1.0 + exp((Cass * 1000.0 - (ryr_ass + 0.02)) / 0.01))
    oinf1 = 1.0 - 1.0 / (1.0 + exp((Cai1 * 1000.0 - (ryr_a1 + 0.22)) / 0.03))
    oinf2 = 1.0 - 1.0 / (1.0 + exp((Cai2 * 1000.0 - (ryr_a2 + 0.22)) / 0.03))
    oinf3 = 1.0 - 1.0 / (1.0 + exp((Cai3 * 1000.0 - (ryr_a3 + 0.22)) / 0.03))
    oinfss = 1.0 - 1.0 / (1.0 + exp((Cass * 1000.0 - (ryr_ass + 0.22)) / 0.03))
    du[26] = d_ryr_a1 = (ainf1 - ryr_a1) / tau_adapt
    du[27] = d_ryr_a2 = (ainf2 - ryr_a2) / tau_adapt
    du[28] = d_ryr_a3 = (ainf3 - ryr_a3) / tau_adapt
    du[29] = d_ryr_ass = (ainfss - ryr_ass) / tau_adapt
    du[30] = d_c1 = (cinf1 - c1) / tau_inact
    du[31] = d_c2 = (cinf2 - c2) / tau_inact
    du[32] = d_c3 = (cinf3 - c3) / tau_inact
    du[33] = d_css = (cinfss - css) / tau_inactss
    du[34] = d_o1 = (oinf1 - o1) / tau_act
    du[35] = d_o2 = (oinf2 - o2) / tau_act
    du[36] = d_o3 = (oinf3 - o3) / tau_act
    du[37] = d_oss = (oinfss - oss) / tau_actss
    Jrel1 = Jrel_multiplier * nu1 * o1 * c1 * SRCa1 * (CaSR1 - Cai1)
    Jrel2 = Jrel_multiplier * nu2 * o2 * c2 * SRCa2 * (CaSR2 - Cai2)
    Jrel3 = Jrel_multiplier * nu3 * o3 * c3 * SRCa3 * (CaSR3 - Cai3)
    Jrelss = Jrel_multiplier * nuss * oss * css * SRCass * (CaSR4 - Cass)
    
    @pack! a = SRCa1, SRCa2, SRCa3, SRCass
    @pack! a = ainf1, ainf2, ainf3, ainfss
    @pack! a = cinf1, cinf2, cinf3, cinfss
    @pack! a = oinf1, oinf2, oinf3, oinfss
    @pack! a = Jrel1, Jrel2, Jrel3, Jrelss
    nothing
end
