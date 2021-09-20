function calc_calcium!(du, u, p, t, a)

    @unpack G_seal, Vnonjunct1, kSRleak, Vnonjunct2, Vnonjunct3, Vss, Aj_nj, xj_nj, DCa, CSQN, KdCSQN, BCa, KdBCa, KdSLhigh, KdSLlow, SLhigh, SLlow, dx, VSR1, DCaSR, VSR2, VSR3, VSR4, DCaBm, Vnonjunct4, F = p
    @unpack INaCa, ICaL, ICab, ICaP, ICa_ghk, J_bulkSERCA1, Jrel1, J_bulkSERCA2, Jrel2, J_bulkSERCA3, Jrel3, J_bulkSERCAss, Jrelss, J_SERCASR1, J_SERCASR2, J_SERCASR3, J_SERCASRss = a
    CaSR1, Cai1, CaSR2, Cai2, CaSR3, Cai3, CaSR4, Cass, Cai4 = u[1], u[5], u[2], u[6], u[3], u[7], u[4], u[9], u[8]
    
    calcium_Cass_i_tot = (-ICaL) - ICaP + 2.0 * INaCa - (ICab + G_seal * ICa_ghk)
    JSRCaleak1 = kSRleak * (CaSR1 - Cai1) * Vnonjunct1
    JSRCaleak2 = kSRleak * (CaSR2 - Cai2) * Vnonjunct2
    JSRCaleak3 = kSRleak * (CaSR3 - Cai3) * Vnonjunct3
    JSRCaleakss = kSRleak * (CaSR4 - Cass) * Vss
    Jj_nj = DCa * Aj_nj / xj_nj * (Cass - Cai4) * 1e-06
    calcium_CaSR1_beta = 1.0 / (1.0 + CSQN * KdCSQN / ((CaSR1 + KdCSQN) ^ 2.0))
    calcium_CaSR2_beta = 1.0 / (1.0 + CSQN * KdCSQN / ((CaSR2 + KdCSQN) ^ 2.0))
    calcium_CaSR3_beta = 1.0 / (1.0 + CSQN * KdCSQN / ((CaSR3 + KdCSQN) ^ 2.0))
    calcium_CaSR4_beta = 1.0 / (1.0 + CSQN * KdCSQN / ((CaSR4 + KdCSQN) ^ 2.0))
    calcium_Cai1_gamma = BCa * KdBCa / ((Cai1 + KdBCa) ^ 2.0)
    calcium_Cai2_gamma = BCa * KdBCa / ((Cai2 + KdBCa) ^ 2.0)
    calcium_Cai3_gamma = BCa * KdBCa / ((Cai3 + KdBCa) ^ 2.0)
    calcium_Cai4_gamma = BCa * KdBCa / ((Cai4 + KdBCa) ^ 2.0)
    calcium_Cass_beta = 1.0 / (1.0 + SLlow * KdSLlow / ((Cass + KdSLlow) ^ 2.0) + SLhigh * KdSLhigh / ((Cass + KdSLhigh) ^ 2.0) + BCa * KdBCa / ((Cass + KdBCa) ^ 2.0))
    JCa1 = (-J_bulkSERCA1) + JSRCaleak1 + Jrel1
    JCa2 = (-J_bulkSERCA2) + JSRCaleak2 + Jrel2
    JCa3 = (-J_bulkSERCA3) + JSRCaleak3 + Jrel3
    JCa4 = Jj_nj
    JCass = (-Jj_nj) + JSRCaleakss - J_bulkSERCAss + Jrelss
    JSRCa1 = J_SERCASR1 - JSRCaleak1 - Jrel1
    JSRCa2 = J_SERCASR2 - JSRCaleak2 - Jrel2
    JSRCa3 = J_SERCASR3 - JSRCaleak3 - Jrel3
    JSRCa4 = J_SERCASRss - JSRCaleakss - Jrelss
    calcium_Cai1_beta = 1.0 / (1.0 + calcium_Cai1_gamma)
    calcium_Cai2_beta = 1.0 / (1.0 + calcium_Cai2_gamma)
    calcium_Cai3_beta = 1.0 / (1.0 + calcium_Cai3_gamma)
    calcium_Cai4_beta = 1.0 / (1.0 + calcium_Cai4_gamma)
    du[1] = d_CaSR1 = calcium_CaSR1_beta * DCaSR * ((CaSR2 - 2.0 * CaSR1 + CaSR1) / ((dx) ^ 2.0) + (CaSR2 - CaSR1) / (2.0 * 1.0 * ((dx) ^ 2.0))) + JSRCa1 / VSR1 * calcium_CaSR1_beta
    du[2] = d_CaSR2 = calcium_CaSR2_beta * DCaSR * ((CaSR3 - 2.0 * CaSR2 + CaSR1) / ((dx) ^ 2.0) + (CaSR3 - CaSR1) / (2.0 * 2.0 * ((dx) ^ 2.0))) + JSRCa2 / VSR2 * calcium_CaSR2_beta
    du[3] = d_CaSR3 = calcium_CaSR3_beta * DCaSR * ((CaSR4 - 2.0 * CaSR3 + CaSR2) / ((dx) ^ 2.0) + (CaSR4 - CaSR2) / (2.0 * 3.0 * ((dx) ^ 2.0))) + JSRCa3 / VSR3 * calcium_CaSR3_beta
    du[4] = d_CaSR4 = calcium_CaSR4_beta * DCaSR * ((CaSR4 - 2.0 * CaSR4 + CaSR3) / ((dx) ^ 2.0) + (CaSR4 - CaSR3) / (2.0 * 4.0 * ((dx) ^ 2.0))) + JSRCa4 / VSR4 * calcium_CaSR4_beta
    du[5] = d_Cai1 = calcium_Cai1_beta * (DCa + calcium_Cai1_gamma * DCaBm) * ((Cai2 - 2.0 * Cai1 + Cai1) / ((dx) ^ 2.0) + (Cai2 - Cai1) / (2.0 * 1.0 * ((dx) ^ 2.0))) - 2.0 * calcium_Cai1_beta * calcium_Cai1_gamma * DCaBm / (KdBCa + Cai1) * (((Cai2 - Cai1) / (2.0 * dx)) ^ 2.0) + JCa1 / Vnonjunct1 * calcium_Cai1_beta
    du[6] = d_Cai2 = calcium_Cai2_beta * (DCa + calcium_Cai2_gamma * DCaBm) * ((Cai3 - 2.0 * Cai2 + Cai1) / ((dx) ^ 2.0) + (Cai3 - Cai1) / (2.0 * 2.0 * ((dx) ^ 2.0))) - 2.0 * calcium_Cai2_beta * calcium_Cai2_gamma * DCaBm / (KdBCa + Cai2) * (((Cai3 - Cai1) / (2.0 * dx)) ^ 2.0) + JCa2 / Vnonjunct2 * calcium_Cai2_beta
    du[7] = d_Cai3 = calcium_Cai3_beta * (DCa + calcium_Cai3_gamma * DCaBm) * ((Cai4 - 2.0 * Cai3 + Cai2) / ((dx) ^ 2.0) + (Cai4 - Cai2) / (2.0 * 3.0 * ((dx) ^ 2.0))) - 2.0 * calcium_Cai3_beta * calcium_Cai3_gamma * DCaBm / (KdBCa + Cai3) * (((Cai4 - Cai2) / (2.0 * dx)) ^ 2.0) + JCa3 / Vnonjunct3 * calcium_Cai3_beta
    du[8] = d_Cai4 = calcium_Cai4_beta * (DCa + calcium_Cai4_gamma * DCaBm) * ((Cai4 - 2.0 * Cai4 + Cai3) / ((dx) ^ 2.0) + (Cai4 - Cai3) / (2.0 * 4.0 * ((dx) ^ 2.0))) - 2.0 * calcium_Cai4_beta * calcium_Cai4_gamma * DCaBm / (KdBCa + Cai4) * (((Cai4 - Cai3) / (2.0 * dx)) ^ 2.0) + JCa4 / Vnonjunct4 * calcium_Cai4_beta
    du[9] = d_Cass = calcium_Cass_beta * (JCass / Vss + calcium_Cass_i_tot / (2.0 * Vss * F))
    
    @pack! a = calcium_Cass_i_tot, JSRCaleak1, JSRCaleak2, JSRCaleak3, JSRCaleakss, Jj_nj, calcium_CaSR1_beta, calcium_CaSR2_beta, calcium_CaSR3_beta, calcium_CaSR4_beta, calcium_Cai1_gamma, calcium_Cai2_gamma, calcium_Cai3_gamma, calcium_Cai4_gamma, calcium_Cass_beta, JCa1, JCa2, JCa3, JCa4, JCass, JSRCa1, JSRCa2, JSRCa3, JSRCa4, calcium_Cai1_beta, calcium_Cai2_beta, calcium_Cai3_beta, calcium_Cai4_beta
    nothing
end
