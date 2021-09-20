function calc_sodium!(du, u, p, t, a)

    @unpack G_seal, Aj_nj, xj_nj_Nai, DNa, BNa, KdBNa, Vnonjunct_Nai, Vss, F = p
    @unpack INaCa, INaK, IfNa, INa, INab, INa_ghk = a
    Nai, Nass = u[42], u[43]
    
    i_ss = INa + 3.0 * INaK + 3.0 * INaCa + IfNa
    i_ss += INab + G_seal * INa_ghk
    JNa = DNa * Aj_nj / xj_nj_Nai * (Nass - Nai) * 1e-06
    betaNass = 1.0 / (1.0 + BNa * KdBNa / ((Nass + KdBNa) ^ 2.0))
    du[42] = d_Nai = JNa / Vnonjunct_Nai
    du[43] = d_Nass = betaNass * ((-JNa) / Vss - i_ss / (Vss * F))
    
    @pack! a = i_ss, JNa, betaNass
    nothing
end
