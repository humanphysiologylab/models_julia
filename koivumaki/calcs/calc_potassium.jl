function calc_potassium!(du, u, p, t, a)

    @unpack G_seal, Vcytosol, F = p
    @unpack i_stim, IKur, INaK, IKs, IKr, IfK, IK1, It, IKb, IK_ghk = a
    
    
    i_tot = It + IKur + IK1 + IKr + IKs - 2.0 * INaK + IfK + i_stim
    i_tot += IKb + G_seal * IK_ghk
    du[25] = d_Ki = (-i_tot) / (Vcytosol * F)
    
    @pack! a = i_tot
    nothing
end
