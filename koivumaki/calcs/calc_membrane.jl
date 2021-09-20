function calc_membrane!(du, u, p, t, a)

    @unpack Cm = p
    @unpack ICaL, IKur, INa, It, INab, IKs, IKr, IK1, INaCa, ICab, INaK, ICaP, If, IKb, Iseal, i_stim = a
    
    i_ion = INa + ICaL + It + IKur + IK1 + IKr + IKs + INab + ICab + INaK + ICaP + INaCa + If + IKb + Iseal
    du[24] = d_V = (-(i_ion + i_stim)) / Cm
    
    @pack! a = i_ion
    nothing
end
