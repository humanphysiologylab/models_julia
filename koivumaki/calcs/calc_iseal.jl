function calc_iseal!(du, u, p, t, a)

    @unpack G_seal = p
    @unpack IK_ghk, INa_ghk, ICl_ghk, IAspartate_ghk, ICa_ghk = a
    V = u[24]
    
    Iseal = G_seal * (IK_ghk + INa_ghk + ICl_ghk + IAspartate_ghk + ICa_ghk)
    
    @pack! a = Iseal
    nothing
end
