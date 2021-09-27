function calc_icab!(du, u, p, t, a)

    @unpack gCab, Ca_ghk_scaler_membrane = p
    @unpack ECa, ICa_ghk = a
    V = u[24]
    
    # ICab = gCab * (V - ECa)
    ICab = Ca_ghk_scaler_membrane * ICa_ghk
    
    @pack! a = ICab
    nothing
end
