function calc_ikb!(du, u, p, t, a)

    @unpack K_ghk_scaler_membrane, gKb = p
    @unpack IK_ghk, EK = a
    V = u[24]
    
    IKb = K_ghk_scaler_membrane * IK_ghk
    # IKb = gKb * (V - EK)
    
    @pack! a = IKb
    nothing
end
