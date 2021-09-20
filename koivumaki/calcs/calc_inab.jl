function calc_inab!(du, u, p, t, a)

    @unpack gNab, Na_ghk_scaler_membrane = p
    @unpack ENa, INa_ghk = a
    V = u[24]
    
    INab = gNab * (V - ENa)
    # INab = Na_ghk_scaler_membrane * INa_ghk
    
    @pack! a = INab
    nothing
end
