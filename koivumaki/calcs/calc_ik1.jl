function calc_ik1!(du, u, p, t, a)

    @unpack Ko, FRT, gK1 = p
    @unpack EK = a
    V = u[24]
    
    IK1 = gK1 * ((Ko * 1.0) ^ 0.4457) * (V - EK) / (1.0 + exp(1.5 * (V - EK + 3.6) * FRT))
    
    @pack! a = IK1
    nothing
end
