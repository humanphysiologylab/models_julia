function calc_if!(du, u, p, t, a)

    @unpack gIf = p
    @unpack EK, ENa = a
    V, y = u[24], u[14]
    
    if_y_inf = 1.0 / (1.0 + exp((V + 97.82874) / 12.48025))
    if_y_tau = 1.0 / (0.00332 * exp((-V) / 16.54103) + 23.71839 * exp(V / 16.54103))
    du[14] = d_y = (if_y_inf - y) / if_y_tau
    IfK = gIf * y * ((1.0 - 0.2677) * (V - EK))
    IfNa = gIf * y * (0.2677 * (V - ENa))
    If = IfK + IfNa
    
    @pack! a = if_y_inf, if_y_tau, IfK, IfNa, If
    nothing
end
