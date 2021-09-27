function calc_ghk(z, mobility, Ci, Co, u, p)

    @unpack F, R, T = p
    @unpack V = u

    zFVRT = z * V * F / R / T

    scaler = 1  # gamma / h
    P = mobility * R * T * scaler

    J = P * zFVRT * (Co - Ci * exp(zFVRT)) / (1 - exp(zFVRT))
    I = z * F * J
end


function calc_ghk_all!(du, u, p, t, a)

    z_K, z_Na, z_Cl, z_Aspartate, z_Ca = +1, +1, -1, -1, +2
    
    @unpack K_mobility, Na_mobility, Cl_mobility, Aspartate_mobility, Ca_mobility = p
    
    @unpack Ko, Nao, Cao, Cl_o, Aspartate_o = p
    
    @unpack Ki, Nai, Cass = u 
    @unpack Cl_i, Aspartate_i = p

    IK_ghk = calc_ghk(z_K, K_mobility, Ki, Ko, u, p)
    INa_ghk = calc_ghk(z_Na, Na_mobility, Nai, Nao, u, p)
    ICa_ghk = calc_ghk(z_Ca, Ca_mobility, Cass, Cao, u, p)
    ICl_ghk = calc_ghk(z_Cl, Cl_mobility, Cl_i, Cl_o, u, p)
    IAspartate_ghk = calc_ghk(z_Aspartate, Aspartate_mobility, Aspartate_i, Aspartate_o, u, p)

    @pack! a = IK_ghk, INa_ghk, ICa_ghk, ICl_ghk, IAspartate_ghk
    
end
