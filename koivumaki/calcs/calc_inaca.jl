function calc_inaca!(du, u, p, t, a)

    @unpack Cao, Nao, FRT, dNaCa, fCaNCX, gam, kNaCa = p
    
    Cass, V, Nass = u[9], u[24], u[43]
    
    INaCa = kNaCa * (exp(gam * V * FRT) * ((Nass) ^ 3.0) * Cao - exp((gam - 1.0) * V * FRT) * ((Nao) ^ 3.0) * Cass * fCaNCX) / (1.0 + dNaCa * (((Nao) ^ 3.0) * Cass * fCaNCX + ((Nass) ^ 3.0) * Cao))
    
    @pack! a = INaCa
    nothing
end
