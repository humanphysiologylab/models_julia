function calc_nernst!(du, u, p, t, a)

    @unpack Cao, RTF, Ko, Nao = p
    
    Cass, Ki, Nass = u[9], u[25], u[43]
    
    ECa = RTF * log(Cao / Cass) / 2.0
    EK = RTF * log(Ko / Ki)
    ENa = RTF * log(Nao / Nass)
    
    @pack! a = ECa, EK, ENa
    nothing
end
