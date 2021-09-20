function calc_inak!(du, u, p, t, a)

    @unpack Ko, INaKmax, kNaKK, kNaKNa = p
    
    Nass, V = u[43], u[24]
    
    Nass15 = ((Nass * 1.0) ^ 1.5)
    INaK = INaKmax * Ko / (Ko + kNaKK) * Nass15 / (Nass15 + ((kNaKNa) ^ 1.5)) * (V + 150.0) / (V + 200.0)
    
    @pack! a = Nass15, INaK
    nothing
end
