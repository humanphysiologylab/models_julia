function calc_icap!(du, u, p, t, a)

    @unpack ICaPmax, kCaP = p
    
    Cass = u[9]
    
    ICaP = ICaPmax * Cass / (kCaP + Cass)
    
    @pack! a = ICaP
    nothing
end
