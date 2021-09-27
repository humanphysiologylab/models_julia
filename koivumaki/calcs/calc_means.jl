function calc_means(u, p)

    CaSR, VSR_sum = 0., 0.
    Cai, V_sum = 0., 0.
    fluo = 0.

    for i in 1:4

        VSRi = p[13 + i]
        CaSRi = u[0 + i]
        CaSR += CaSRi * VSRi
        VSR_sum += VSRi

        Vnonjuncti = p[7 + i]
        Caii = u[4 + i]
        fluo_i = u[43 + i]
        Cai += Caii * Vnonjuncti
        fluo += fluo_i * Vnonjuncti
        V_sum += Vnonjuncti

    end

    @unpack Cass, fluo_ss = u
    @unpack Vss = p

    Cai += Cass * Vss
    fluo += fluo_ss * Vss
    V_sum += Vss

    CaSR /= VSR_sum
    Cai  /= V_sum
    fluo /= V_sum

    return [CaSR, Cai, fluo]

end


function rescale_CaSRs!(u, CaSR_desired)
    
    @unpack CaSR = u
    
    for i in 1:5
        u[i] *= CaSR_desired / CaSR
    end

end
