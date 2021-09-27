function calc_fluo!(du, u, p, t, a)

    @unpack fluo_tot, fluo_k_on, fluo_k_off = p

    for i in 1:5

        i_fluo = 43 + i
        i_Cai  = 4 + i

        fluo = u[i_fluo]
        Cai  = u[i_Cai]

        du[i_fluo] = fluo_k_on * Cai * (fluo_tot - fluo) - fluo_k_off * fluo
        du[i_Cai] -= du[i_fluo]

    end

end