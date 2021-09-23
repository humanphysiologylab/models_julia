using Parameters

include("./calcs/calc_calcium.jl")
# include("./calcs/calc_fluo.jl")
# include("./calcs/calc_ghk_Aspartate.jl")
# include("./calcs/calc_ghk_Ca.jl")
# include("./calcs/calc_ghk_Cl.jl")
# include("./calcs/calc_ghk_K.jl")
# include("./calcs/calc_ghk_Na.jl")
include("./calcs/calc_icab.jl")
include("./calcs/calc_ical.jl")
include("./calcs/calc_icap.jl")
include("./calcs/calc_if.jl")
include("./calcs/calc_ik1.jl")
include("./calcs/calc_ikb.jl")
include("./calcs/calc_ikr.jl")
include("./calcs/calc_iks.jl")
include("./calcs/calc_ikur.jl")
include("./calcs/calc_ina.jl")
include("./calcs/calc_inab.jl")
include("./calcs/calc_inaca.jl")
include("./calcs/calc_inak.jl")
include("./calcs/calc_iseal.jl")
include("./calcs/calc_it.jl")
# include("./calcs/calc_means.jl")
include("./calcs/calc_membrane.jl")
include("./calcs/calc_nernst.jl")
include("./calcs/calc_potassium.jl")
include("./calcs/calc_ryr.jl")
include("./calcs/calc_serca.jl")
include("./calcs/calc_sodium.jl")
include("./calcs/calc_stimulus.jl")


function compute_rates_algebraic!(du, u, p, t, a)

    # calc_ghk_K!(du, u, p, t, a)
    # calc_ghk_Na!(du, u, p, t, a)
    # calc_ghk_Cl!(du, u, p, t, a)
    # calc_ghk_Aspartate!(du, u, p, t, a)
    # calc_ghk_Ca!(du, u, p, t, a)

    # calc_iseal!(du, u, p, t, a)

    du .= 0.

    map!(x -> 0., values(a))
    # a .= 0.

    # calc_stimulus!(du, u, p, t, a)  # if moved to callback
    calc_stimulus_v2!(du, u, p, t, a)

    calc_nernst!(du, u, p, t, a)

    calc_ikb!(du, u, p, t, a)
    calc_inab!(du, u, p, t, a)
    calc_icab!(du, u, p, t, a)

    calc_ical!(du, u, p, t, a)
    calc_icap!(du, u, p, t, a)
    calc_inaca!(du, u, p, t, a)
    calc_inak!(du, u, p, t, a)
    calc_ryr!(du, u, p, t, a)
    calc_serca!(du, u, p, t, a)
    calc_if!(du, u, p, t, a)
    calc_ik1!(du, u, p, t, a)
    calc_ikr!(du, u, p, t, a)
    calc_iks!(du, u, p, t, a)
    calc_ikur!(du, u, p, t, a)
    calc_ina!(du, u, p, t, a)
    calc_it!(du, u, p, t, a)

    calc_calcium!(du, u, p, t, a)
    calc_potassium!(du, u, p, t, a)
    calc_sodium!(du, u, p, t, a)

    calc_membrane!(du, u, p, t, a)

    # calc_fluo!(du, u, p, t, a)
    # calc_means!(du, u, p, t, a)

end
