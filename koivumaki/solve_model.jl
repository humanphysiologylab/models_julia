## This is used for python integration

using DifferentialEquations
using Sundials
using LabelledArrays

include("compute_rates_algebraic.jl")
include("csvs_v2.jl")
include("prepare_ivp.jl")
include("callbacks.jl")

include("calcs/calc_means.jl")

dirname_legends = "legends/"
legends = read_legends(dirname_legends)

solver = Rodas5()
dt = 1e-6
abstol = 1e-9
reltol = 1e-3
saveat = 1e-3
kw_solve = (;solver, dt, abstol, reltol, saveat)

function solve_model!(initial_state::Vector{Float64},
                      constants::Vector{Float64}, 
                      n_beats::Int, saveat::Float64;
                      # sol_output::Matrix{Float64},

                      legends::Dict=legends,
                      kw_solve::NamedTuple=kw_solve)

    legend_states    = legends["states"]
    legend_constants = legends["constants"]
    legend_algebraic = legends["algebraic"]

    @unpack u, p, a = prepare_ivp(legends)
    u .= initial_state
    p .= constants

    CaSR = u[end-2]  # desired
    u[end-2: end] = calc_means(u, p)
    rescale_CaSRs!(u, CaSR)

    tspan = (0., n_beats * p.STIM_PERIOD)

    rhs = ODEFunction((du, u, p, t) -> compute_rates_algebraic!(du, u, p, t, a),
                      syms=keys(u))
    cb_set = create_cb_set()                  
    prob = ODEProblem(rhs, u, tspan, p, callback=cb_set)

    sol = solve(prob; kw_solve...)

    output = hcat(sol.u...)

    m = map(u -> calc_means(u, p), sol)
    m = hcat(m...)
    output[end-2: end, :] = m

    output = Matrix(output)

    return output

end


##

# example: 
# output = solve_model!(Vector(u), Vector(p), 2, 1e-3)