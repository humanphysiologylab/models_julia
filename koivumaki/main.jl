cd("koivumaki/")

include("csvs.jl")
include("compute_rates_algebraic.jl")

##

u = deepcopy(legend_states.value)
du = similar(u)

p = create_LVector_from_legend(legend_constants)
a = create_LVector_from_legend(legend_algebraic)

t = 0.
compute_rates_algebraic!(du, u, p, t, a)

##

using DifferentialEquations
using Sundials
using BenchmarkTools

u₀ = copy(u)

CL = 1.
p.STIM_PERIOD = CL
n_beats = 1
tspan = (0., n_beats * CL)

function reset_dt!(integrator)
    set_proposed_dt!(integrator, 1e-6)
end

function stimulate!(integrator)
    @unpack STIM_PERIOD, STIM_DURATION = integrator.p
    t = integrator.t
    if (t % STIM_PERIOD ≤ STIM_DURATION)
        integrator.p[83] = 1.
        # integrator.opts.dtmax = 5e-5
        # set_proposed_dt!(integrator, 1e-6)
    else
        integrator.p[83] = 0.
        # integrator.opts.dtmax = Inf
    end
end

cb = PresetTimeCallback(0: CL: tspan[2], reset_dt!)
# cb = DiscreteCallback((u, t, integrator) -> true, stimulate!, save_positions=(false, false))

##


##
# a = create_LVector_from_legend(legend_algebraic)
rhs = ODEFunction((du, u, p, t) -> compute_rates_algebraic!(du, u, p, t, a),
                  syms=keys(u))
prob = ODEProblem(rhs, u₀, tspan, p)


##

kw_solve = (solver = Rodas5(),
            saveat = 0: 1e-3: tspan[2],
            tstops = 0: CL: tspan[end])

sol = solve(prob, kw_solve...)