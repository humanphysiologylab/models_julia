cd("/home/andrey/WORK/HPL/Code/ga/models_jl/koivumaki")

##

include("csvs_v2.jl")
include("compute_rates_algebraic.jl")
include("callbacks.jl")
include("prepare_ivp.jl")

##

using Parameters: @unpack

legends = read_legends("./legends/")
@unpack u, p, a = prepare_ivp(legends)

du = similar(u)
t = 0.

compute_rates_algebraic!(du, u, p, t, a)

##

using DifferentialEquations
using Sundials
using BenchmarkTools

##

u₀ = copy(u)

CL = 0.344
p.STIM_PERIOD = CL
n_beats = 10.
tspan = (0., n_beats * CL)

cb_set = create_cb_set()

##

rhs = ODEFunction((du, u, p, t) -> compute_rates_algebraic!(du, u, p, t, a),
                  syms=keys(u))
prob = ODEProblem(rhs, u₀, tspan, p, callback=cb_set)

##

using ModelingToolkit
jac = eval(generate_jacobian(modelingtoolkitize(prob)))
prob = ODEProblem(rhs, u₀, tspan, p, callback=cb_set, jac=jac)

##

solver = Rodas5()

saveat = 1e-3
dt = 1e-6

abstol = 1e-9
reltol = 1e-3

kw_solve = (;solver, saveat, dt, abstol, reltol)
# kw_solve = (;solver, dt)

sol = solve(prob; kw_solve...)

##
