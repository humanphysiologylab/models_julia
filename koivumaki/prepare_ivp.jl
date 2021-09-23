using LabelledArrays: LVector

function create_NamedTuple_from_legend(legend)
    name = Symbol.(legend.name)
    value = Vector{Float64}(legend.value)
    result = (; zip(name, value)...)
end


function create_Dict_from_legend(legend)
    name = legend.name
    value = Vector{Float64}(legend.value)
    result = Dict(zip(name, value))
end


function create_LVector_from_legend(legend)
    result = LVector(create_NamedTuple_from_legend(legend))
end


function Dict_from_NamedTuple(nt)
    k = String.(keys(nt))
    v = values(nt)
    Dict(zip(k, v))
end


function prepare_ivp(legends)
    u = LVector(legends["states"])
    p = LVector(legends["constants"])
    a = Dict{String, Real}(Dict_from_NamedTuple(legends["algebraic"]))
    return (;u, p, a)
end


# function create_a(legend=legend_algebraic)
#     name = Symbol.(legend.name)
#     value = Vector{Real}(legend.value)
#     @LVector value Tuple(name)
# end