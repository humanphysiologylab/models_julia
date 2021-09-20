using DataFrames: DataFrame
using CSV: File as CSVFile
using LabelledArrays

read_csv(filename) = DataFrame(CSVFile(filename))

dirname_data = "./"
dirname_legends = joinpath(dirname_data, "legends")

filename_legend_constants = joinpath(dirname_legends, "legend_constants.csv")
filename_legend_states = joinpath(dirname_legends, "legend_states.csv")
filename_legend_algebraic = joinpath(dirname_legends, "legend_algebraic.csv")


legend_states = read_csv(filename_legend_states)
legend_constants = read_csv(filename_legend_constants)
legend_algebraic = read_csv(filename_legend_algebraic)


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


# dirname_protocols = joinpath(dirname_data, "protocols")
# filename_protocol = joinpath(dirname_protocols, "protocol_sparse.csv")
# protocol = read_csv(filename_protocol)
# find_step(t, protocol=protocol) = protocol.v[findfirst(x -> x >= t, protocol.t)]