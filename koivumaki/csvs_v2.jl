include("../read_csv.jl")

function read_legends(dirname_legends)
    filename_legend_constants = joinpath(dirname_legends, "legend_constants.csv")
    filename_legend_states = joinpath(dirname_legends, "legend_states.csv")
    filename_legend_algebraic = joinpath(dirname_legends, "legend_algebraic.csv")

    legend_states = read_csv(filename_legend_states)
    legend_constants = read_csv(filename_legend_constants)
    legend_algebraic = read_csv(filename_legend_algebraic)

    return Dict("states"    => legend_states,
                "constants" => legend_constants,
                "algebraic" => legend_algebraic)
end
