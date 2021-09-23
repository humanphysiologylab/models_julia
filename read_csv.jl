sep = ","

filename = "koivumaki/legends/legend_states.csv"

function read_csv(filename::String, sep=",")::NamedTuple

    lines = readlines(filename)

    header = split(lines[1], sep)
    @assert all(header .== ["name", "value"])

    n = length(lines)

    columns = [Vector{String}(undef, n - 1),
               Vector{Float64}(undef, n - 1)]

    for (i, line) in enumerate(lines[2:end])
        x, y = split(line, sep)
        columns[1][i] = x
        columns[2][i] = parse(Float64, y)
    end

    nt = (; zip(Symbol.(columns[1]), columns[2])...)
end
