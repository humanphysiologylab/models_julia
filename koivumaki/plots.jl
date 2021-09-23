using PyPlot
using PyCall
pygui(true)

##

plot(sol.t, sol[:V], ".-")
# plot(sol[:V], sol[:m], ".-")
