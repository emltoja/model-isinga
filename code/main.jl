using .Ising
using Plots
using JLD2

mags = simulate(Coordinate(100); temp=1.0)
scatter(mags[1:100:end], markersize=0.1)

anim = animation("./data/data.jld2", :imola10)

gif(anim, "./data/anim_imola.gif")


# Zadanie 1

# L = 10, T = 1
mags10 = simulate(
    UInt8(10);
    runtime=100_000,
    serialize=false,
    dumpdest="./data/test.jld2"
)

scatter(mags10)