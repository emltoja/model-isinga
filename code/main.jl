# IMPORTY 
# ===================
using .Ising
using Plots
using JLD2
using LaTeXStrings
using Statistics
# ===================


# FUNKCJE DO WYKONANIA WYKRESÓW DO RAPORTU 
# ===================================================================

function zad1(size, temp)

    mags = simulate(
        Coordinate(size);
        temp=temp,
        runtime=100,
        dumpdest="./data/$(size)x$(size)_10MCS_T$(temp).jld2"
        )

        p = plot(
            mags,
            markershape=:circle,
            legend=nothing,
            xlabel="Kroki MCS",
            ylabel=L"m"    
        )

        savefig(p, "./data/$(size)x$(size)_10MCS_T$(temp).png")

        return p
end

function snapshot(size, temp, colorscheme, moment)
   
    snap = load("./data/$(size)x$(size)_10MCS_T$(temp).jld2", "$moment")

    heatmap(snap, legend=nothing, c=colorscheme, axis=nothing, aspect_ratio=1, size=(400, 400))

end

function get_snaps(size, temp, colorscheme)
    for moment in (1, 50, 100)
        p = snapshot(size, temp, colorscheme, moment)
        savefig(p, "./data/snaps/snap$(size)_temp$(temp)_moment$(moment).png")
    end
end

function trajectories(size, temp, runs)

    p = plot(
        legend=nothing,
        xlabel="Kroki MCS",
        ylabel=L"m",
        ylims = (-1.1, 1.1)
    )

    for _ in 1:10
        
        mags = simulate(
            Coordinate(size);
            temp=temp,
            runtime=runs,
            serialize=false
            )

        plot!(mags)

    end

        savefig(p, "./data/magnetizations/$(size).png")

        return p
        
end
    
function mean_mag(size, runs)

    ts = 0.5:0.2:3.5
    result = Vector{Float64}(undef, length(ts))
    for (i, t) in enumerate(ts)
        mean_vec = Vector{Float64}(undef, 20)
        for j in 1:20
            mags = simulate(
                Coordinate(size);
                temp=t,
                runtime = runs,
                setupmode=:ordered,
                serialize=false
            )
            mean_vec[j] = abs(mags[end])
        end
        result[i] = mean(mean_vec)
    end
    return result
end
    
function mean_mag_time(size, runs, sample_interval)

    ts = 0.5:0.1:3.5
    result = Vector{Float64}(undef, length(ts))
    for (i, t) in enumerate(ts)

        mags = simulate(
            Coordinate(size);
            temp=t,
            runtime=runs,
            setupmode=:ordered,
            serialize=false
        )

        result[i] = mean(abs.(mags[end-sample_interval:end]))
    end

    return result
end

# ===================================================================



# Zadanie 1
zad1(10, 1.0)   # L = 10, T = 1
zad1(80, 1.0)   # L = 80, T = 1
zad1(10, 2.26)  # L = 10, T = 2.26
zad1(80, 2.26)  # L = 80, T = 2.26
zad1(10, 4.0)   # L = 10, T = 4
zad1(80, 4.0)   # L = 80, T = 4

# Zadanie 2
snapshot(80, 1.0, :imola, 50)
get_snaps(10, 1.0, :imola)
get_snaps(10, 2.26, :imola)
get_snaps(10, 4.0, :imola)
get_snaps(80, 1.0, :imola)
get_snaps(80, 2.26, :imola)
get_snaps(80, 4.0, :imola)

trajectories(10, 1.0, 250)
trajectories(20, 1.0, 1000)
trajectories(40, 1.0, 5000)
trajectories(80, 1.0, 10000)

# Zadanie 3
for t in (1.3, 2.26, 3.6)
    trajectories(10, t, 250)
    trajectories(20, t, 500)
    trajectories(40, t, 1000)
    trajectories(80, t, 5000)
end

# Zad 4 
mags_10 = mean_mag(10, 100)
mags_20 = mean_mag(20, 500)
mags_40 = mean_mag(40, 1000)
mags_80 = mean_mag(80, 5000)

mag_temp_p = plot(
    0.5:0.2:3.5,
    [mags_10 mags_20 mags_40 mags_80],
    markershape=[:utriangle :dtriangle :diamond :circle],
    label=["L=10" "L=20" "L=40" "L=80"],
    xlabel=L"T^*",
    ylabel=L"<m>",
    title="Magnetyzacja uśredniana po zespole"
)

savefig(mag_temp_p, "./data/magnetizations/mag_temp.png")

mags_10_time = mean_mag_time(10, 3000, 2000)
mags_20_time = mean_mag_time(20, 2000, 500)
mags_30_time = mean_mag_time(40, 3000, 500)
mags_80_time = mean_mag_time(80, 5100, 100)

mag_time_p = plot(
    0.5:0.1:3.5,
    [mags_10_time mags_20_time mags_30_time mags_80_time],
    markershape=[:utriangle :dtriangle :diamond :circle],
    label=["L=10" "L=20" "L=40" "L=80"],
    xlabel=L"T^*",
    ylabel=L"<m>",
    title="Magnetyzacja uśredniana po czasie"
)

savefig(mag_time_p, "./data/magnetizations/mag_time.png")