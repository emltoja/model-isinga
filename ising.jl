# Metropolis algorithm

# IMPORTY 
####################### 
using Statistics
using Plots
#######################

# STAŁE GLOBALNE
#######################
Coordinate = UInt8
const J  = 1
const kb = 1
const β  = 1/kb
#######################

"""
Zapisz siatkę do pliku 
"""
function _writelattice(lattice :: Matrix{Int8}, file :: IOStream)
    for row in eachrow(lattice)
        write(file, repr(row) * "\n")
    end
    write(file, "\n\n")
end




"""

Zwróć siatkę losowo ustawionych spinów (-1 lub 1) dla `setupmode = :random`
albo spinów równych 1 dla `setupmode = :ordered`

"""
function simsetup(size :: Coordinate, setupmode :: Symbol = :random) 
    
    if setupmode === :random 

        return rand([-Int8(1), Int8(1)], size, size)

    elseif setupmode === :ordered

        return ones(Int8, size, size)

    else error("setupmode should be either :random or :ordered")

    end

end




"""
Zmień spin jeśli spowoduje to spadek energii układu. Jeśli spowoduje wzrost to zmień
z prawdopodobieństwem `exp(-β * ΔE / T)`
"""
function metro!(lat :: Matrix{Int8}, i :: Coordinate, j :: Coordinate, temp :: Float64)

    ΔE = getenergydifference(lat, i, j)
    if ΔE <= 0
        return lat[i, j] *= -1
    end
    return rand() <= exp(-β * ΔE / temp) && (lat[i, j] *= -1) 

end




"""
Zwróć sumaryczne namagentyzowanie układu 
"""
function magnetization(lat :: Matrix{Int8})
    return mean(lat)
end




"""
Zwróć energię odzdziaływań z sąsiadami (siatka kwadratowa) dla danego agenta 
"""
function getenergy(lat :: Matrix{Int8}, i :: Coordinate, j :: Coordinate)
    
    topneighbor    = 0
    leftneighbor   = 0
    rightneighbor  = 0
    bottomneighbor = 0

    # Przypisz wartości do poszczegołnych sąsiadów (jeśli istnieją)
    try 
        topneighbor = lat[i, j - 1]
    catch
    end

    try 
        leftneighbor = lat[i - 1, j]
    catch
    end

    try 
        rightneighbor = lat[i + 1, j]
    catch
    end

    try
        bottomneighbor = lat[i, j + 1]
    catch
    end


    return - J * lat[i, j] * sum([topneighbor, leftneighbor, rightneighbor, bottomneighbor])
        
end




"""
Zwróć przyrost energii wynikający ze zmiany spinu danego agenta 
"""
function getenergydifference(lat :: Matrix{Int8}, i :: Coordinate, j :: Coordinate)
    return 2 * getenergy(lat, i, j)
end




"""
Zwróć sumaryczną energię całego układu 
"""
function getenergysystem(lat :: Matrix{Int8})
    totalenergy = 0
    for i = 1:size(lat, 1)
        for j = 1:size(lat, 2)
            totalenergy += getenergy(lat, i ,j)
        end
    end
    return totalenergy
end




"""
Przeprowadź symulację modelu Isinga 
"""
function simulation(
    size :: Coordinate;                 # Wymiary siatki
    runtime :: Int = 100_000,           # Długość symulacji 
    setupmode :: Symbol = :random,      # Stan początkowy
    animdest :: String = "./anim.gif",  # Ścieżka do animacji
    serialize :: Bool = true,           # True - zapisz dane do pliku
    dumpdest :: String = "./data.txt"   # Ścieżka do zapisywania danych
)

    lattice = simsetup(size, setupmode)

    if serialize
        touch(dumpdest)
        file = open(dumpdest, "w")
    end
    
    anim = @animate for MCS in 1:runtime 
        
        if !isnothing(file)
            _writelattice(lattice, file)
        end
        indexes = ceil.(Coordinate, size .* rand(size^2, 2))
        
        for point in eachrow(indexes)
            metro!(lattice, point[1], point[2], 1.0)
        end
        
        heatmap(lattice, legend=false, axis_ratio=1)
        MCS += 1
    end
    
    if !isnothing(file)
        close(file)
    end
    
    gif(anim, animdest, fps=25)

end

    
    
simulation(
    Coordinate(100),
    setupmode = :random,
    animdest = "./anim.gif",
    serialize = true,
    dumpdest = "./data.txt"
)