module Ising
    
    # IMPORTY 
    # ========================
    using Statistics
    using Plots
    using ProgressBars
    using Random
    using JLD2
    # ========================
    
    # STAŁE GLOBALNE
    # ========================
    const Coordinate = UInt8
    const J  = 1
    const kb = 1
    const β  = 1/kb
    # ========================
    
    export simulate, animation, Coordinate


    """
    Zwróć siatkę losowo ustawionych spinów (-1 lub 1) dla `setupmode = :random`
    albo spinów równych 1 dla `setupmode = :ordered`
    """
    function simsetup(size :: Coordinate, setupmode :: Symbol = :random) 
        
        if setupmode === :random 
             return rand([Int8(-1), Int8(1)], size, size)
        end

        if setupmode === :ordered
            return ones(Int8, size, size)
        end

        throw(ArgumentError("`setupmode` should be either :ordered or :random"))

    end

    """
    Zwróć sumaryczne namagentyzowanie układu 
    """
    function magnetization(lat :: Matrix{Int8})
        return mean(lat)
    end

    """
    Zwróć przyrost energii wynikający ze zmiany spinu danego agenta 
    """
    function getenergydifference(lat :: Matrix{Int8}, i :: Coordinate, j :: Coordinate)
        
        L = size(lat, 1)
        return 2 * J * lat[i, j] * (lat[i, mod1(j - 1, L)] + lat[mod1(i - 1, L), j] + lat[mod1(i + 1, L), j] + lat[i, mod1(j + 1, L)])
            
    end
    
    """
    Zmień spin jeśli spowoduje to spadek energii układu. Jeśli spowoduje wzrost to zmień
    z prawdopodobieństwem `exp(-β * ΔE / T)`
    """
    function metro!(lat :: Matrix{Int8}, i :: Coordinate, j :: Coordinate, temp :: Float64)

        ΔE = getenergydifference(lat, i, j)
        if ΔE <= 0
            return lat[i, j] *= -1
        return (temp * randexp() > β * ΔE) && (lat[i, j] *= -1)

    end

    """
    Symulacja modelu Isinga
    Zwraca wektor magentyzacji i zapisuje siatkę po każdym kroku (dla serialize == true)
    """
    function simulate(
        size :: Coordinate;                      # Wymiary siatki
        temp :: Float64 = 1.0,                   # Temperatura 
        runtime :: Int = 100,                    # Długość symulacji 
        setupmode :: Symbol = :random,           # Stan początkowy
        serialize :: Bool = true,                # True - zapisz dane do pliku
        dumpdest :: String = "./data/data.jld2", # Ścieżka do zapisywania danych
    )

        lattice = simsetup(size, setupmode)

        # Magnetyzacje po każdym kroku 
        mags = Vector{Float64}(undef, runtime)

        if serialize
            touch(dumpdest)
            file = jldopen(dumpdest, "w")
        end

        # Symulacja Monte-Carlo
        @simd for MCS in ProgressBar(1:runtime)
            
            mags[MCS] = magnetization(lattice)

            if serialize
                file["$MCS"] = lattice
            end

            # Wylosowanie L² indeksów do symulacji
            indexes = ceil.(Coordinate, size .* rand(Int64(size)^2, 2))
            
            # Wywołanie algorytmu Metropolisa dla wylosowanyhch indeksów
            @simd for point in eachrow(indexes)
                metro!(lattice, point[1], point[2], temp)
            end
            
        end

        if serialize
            close(file)
        end

        return mags 

    end

    """
    Zwróć animację ewolucji układu utworzoną z danych z pliku źródłowego.
    """
    function animation(src, colorscheme)
        snaps = load(src)

        anim = @animate for idx in ProgressBar(1:length(snaps))
            heatmap(snaps["$idx"], legend=nothing, c=colorscheme)
        end

        return anim
    end

end