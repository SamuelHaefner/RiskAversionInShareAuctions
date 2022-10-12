# Contains the script to read in the estimated monotonicity violations 
# obtained and saved with TestMonGeneric.jl and produces the respective tables.

using LaTeXStrings
using Latexify
using BSON: @save, @load  ## @save and @load

include("Auxiliary.jl")

rhovec = [exp(x) for x in [-10:0.5:0;]]

## loads and returns the data in the file TestMonA[auction]R[R]rho[rhodinex].dat
function TestData(auction, R, rhoindex)
    eval(Meta.parse(string(
        "@load \"TestMonA",
        auction,
        "R",
        R,
        "rho",
        rhoindex,
        ".dat\" T",
    )))
    return T
end

# Returns Monotonicity Violations for R=50,100,200
function MonTestRead()
    rhovec = [exp(x) for x in [-10:0.5:0;]]
    MonTest = []
    for R in [50,100,200]
        Ma = Matrix{Float64}(undef,39,length(rhovec))
        for auction in [1:1:39;]
            for rhoindex in [1:1:length(rhovec)-1;]
                T=TestData(auction, R, rhoindex)
                Ma[auction,rhoindex] = sum([T[1][x][2].>0 for x in [1:1:activebidders[auction];]])
            end
        end
        MonTest=push!(MonTest,Ma)
    end
    return MonTest
end

## Construct the data for the three tables (R=50,100,200) in Table 4 of the Main Text 
T = MonTestRead()
T50 = DataFrame(
    rho = rhovec[1:20],
    av = [mean(T[1][:, x] ./ (72 * 4.4)) for x in [1:1:20;]],
    se = [std(T[1][:, x] ./ (72 * 4.4)) for x in [1:1:20;]],
)
T100 = DataFrame(
    rho = rhovec[1:20],
    av = [mean(T[2][:, x] ./ (72 * 4.4)) for x in [1:1:20;]],
    se = [std(T[2][:, x] ./ (72 * 4.4)) for x in [1:1:20;]],
)
T200 = DataFrame(
    rho = rhovec[1:20],
    av = [mean(T[3][:, x] ./ (72 * 4.4)) for x in [1:1:20;]],
    se = [std(T[3][:, x] ./ (72 * 4.4)) for x in [1:1:20;]],
)


## Return the data for the three tables (R=50,100,200) in Table 4 of the Main Text 
t1=latexify(
    T50[:, [:rho, :av, :se]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
write("MonotonicityTest50.txt",t1)
t2=latexify(
    T100[:, [:rho, :av, :se]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
write("MonotonicityTest100.txt",t2)
t3=latexify(
    T200[:, [:rho, :av, :se]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
write("MonotonicityTest200.txt",t3)
