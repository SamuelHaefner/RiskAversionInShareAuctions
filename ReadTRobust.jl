## This file contains the script to read in estimates of T and produce the plots and tables.

include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

using LaTeXStrings
using Latexify
using Plots
using StatsPlots

# number of runs
runs = 40

# number of boostrap rounds per run
m = 5

# rhovector used
rhovec = sort!(append!([exp(x) for x in [-10:0.5:0;]], 0))

for runno in [1:1:runs;]
    eval(Meta.parse(string("@load \"ThetaRobust", runno, ".dat\" Theta", runno)))
end

Tested = []
ViolatedPAll = []
ViolatedQAll = []
for runno in [1:1:runs;]
    for bootstrapround in [1:1:m;]
        global TestedPQBootstrap = []
        global ViolatedQBootstrap = []
        global ViolatedPBootstrap = []
        for rho in [1:1:length(rhovec);]
            global TestedPQ = zeros(length(unique(bidderassignment)), 1)
            global ViolatedP = zeros(length(unique(bidderassignment)), 1)
            global ViolatedQ = zeros(length(unique(bidderassignment)), 1)
            for g in [1:1:length(group);]
                for auction in [1:1:length(group[g]);]
                    eval(Meta.parse(string(
                        "TestedPQ.+=Theta",
                        runno,
                        "[",
                        g,
                        "][",
                        auction,
                        "][",
                        rho,
                        "][1][:,",
                        bootstrapround,
                        "]",
                    )))
                    eval(Meta.parse(string(
                        "ViolatedQ.+=Theta",
                        runno,
                        "[",
                        g,
                        "][",
                        auction,
                        "][",
                        rho,
                        "][2][:,",
                        bootstrapround,
                        "]",
                    )))
                    eval(Meta.parse(string(
                        "ViolatedP.+=Theta",
                        runno,
                        "[",
                        g,
                        "][",
                        auction,
                        "][",
                        rho,
                        "][3][:,",
                        bootstrapround,
                        "]",
                    )))
                end
            end
            push!(TestedPQBootstrap, TestedPQ)
            push!(ViolatedQBootstrap, ViolatedQ)
            push!(ViolatedPBootstrap, ViolatedP)
        end
        push!(Tested, TestedPQBootstrap)
        push!(ViolatedQAll, ViolatedQBootstrap)
        push!(ViolatedPAll, ViolatedPBootstrap)
    end
end

TDataQ = DataFrame(
    [rhovec mean([
        [
            sum((ViolatedQAll[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
            for x in [1:1:length(rhovec);]
        ] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataQ, (:x1 => :rho))
rename!(TDataQ, (:x2 => :meanest))
TDataAdd = DataFrame([
    [
        sum((ViolatedQAll[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
        for x in [1:1:length(rhovec);]
    ] for y in [1:1:m*runs;]],:auto)
TDataQ = hcat(TDataQ, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataQ, (:x1 => :std))
TDataQ = hcat(TDataQ, TDataAdd)

TDataP = DataFrame(
    [rhovec mean([
        [
            sum((ViolatedPAll[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
            for x in [1:1:length(rhovec);]
        ] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataP, (:x1 => :rho))
rename!(TDataP, (:x2 => :meanest))
TDataAdd = DataFrame([
    [
        sum((ViolatedPAll[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
        for x in [1:1:length(rhovec);]
    ] for y in [1:1:m*runs;]],:auto)
TDataP = hcat(TDataP, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataP, (:x1 => :std))
TDataP = hcat(TDataP, TDataAdd)

TDataSum= DataFrame(
    [rhovec mean([
        [
            sum(((ViolatedPAll[y][x]+ViolatedQAll[y][x]) ./ Tested[y][x]) .* [105, 15, 3]) / 123
            for x in [1:1:length(rhovec);]
        ] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataSum, (:x1 => :rho))
rename!(TDataSum, (:x2 => :meanest))
TDataAdd = DataFrame([
    [
        sum(((ViolatedPAll[y][x]+ViolatedQAll[y][x]) ./ Tested[y][x]) .* [105, 15, 3]) / 123
        for x in [1:1:length(rhovec);]
    ] for y in [1:1:m*runs;]],:auto)
TDataSum = hcat(TDataSum, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataSum, (:x1 => :std))
TDataSum = hcat(TDataSum, TDataAdd)

DataMeanest=vcat(TDataQ.meanest',TDataP.meanest')'
ctg = repeat([L"\Theta_{Q}", L"\Theta_{P}"], inner = 22)
nam = repeat(pushfirst!(log.(TDataP.rho[2:end]), -11),outer=2)
p=groupedbar(
    nam,
    xlabel = L"\ln(\rho)",
    DataMeanest,
    bar_position=:stack,
    group = ctg, 
    bar_width=0.5,
    title="Gamma distribution",
    legend = :topleft
)
p=plot(p,size = [500, 400])
plot(p)
savefig(p,"BRViolationsLogNorm.pdf")

data = hcat(
    TDataQ[:, [:rho, :meanest]],
    TDataP[:, [:meanest]],
    TDataSum[:,[:meanest, :std]],
    makeunique=true
)

t=latexify(
    data,
    env=:tabular,
    fmt = x -> round(x, sigdigits = 3),
)
write("BRViolationsGamma.txt",t)



