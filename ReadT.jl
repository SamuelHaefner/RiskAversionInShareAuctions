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
    eval(Meta.parse(string("@load \"Theta", runno, ".dat\" Theta", runno)))
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
savefig(p,"BRViolationsGamma.pdf")

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



#### Thetas by Bidder Assignment

TDataQG1 = DataFrame(
    [rhovec mean([
        [ViolatedQAll[y][x][1] / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataQG1, (:x1 => :rho))
rename!(TDataQG1, (:x2 => :meanest))
TDataAdd = DataFrame([
    [ViolatedQAll[y][x][1] / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataQG1 = hcat(TDataQG1, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataQG1, (:x1 => :std))
TDataQG1 = hcat(TDataQG1, TDataAdd)

TDataPG1 = DataFrame(
    [rhovec mean([
        [ViolatedPAll[y][x][1] / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataPG1, (:x1 => :rho))
rename!(TDataPG1, (:x2 => :meanest))
TDataAdd = DataFrame([
    [ViolatedPAll[y][x][1] / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataPG1 = hcat(TDataPG1, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataPG1, (:x1 => :std))
TDataPG1 = hcat(TDataPG1, TDataAdd)

TDataSumG1 = DataFrame(
    [rhovec mean([
        [(ViolatedPAll[y][x][1]+ViolatedQAll[y][x][1]) / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataSumG1, (:x1 => :rho))
rename!(TDataSumG1, (:x2 => :meanest))
TDataAdd = DataFrame([
    [(ViolatedPAll[y][x][1]+ViolatedQAll[y][x][1]) / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataSumG1 = hcat(TDataSumG1, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataSumG1, (:x1 => :std))
TDataSumG1 = hcat(TDataSumG1, TDataAdd)

TDataQG2 = DataFrame(
    [rhovec mean([
        [ViolatedQAll[y][x][2] / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataQG2, (:x1 => :rho))
rename!(TDataQG2, (:x2 => :meanest))
TDataAdd = DataFrame([
    [ViolatedQAll[y][x][2] / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataQG2 = hcat(TDataQG2, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataQG2, (:x1 => :std))
TDataQG2 = hcat(TDataQG2, TDataAdd)

TDataPG2 = DataFrame(
    [rhovec mean([
        [ViolatedPAll[y][x][2] / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataPG2, (:x1 => :rho))
rename!(TDataPG2, (:x2 => :meanest))
TDataAdd = DataFrame([
    [ViolatedPAll[y][x][2] / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataPG2 = hcat(TDataPG2, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataPG2, (:x1 => :std))
TDataPG2 = hcat(TDataPG2, TDataAdd)

TDataSumG2 = DataFrame(
    [rhovec mean([
        [(ViolatedPAll[y][x][2]+ViolatedQAll[y][x][2]) / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataSumG2, (:x1 => :rho))
rename!(TDataSumG2, (:x2 => :meanest))
TDataAdd = DataFrame([
    [(ViolatedPAll[y][x][2]+ViolatedQAll[y][x][2]) / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataSumG2 = hcat(TDataSumG2, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataSumG2, (:x1 => :std))
TDataSumG2 = hcat(TDataSumG2, TDataAdd)

TDataQG3 = DataFrame(
    [rhovec mean([
        [ViolatedQAll[y][x][3] / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataQG3, (:x1 => :rho))
rename!(TDataQG3, (:x2 => :meanest))
TDataAdd = DataFrame([
    [ViolatedQAll[y][x][3] / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataQG3 = hcat(TDataQG3, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataQG3, (:x1 => :std))
TDataQG3 = hcat(TDataQG3, TDataAdd)

TDataPG3 = DataFrame(
    [rhovec mean([
        [ViolatedPAll[y][x][3] / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataPG3, (:x1 => :rho))
rename!(TDataPG3, (:x2 => :meanest))
TDataAdd = DataFrame([
    [ViolatedPAll[y][x][3] / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataPG3 = hcat(TDataPG3, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataPG3, (:x1 => :std))
TDataPG3 = hcat(TDataPG3, TDataAdd)

TDataSumG3 = DataFrame(
    [rhovec mean([
        [(ViolatedPAll[y][x][3]+ViolatedQAll[y][x][3]) / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],:auto)
rename!(TDataSumG3, (:x1 => :rho))
rename!(TDataSumG3, (:x2 => :meanest))
TDataAdd = DataFrame([
    [(ViolatedPAll[y][x][3]+ViolatedQAll[y][x][3]) / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]],:auto)
TDataSumG3 = hcat(TDataSumG3, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataSumG3, (:x1 => :std))
TDataSumG3 = hcat(TDataSumG3, TDataAdd)

DataMeanest1=vcat(TDataQG1.meanest',TDataPG1.meanest')'
DataMeanest2=vcat(TDataQG2.meanest',TDataPG2.meanest')'
DataMeanest3=vcat(TDataQG3.meanest',TDataPG3.meanest')'
ctg = repeat([L"\Theta_{Q,g}", L"\Theta_{P,g}"], inner = 22)
nam = repeat(pushfirst!(log.(TDataP.rho[2:end]), -11),outer=2)
p1=groupedbar(
    nam,
    xlabel = L"\ln(\rho)",
    DataMeanest1,
    bar_position=:stack,
    group = ctg, 
    bar_width=0.5,
    title = "Bidder Group g=1",
    legend = :topleft
)
p2=groupedbar(
    nam,
    xlabel = L"\ln(\rho)",
    DataMeanest2,
    bar_position=:stack,
    group = ctg, 
    bar_width=0.5,
    title = "Bidder Group g=2",
    legend = :topleft
)
p3=groupedbar(
    nam,
    xlabel = L"\ln(\rho)",
    DataMeanest3,
    bar_position=:stack,
    group = ctg, 
    bar_width=0.5,
    title = "Bidder Group g=3",
    legend = :topleft
)
p=plot(p1,p2,p3,layout=(1,3),size = [1600, 400])
plot(p)
savefig(p,"BRViolationsBidderGroups.pdf")

data = hcat(
    TDataQG1[:, [:rho, :meanest]],
    TDataPG1[:, [:meanest]],
    TDataSumG1[:,[:meanest, :std]],
    TDataQG2[:, [:meanest]],
    TDataPG2[:, [:meanest]],
    TDataSumG2[:,[:meanest, :std]],
    TDataQG3[:, [:meanest]],
    TDataPG3[:, [:meanest]],
    TDataSumG3[:,[:meanest, :std]],
    makeunique=true
)

t=latexify(
    data,
    env=:tabular,
    fmt = x -> round(x, sigdigits = 3),
)
write("BRViolationsBidderGroups.txt",t)


## By Auction Group
# compute for group 1
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
            for g in 1
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
TDataQG1 = hcat(TDataQ, TDataAdd)

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
TDataPG1 = hcat(TDataP, TDataAdd)

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
TDataSumG1 = hcat(TDataSum, TDataAdd)

# compute for group 2
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
            for g in 2
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
TDataQG2 = hcat(TDataQ, TDataAdd)

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
TDataPG2 = hcat(TDataP, TDataAdd)

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
TDataSumG2 = hcat(TDataSum, TDataAdd)

# compute for group 3
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
            for g in 3
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
TDataQG3 = hcat(TDataQ, TDataAdd)

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
TDataPG3 = hcat(TDataP, TDataAdd)

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
TDataSumG3 = hcat(TDataSum, TDataAdd)

DataMeanest1=vcat(TDataQG1.meanest',TDataPG1.meanest')'
DataMeanest2=vcat(TDataQG2.meanest',TDataPG2.meanest')'
DataMeanest3=vcat(TDataQG3.meanest',TDataPG3.meanest')'
ctg = repeat([L"\Theta_{Q,a}", L"\Theta_{P,a}"], inner = 22)
nam = repeat(pushfirst!(log.(TDataP.rho[2:end]), -11),outer=2)
p1=groupedbar(
    nam,
    xlabel = L"\ln(\rho)",
    DataMeanest1,
    bar_position=:stack,
    group = ctg, 
    bar_width=0.5,
    title = "Auction Group a=1",
    legend = :topleft
)
p2=groupedbar(
    nam,
    xlabel = L"\ln(\rho)",
    DataMeanest2,
    bar_position=:stack,
    group = ctg, 
    bar_width=0.5,
    title = "Auction Group a=2",
    legend = :topleft
)
p3=groupedbar(
    nam,
    xlabel = L"\ln(\rho)",
    DataMeanest3,
    bar_position=:stack,
    group = ctg, 
    bar_width=0.5,
    title = "Bidder Group a=3",
    legend = :topleft
)
p=plot(p1,p2,p3,layout=(1,3),size = [1600, 400])
plot(p)
savefig(p,"BRViolationsAuctionGroups.pdf")

data = hcat(
    TDataQG1[:, [:rho, :meanest]],
    TDataPG1[:, [:meanest]],
    TDataSumG1[:,[:meanest, :std]],
    TDataQG2[:, [:meanest]],
    TDataPG2[:, [:meanest]],
    TDataSumG2[:,[:meanest, :std]],
    TDataQG3[:, [:meanest]],
    TDataPG3[:, [:meanest]],
    TDataSumG3[:,[:meanest, :std]],
    makeunique=true
)

t=latexify(
    data,
    env=:tabular,
    fmt = x -> round(x, sigdigits = 3),
)
write("BRViolationsAuctionGroups.txt",t)

## CE-Computation
CE(x, rho) = -1 / rho * log(0.5 * (1 + exp(-rho * x)))
