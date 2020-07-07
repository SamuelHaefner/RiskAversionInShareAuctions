## This file contains the script to read in estimates of T and produce the plots and tables.

include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

using LaTeXStrings
using Latexify
using Plots

# number of runs
runs = 40

# number of boostrap rounds per runf
m = 5

# rhovector used
rhovec = sort!(append!([exp(x) for x in [-10:0.5:0;]], 0))

for runno in [1:1:runs;]
    eval(Meta.parse(string("@load \"Theta", runno, ".dat\" Theta", runno)))
end

Tested = []
Violated = []
for runno in [1:1:runs;]
    for bootstrapround in [1:1:m;]
        global TestedRhoBootstrap = []
        global ViolatedRhoBootstrap = []
        for rho in [1:1:length(rhovec);]
            global TestedRho = zeros(length(unique(bidderassignment)), 1)
            global ViolatedRho = zeros(length(unique(bidderassignment)), 1)
            for g in [1:1:length(group);]
                for auction in [1:1:length(group[g]);]
                    eval(Meta.parse(string(
                        "TestedRho.+=Theta",
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
                        "ViolatedRho.+=Theta",
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
                end
            end
            push!(TestedRhoBootstrap, TestedRho)
            push!(ViolatedRhoBootstrap, ViolatedRho)
        end
        push!(Tested, TestedRhoBootstrap)
        push!(Violated, ViolatedRhoBootstrap)
    end
end

TData = DataFrame(
    [rhovec mean([
        [
            sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
            for x in [1:1:length(rhovec);]
        ] for y in [1:1:m*runs;]
    ])],
)
rename!(TData, (:x1 => :rho))
rename!(TData, (:x2 => :meanest))
TDataAdd = DataFrame([
    [
        sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
        for x in [1:1:length(rhovec);]
    ] for y in [1:1:m*runs;]
])
TData = hcat(TData, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TData, (:x1 => :std))
TData = hcat(TData, TDataAdd)

latexify(
    TData[:, [:rho, :meanest, :std]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)

default(lab = "")
ticks = pushfirst!(log.(TData.rho[2:end]), -11)
ticklabels = string.(pushfirst!(log.(TData.rho[2:end]), -Inf))
toplot = plot(
    pushfirst!(log.(TData.rho[2:end]), -11),
    TData.meanest,
    yerror=1.96*TData.std,
    ylim = (0.1, 0.8),
    seriestype = :scatter,
    #xticks = (ticks, ticklabels),
    grid=true,
    label = "Fraction of Best-Response Violations",
    legend = :topleft,
    xlabel = L"\ln(\rho)"
)
plot(toplot)
savefig("ThetaEst.pdf")

#### Thetas by Bidder Assignment

TDataG1 = DataFrame(
    [rhovec mean([
        [Violated[y][x][1] / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],
)
rename!(TDataG1, (:x1 => :rho))
rename!(TDataG1, (:x2 => :meanest))
TDataAdd = DataFrame([
    [Violated[y][x][1] / Tested[y][x][1] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
])
TDataG1 = hcat(TDataG1, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataG1, (:x1 => :std))
TDataG1 = hcat(TDataG1, TDataAdd)

TDataG2 = DataFrame(
    [rhovec mean([
        [Violated[y][x][2] / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],
)
rename!(TDataG2, (:x1 => :rho))
rename!(TDataG2, (:x2 => :meanest))
TDataAdd = DataFrame([
    [Violated[y][x][2] / Tested[y][x][2] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
])
TDataG2 = hcat(TDataG2, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataG2, (:x1 => :std))
TDataG2 = hcat(TDataG2, TDataAdd)

TDataG3 = DataFrame(
    [rhovec mean([
        [Violated[y][x][3] / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
    ])],
)
rename!(TDataG3, (:x1 => :rho))
rename!(TDataG3, (:x2 => :meanest))
TDataAdd = DataFrame([
    [Violated[y][x][3] / Tested[y][x][3] for x in [1:1:length(rhovec);]] for y in [1:1:m*runs;]
])
TDataG3 = hcat(TDataG3, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TDataG3, (:x1 => :std))
TDataG3 = hcat(TDataG3, TDataAdd)

default(lab = "")
ticks = pushfirst!(log.(TDataG1.rho[2:end]), -11)
ticklabels = string.(pushfirst!(log.(TDataG1.rho[2:end]), -Inf))
toplot = plot(
    pushfirst!(log.(TDataG1.rho[2:end]) .- 0.12, -11.12),
    TDataG1.meanest,
    yerror=1.96*TDataG1.std,
    ylim = (0.1, 0.8),
    seriestype = :scatter,
    #xticks = (ticks, ticklabels),
    label = "Fraction of BR Violations Bidder Group 1",
    legend = :topleft,
)
plot!(
    pushfirst!(log.(TDataG2.rho[2:end]), -11),
    TDataG2.meanest,
    ylim = (0.1, 0.8),
    yerror=1.96*TDataG2.std,
    seriestype = :scatter,
    #xticks = (ticks, ticklabels),
    label = "Fraction of BR Violations Bidder Group 2",
    legend = :topleft,
)
plot!(
    pushfirst!(log.(TDataG3.rho[2:end]) .+ 0.12, -11 + 0.12),
    TDataG3.meanest,
    yerror=1.96*TDataG3.std,
    ylim = (0.1, 0.8),
    seriestype = :scatter,
    #xticks = (ticks, ticklabels),
    label = "Fraction of BR Violations Bidder Group 3",
    legend = :topleft,
)

plot!(xlabel = L"\ln(\rho)")
plot(toplot)
savefig("ThetaEstBidderGroups.pdf")

latexify(
    TDataG1[:, [:rho, :meanest, :std]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
latexify(
    TDataG2[:, [:rho, :meanest, :std]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
latexify(
    TDataG3[:, [:rho, :meanest, :std]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)

## By Auction Group
Tested = []
Violated = []
for runno in [1:1:runs;]
    for bootstrapround in [1:1:m;]
        global TestedRhoBootstrap = []
        global ViolatedRhoBootstrap = []
        for rho in [1:1:length(rhovec);]
            global TestedRho = zeros(length(unique(bidderassignment)), 1)
            global ViolatedRho = zeros(length(unique(bidderassignment)), 1)
            for g in 1
                for auction in [1:1:length(group[g]);]
                    eval(Meta.parse(string(
                        "TestedRho.+=Theta",
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
                        "ViolatedRho.+=Theta",
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
                end
            end
            push!(TestedRhoBootstrap, TestedRho)
            push!(ViolatedRhoBootstrap, ViolatedRho)
        end
        push!(Tested, TestedRhoBootstrap)
        push!(Violated, ViolatedRhoBootstrap)
    end
end

TData = DataFrame(
    [rhovec mean([
        [
            sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
            for x in [1:1:length(rhovec);]
        ] for y in [1:1:m*runs;]
    ])],
)
rename!(TData, (:x1 => :rho))
rename!(TData, (:x2 => :meanest))
TDataAdd = DataFrame([
    [
        sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
        for x in [1:1:length(rhovec);]
    ] for y in [1:1:m*runs;]
])
TData = hcat(TData, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TData, (:x1 => :std))
TDataG1 = hcat(TData, TDataAdd)

Tested = []
Violated = []
for runno in [1:1:runs;]
    for bootstrapround in [1:1:m;]
        global TestedRhoBootstrap = []
        global ViolatedRhoBootstrap = []
        for rho in [1:1:length(rhovec);]
            global TestedRho = zeros(length(unique(bidderassignment)), 1)
            global ViolatedRho = zeros(length(unique(bidderassignment)), 1)
            for g in 2
                for auction in [1:1:length(group[g]);]
                    eval(Meta.parse(string(
                        "TestedRho.+=Theta",
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
                        "ViolatedRho.+=Theta",
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
                end
            end
            push!(TestedRhoBootstrap, TestedRho)
            push!(ViolatedRhoBootstrap, ViolatedRho)
        end
        push!(Tested, TestedRhoBootstrap)
        push!(Violated, ViolatedRhoBootstrap)
    end
end

TData = DataFrame(
    [rhovec mean([
        [
            sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
            for x in [1:1:length(rhovec);]
        ] for y in [1:1:m*runs;]
    ])],
)
rename!(TData, (:x1 => :rho))
rename!(TData, (:x2 => :meanest))
TDataAdd = DataFrame([
    [
        sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
        for x in [1:1:length(rhovec);]
    ] for y in [1:1:m*runs;]
])
TData = hcat(TData, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TData, (:x1 => :std))
TDataG2 = hcat(TData, TDataAdd)

Tested = []
Violated = []
for runno in [1:1:runs;]
    for bootstrapround in [1:1:m;]
        global TestedRhoBootstrap = []
        global ViolatedRhoBootstrap = []
        for rho in [1:1:length(rhovec);]
            global TestedRho = zeros(length(unique(bidderassignment)), 1)
            global ViolatedRho = zeros(length(unique(bidderassignment)), 1)
            for g in 3
                for auction in [1:1:length(group[g]);]
                    eval(Meta.parse(string(
                        "TestedRho.+=Theta",
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
                        "ViolatedRho.+=Theta",
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
                end
            end
            push!(TestedRhoBootstrap, TestedRho)
            push!(ViolatedRhoBootstrap, ViolatedRho)
        end
        push!(Tested, TestedRhoBootstrap)
        push!(Violated, ViolatedRhoBootstrap)
    end
end

TData = DataFrame(
    [rhovec mean([
        [
            sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
            for x in [1:1:length(rhovec);]
        ] for y in [1:1:m*runs;]
    ])],
)
rename!(TData, (:x1 => :rho))
rename!(TData, (:x2 => :meanest))
TDataAdd = DataFrame([
    [
        sum((Violated[y][x] ./ Tested[y][x]) .* [105, 15, 3]) / 123
        for x in [1:1:length(rhovec);]
    ] for y in [1:1:m*runs;]
])
TData = hcat(TData, [std(TDataAdd[x, :]) for x in [1:1:length(rhovec);]])
rename!(TData, (:x1 => :std))
TDataG3 = hcat(TData, TDataAdd)

default(lab = "")
ticks = pushfirst!(log.(TDataG1.rho[2:end]), -11)
ticklabels = string.(pushfirst!(log.(TDataG1.rho[2:end]), -Inf))
toplot = plot(
    pushfirst!(log.(TDataG1.rho[2:end]) .- 0.12, -11.12),
    TDataG1.meanest,
    ylim = (0.1, 0.8),
    yerror=1.96*TDataG1.std,
    seriestype = :scatter,
    #xticks = (ticks, ticklabels),
    label = "Fraction of BR Violations Auction Group 1",
    legend = :topleft,
)
plot!(
    pushfirst!(log.(TDataG2.rho[2:end]), -11),
    TDataG2.meanest,
    yerror=1.96*TDataG2.std,        
    ylim = (0.1, 0.8),
    seriestype = :scatter,
    #xticks = (ticks, ticklabels),
    label = "Fraction of BR Violations Auction Group 2",
    legend = :topleft,
)
plot!(
    pushfirst!(log.(TDataG3.rho[2:end]) .+ 0.12, -11 + 0.12),
    TDataG3.meanest,
    yerror=1.96*TDataG3.std,
    ylim = (0.1, 0.8),
    seriestype = :scatter,
    #xticks = (ticks, ticklabels),
    label = "Fraction of BR Violations Auction Group 3",
    legend = :topleft,
)

plot!(xlabel = L"\ln(\rho)")
plot(toplot)
savefig("ThetaEstAuctionGroups.pdf")

latexify(
    TDataG1[:, [:rho, :meanest, :std]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
latexify(
    TDataG2[:, [:rho, :meanest, :std]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
latexify(
    TDataG3[:, [:rho, :meanest, :std]],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)


## CE-Computation
CE(x, rho) = -1 / rho * log(0.5 * (1 + exp(-rho * x)))
