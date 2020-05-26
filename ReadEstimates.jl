## reads in the estimates of the different AvP-variables and produces the tables

include("Auxiliary.jl")
include("Grouping.jl")

using Latexify
using Plots


#########################################################################
## compute AvgP_\ell^{pre} and AvgP_u^{pre} for [auction] using the 
## [j]-th bootstrap round in [bounds]
#########################################################################

function AvgPpre(auction, bounds, j)
    AvgPpreUInd = []
    AvgPpreLInd = []
    for i in [1:1:length(bounds);]
        push!(
            AvgPpreUInd,
            IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][j][2],
            ),
        )
        push!(
            AvgPpreLInd,
            IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][j][1],
            ),
        )
    end
    return ([
        sum(AvgPpreLInd) / quotas[auction],
        sum(AvgPpreUInd) / quotas[auction],
    ])
end

#########################################################################
## compute AvgP_\ell^{post} and AvgP_u^{post} for [auction] using the 
## [j]-th bootstrap round in [bounds]
#########################################################################

function AvgPpost(auction, bounds, j)
    AvgPpostUInd = []
    AvgPpostLInd = []
    for i in [1:1:length(bounds);]
        push!(
            AvgPpostLInd,
            IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][j][1],
            ) - IntBid(
                0,
                qRec(activebidderindeces[auction][i], auction),
                qpBid(activebidderindeces[auction][i], auction),
            ),
        )
        push!(
            AvgPpostUInd,
            IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][j][2],
            ) - IntBid(
                0,
                qRec(activebidderindeces[auction][i], auction),
                qpBid(activebidderindeces[auction][i], auction),
            ),
        )
    end
    return ([
        sum(AvgPpostLInd) / quotas[auction],
        sum(AvgPpostUInd) / quotas[auction],
    ])
end

#########################################################################
## compute AvgP_\ell^{ratio} and AvgP_u^{ratio} for [auction] using the 
## [j]-th bootstrap round in [bounds]
#########################################################################

function AvgPrat(auction, bounds, j)
    n = length(findall(
        [
            qRec(activebidderindeces[auction][i], auction)
            for i in [1:1:length(bounds);]
        ] .> 0,
    ))
    AvgPratUInd = []
    AvgPratLInd = []
    for i in [1:1:length(bounds);]
        push!(
            AvgPratLInd,
            (
                IntV(
                    0,
                    qRec(activebidderindeces[auction][i], auction),
                    bounds[i][j][1],
                ) - IntBid(
                    0,
                    qRec(activebidderindeces[auction][i], auction),
                    qpBid(activebidderindeces[auction][i], auction),
                )
            ) / IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][j][1],
            ),
        )
        push!(
            AvgPratUInd,
            (
                IntV(
                    0,
                    qRec(activebidderindeces[auction][i], auction),
                    bounds[i][j][2],
                ) - IntBid(
                    0,
                    qRec(activebidderindeces[auction][i], auction),
                    qpBid(activebidderindeces[auction][i], auction),
                )
            ) / IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][j][2],
            ),
        )
    end
    return ([
        sum(AvgPratLInd[findall(convert(
            BitArray{1},
            -isnan.(AvgPratUInd) .+ 1,
        ))]) / n,
        sum(AvgPratUInd[findall(convert(
            BitArray{1},
            -isnan.(AvgPratLInd) .+ 1,
        ))]) / n,
    ])
end

#########################################################################
## returns the bounds for all bootstrap rounds and all auctions in 
## [auctiongroup], computed from [Bounds]
#########################################################################

function ComputeBounds(auctiongroup, Bounds)
    m=length(Bounds[1][1][1][1])
    Est=[]
    for j in [1:1:m;]
        AvgPpreL0 = []
        AvgPpreL1 = []
        AvgPpreL2 = []
        AvgPpreU0 = []
        AvgPpreU1 = []
        AvgPpreU2 = []

        AvgPpostL0 = []
        AvgPpostL1 = []
        AvgPpostL2 = []
        AvgPpostU0 = []
        AvgPpostU1 = []
        AvgPpostU2 = []

        AvgPratioL0 = []
        AvgPratioL1 = []
        AvgPratioL2 = []
        AvgPratioU0 = []
        AvgPratioU1 = []
        AvgPratioU2 = []

        for auction in auctiongroup
            # get group no. of auction
            g = findall([in(auction, group[x]) for x in [1:1:length(group);]])
            auctionindexingroup = findall(x -> x == auction, group[g[1]])
            t0 = Bounds[g[1]][auctionindexingroup[1]][1]
            t1 = Bounds[g[1]][auctionindexingroup[1]][2][1]
            t2 = Bounds[g[1]][auctionindexingroup[1]][2][2]

            Pre0 = AvgPpre(auction, t0, j)
            Pre1 = AvgPpre(auction, t1, j)
            Pre2 = AvgPpre(auction, t2, j)

            Post0 = AvgPpost(auction, t0, j)
            Post1 = AvgPpost(auction, t1, j)
            Post2 = AvgPpost(auction, t2, j)

            Ratio0 = AvgPrat(auction, t0, j)
            Ratio1 = AvgPrat(auction, t1, j)
            Ratio2 = AvgPrat(auction, t2, j)

            push!(AvgPpreL0, Pre0[1])
            push!(AvgPpreL1, Pre1[1])
            push!(AvgPpreL2, Pre2[1])

            push!(AvgPpreU0, Pre0[2])
            push!(AvgPpreU1, Pre1[2])
            push!(AvgPpreU2, Pre2[2])

            push!(AvgPpostL0, Post0[1])
            push!(AvgPpostL1, Post1[1])
            push!(AvgPpostL2, Post2[1])

            push!(AvgPpostU0, Post0[2])
            push!(AvgPpostU1, Post1[2])
            push!(AvgPpostU2, Post2[2])

            push!(AvgPratioL0, Ratio0[1])
            push!(AvgPratioL1, Ratio1[1])
            push!(AvgPratioL2, Ratio2[1])

            push!(AvgPratioU0, Ratio0[2])
            push!(AvgPratioU1, Ratio1[2])
            push!(AvgPratioU2, Ratio2[2])
        end
        push!(Est,DataFrame([
        AvgPpreL0,
        AvgPpreL1,
        AvgPpreL2,
        AvgPpreU0,
        AvgPpreU1,
        AvgPpreU2,
        AvgPpostL0,
        AvgPpostL1,
        AvgPpostL2,
        AvgPpostU0,
        AvgPpostU1,
        AvgPpostU2,
        AvgPratioL0,
        AvgPratioL1,
        AvgPratioL2,
        AvgPratioU0,
        AvgPratioU1,
        AvgPratioU2,
        ]))  
    end
    return(Est)
end

##########################################################################
## performs the actual computation of the bounds (returning means and 
## se) using the [Runs] runs of estimate filed Bounds[run].dat
###########################################################################

function ComputeBoundsMeanStd(Runs)
    Estimates = []
    for i in Runs
        eval(Meta.parse(string("@load \"Bounds", i, ".dat\" Bounds", i)))
        eval(Meta.parse(string("E = ComputeBounds([1:1:39;],Bounds", i, ")")))
        append!(Estimates, E)
    end

    Means = Matrix(undef, 39, 18)
    Std = Matrix(undef, 39, 18)

    for a in [1:1:39;]
        for e in [1:1:18;]
            Means[a, e] =
                mean([Estimates[x][a, e] for x in [1:1:length(Estimates);]])
            Std[a, e] = std([Estimates[x][a, e] for x in [1:1:length(Estimates);]])
        end
    end
    return ([DataFrame(Means), DataFrame(Std)])
end

######################################################################################
# read the actual estimates and construct the tables
######################################################################################

Bds = ComputeBoundsMeanStd(40)

# table with means
latexify(
    Bds[1],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )

# table with se
latexify(
    Bds[2],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )

# table with estimate overview
latexify(
    describe(Bds[1])[:,1:5],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )
