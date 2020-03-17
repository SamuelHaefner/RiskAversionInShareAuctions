###########################################################################
## Plot bounds bds obtained from Algorithm (TighterBounds())
###########################################################################
function PlotTighterBounds(bds, initvl, initvu, bidder, auction, f)
    plot(
        pushfirst!(copy(initvl.qval), 0),
        [
            pushfirst!(copy(initvl.vval), initvl.vval[1]),
            pushfirst!(copy(initvu.vval), initvu.vval[1]),
        ],
        linetype = :steppre,
        color = :orange,
        label = ["Initial Conditions" ""],
    )
    plot!(
        pushfirst!(copy(bds[1].qval), 0),
        #pushfirst!(bds[1].vval, bds[1].vval[1]),
        [
            pushfirst!(copy(bds[1].vval), bds[1].vval[1]),
            pushfirst!(copy(bds[2].vval), bds[2].vval[1]),
        ],
        linetype = :steppre,
        color = :red,
        label = ["Estimated Bounds" ""],
        title = join(["Bidder ", bidder, " in Auction ", auction]),
        xlabel = "q",
        ylabel = "p",
    )
    bid = qpBid(bidder, auction)
    plot!(
        pushfirst!(copy(bid.cumqb), 0),
        pushfirst!(copy(bid.pb), bid.pb[1]),
        linetype = :steppre,
        color = :black,
        label = "Bid Schedule",
    )
    savefig(f)
end


#########################################################################
## compute AvgP_\ell^{pre} and AvgP_u^{pre}
#########################################################################

function AvgPpre(auction, bounds)
    AvgPpreUInd = []
    AvgPpreLInd = []
    for i in [1:1:length(bounds);]
        push!(
            AvgPpreUInd,
            IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][1][2],
            ),
        )
        push!(
            AvgPpreLInd,
            IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][1][1],
            ),
        )
    end
    return ([
        sum(AvgPpreLInd) / quotas[auction],
        sum(AvgPpreUInd) / quotas[auction],
    ])
end

function AvgPpost(auction, bounds)
    AvgPpostUInd = []
    AvgPpostLInd = []
    for i in [1:1:length(bounds);]
        push!(
            AvgPpostLInd,
            IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][1][1],
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
                bounds[i][1][2],
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

function AvgPrat(auction, bounds)
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
                    bounds[i][1][1],
                ) - IntBid(
                    0,
                    qRec(activebidderindeces[auction][i], auction),
                    qpBid(activebidderindeces[auction][i], auction),
                )
            ) / IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][1][1],
            ),
        )
        push!(
            AvgPratUInd,
            (
                IntV(
                    0,
                    qRec(activebidderindeces[auction][i], auction),
                    bounds[i][1][2],
                ) - IntBid(
                    0,
                    qRec(activebidderindeces[auction][i], auction),
                    qpBid(activebidderindeces[auction][i], auction),
                )
            ) / IntV(
                0,
                qRec(activebidderindeces[auction][i], auction),
                bounds[i][1][2],
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


function ComputeBounds(auctiongroup, Bounds)
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

        Pre0 = AvgPpre(auction, t0)
        Pre1 = AvgPpre(auction, t1)
        Pre2 = AvgPpre(auction, t2)

        Post0 = AvgPpost(auction, t0)
        Post1 = AvgPpost(auction, t1)
        Post2 = AvgPpost(auction, t2)

        Ratio0 = AvgPrat(auction, t0)
        Ratio1 = AvgPrat(auction, t1)
        Ratio2 = AvgPrat(auction, t2)

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
    return (DataFrame([
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


function ComputeBoundsMeanStd(Runs)
    Estimates = []
    for i in Runs
        eval(Meta.parse(string("@load \"Bounds", i, ".dat\" Bounds", i)))
        eval(Meta.parse(string("E = ComputeBounds([1:1:39;],Bounds", i, ")")))
        push!(Estimates, E)
    end

    Means = Matrix(undef, 39, 18)
    Std = Matrix(undef, 39, 18)

    for a in [1:1:39;]
        for e in [1:1:18;]
            Means[a, e] =
                mean([Estimates[x][a, e] for x in [1:1:length(Runs);]])
            Std[a, e] = std([Estimates[x][a, e] for x in [1:1:length(Runs);]])
        end
    end
    return ([DataFrame(Means), DataFrame(Std)])
end
