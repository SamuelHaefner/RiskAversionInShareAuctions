## Reads in the estimates of the different AvP-variables and produces the respective tables

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
    QReceived = []
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
        push!(QReceived,qRec(activebidderindeces[auction][i], auction))
    end

    AvgPpre = []
    push!(AvgPpre,[sum(AvgPpreLInd) / quotas[auction],sum(AvgPpreUInd) / quotas[auction]])
    for i in unique(bidderassignment)
        push!(AvgPpre,
            [sum(AvgPpreLInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))])/
            sum(QReceived[findall(in.(bidderassignment[activebidderindeces[auction]],i))]),
            sum(AvgPpreUInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))])/
            sum(QReceived[findall(in.(bidderassignment[activebidderindeces[auction]],i))])])
    end
    return(AvgPpre)
end

#########################################################################
## compute AvgP_\ell^{post} and AvgP_u^{post} for [auction] using the 
## [j]-th bootstrap round in [bounds]
#########################################################################

function AvgPpost(auction, bounds, j)
    AvgPpostUInd = []
    AvgPpostLInd = []
    QReceived = []
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
        push!(QReceived,qRec(activebidderindeces[auction][i], auction))
    end
    
    AvgPpost = []
    push!(AvgPpost,[sum(AvgPpostLInd) / quotas[auction],sum(AvgPpostUInd) / quotas[auction]])
    for i in unique(bidderassignment)
        push!(AvgPpost,
            [sum(AvgPpostLInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))])/
            sum(QReceived[findall(in.(bidderassignment[activebidderindeces[auction]],i))]),
            sum(AvgPpostUInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))])/
            sum(QReceived[findall(in.(bidderassignment[activebidderindeces[auction]],i))])])
    end
    return(AvgPpost)
end

#########################################################################
## compute AvgP_\ell^{ratio} and AvgP_u^{ratio} for [auction] using the 
## [j]-th bootstrap round in [bounds]
#########################################################################

function AvgPrat(auction, bounds, j)
    ### need to do this for every group!
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

    AvgPratList = []
    push!(AvgPratList,[
        sum(AvgPratLInd[findall(convert(
            BitArray{1},
            -isnan.(AvgPratLInd) .+ 1,
        ))]) / length(AvgPratLInd[findall(convert(
            BitArray{1},
            -isnan.(AvgPratLInd) .+ 1,
        ))]),
        sum(AvgPratUInd[findall(convert(
            BitArray{1},
            -isnan.(AvgPratUInd) .+ 1,
        ))]) / length(AvgPratUInd[findall(convert(
            BitArray{1},
            -isnan.(AvgPratUInd) .+ 1,
        ))]),
    ])

    for i in unique(bidderassignment)
        if length(findall(convert(BitArray{1}, 
            -isnan.(AvgPratLInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))]) .+ 1,
            ))) == 0
            push!(AvgPratList,[NaN,NaN])
        else
            push!(AvgPratList,[
            sum(AvgPratLInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))][findall(convert(
                BitArray{1},
                -isnan.(AvgPratLInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))]) .+ 1,
            ))]) / length(AvgPratLInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))][findall(convert(
                BitArray{1},
                -isnan.(AvgPratLInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))]) .+ 1,
            ))]),
            sum(AvgPratUInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))][findall(convert(
                BitArray{1},
                -isnan.(AvgPratUInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))]) .+ 1,
            ))]) / length(AvgPratUInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))][findall(convert(
                BitArray{1},
                -isnan.(AvgPratUInd[findall(in.(bidderassignment[activebidderindeces[auction]],i))]) .+ 1,
            ))])])
        end
        
    end
    return(AvgPratList)
end

#########################################################################
## returns the bounds for all bootstrap rounds and all auctions in 
## [auctiongroup], computed from [Bounds]
#########################################################################

function ComputeBounds(auctiongroup, Bounds,subg)
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

            push!(AvgPpreL0, Pre0[subg][1])
            push!(AvgPpreL1, Pre1[subg][1])
            push!(AvgPpreL2, Pre2[subg][1])

            push!(AvgPpreU0, Pre0[subg][2])
            push!(AvgPpreU1, Pre1[subg][2])
            push!(AvgPpreU2, Pre2[subg][2])

            push!(AvgPpostL0, Post0[subg][1])
            push!(AvgPpostL1, Post1[subg][1])
            push!(AvgPpostL2, Post2[subg][1])

            push!(AvgPpostU0, Post0[subg][2])
            push!(AvgPpostU1, Post1[subg][2])
            push!(AvgPpostU2, Post2[subg][2])

            push!(AvgPratioL0, Ratio0[subg][1])
            push!(AvgPratioL1, Ratio1[subg][1])
            push!(AvgPratioL2, Ratio2[subg][1])

            push!(AvgPratioU0, Ratio0[subg][2])
            push!(AvgPratioU1, Ratio1[subg][2])
            push!(AvgPratioU2, Ratio2[subg][2])
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
        ],:auto))  
    end
    return(Est)
end

##########################################################################
## performs the actual computation of the bounds (returning means and 
## se) using the [Runs] runs of estimate filed Bounds[run].dat
###########################################################################

function ComputeBoundsMeanStd(Runs,subg)
    Estimates = []
    for i in Runs
        eval(Meta.parse(string("@load \"Bounds", i, ".dat\" Bounds", i)))
        eval(Meta.parse(string("E = ComputeBounds([1:1:39;],Bounds", i, ",",subg,")")))
        eval(Meta.parse(string("Bounds", i, "= nothing")))
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
    return ([DataFrame(Means,:auto), DataFrame(Std,:auto)])
end


######################################################################################
# read the actual estimates and construct the tables
######################################################################################

# Compute bound estimates and std across all bidders (subg=1) 
Bds=0
Bds = ComputeBoundsMeanStd([1:1:40;],1)

# Construct data for the first table in Appendix D of the Supplementary Appendix (Estimates)
alldata_means=latexify(
    Bds[1],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )
write("AllEstMean.txt",alldata_means)

# Construct data for the second table in Appendix D of the Supplementary Appendix (std's)
alldata_se = latexify(
    Bds[2],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )
write("AllEstSE.txt",alldata_se)

# Construct data for Table 3 in the Main Text
alldata_summary=latexify(
    describe(Bds[1])[:,1:5],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )
write("AllEstSummary.txt",alldata_summary)

# Compute bound estimates and std across the different bidder groups
# (The different *.txt files below are used to construct Table 5 in the Main Text) 

# Compute bound estimates and std across all bidders of bidder group 1 (subg=2)
Bds = ComputeBoundsMeanStd([1:1:40;],2)

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
summary1 = latexify(
    describe(Bds[1][setdiff(1:end,31),:])[:,1:5],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )
write("Group1Summary.txt",summary1)

# Compute bound estimates and std across all bidders of bidder group 2 (subg=3)
Bds=0
Bds = ComputeBoundsMeanStd([1:1:40;],3)

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
summary2=latexify(
    describe(Bds[1][setdiff(1:end,31),:])[:,1:5],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
    )
write("Group2Summary.txt",summary2)

# Compute bound estimates and std across all bidders of bidder group 3 (subg=4)
Bds=0
Bds = ComputeBoundsMeanStd([1:1:40;],4)

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
summary3=latexify(
    describe(Bds[1][setdiff(1:end,31),:])[:,1:5],
    env = :tabular,
    fmt = x -> round(x, sigdigits = 4),
)
write("Group3Summary.txt",summary3)  