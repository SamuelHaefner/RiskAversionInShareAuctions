# Estimation of WValuesTestMon.dat and BoundsTestMon[auction].dat 
# required for TestMonGenereic.jl

include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")
include("TestMon.jl")

n = AvgNoBidders(bidderassignment)
rhovec = [[exp(x),exp(x),exp(x)] for x in sort!(append!([-10:0.5:0;],-Inf))]

m = 1
P = 500
W = []
for i in [1:1:length(group);]
    auctionset = group[i]
    prices = PriceBids(auctionset)
    @time push!(W, Wgamma(prices, auctionset, bidderassignment, n, m, P))
end
@save "WValuesTestMon.dat" W


for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    @time for auction in auctionset
        global BoundsAuction = []
        for i in [1:1:length(rhovec);]
            push!(BoundsAuction, EstimateSimpleBoundsRobust(auction, W[g], bidderassignment, prices, rhovec[i], m))
        end
        eval(Meta.parse(string("@save \"BoundsTestMon", auction, ".dat\" BoundsAuction")))
    end
end