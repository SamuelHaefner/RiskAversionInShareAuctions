include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

n=AvgNoBidders(bidderassignment)
rhovec = [exp(x) for x in [-10:0.5:0;]]

# Estimation of WValuesTestMon.dat and BoundsTestMon[auction].dat only needs to be done once
m = 1
P = 500
W = []
for i in [1:1:length(group);]
    auctionset = group[i]
    prices = PriceBids(auctionset)
    @time push!(W,Wgamma(prices, auctionset, bidderassignment, n, m, P))
end
@save "WValuesTestMon.dat" W


for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    @time for auction in auctionset
        global BoundsAuction = []
        for i in [1:1:length(rhovec);]
            push!(BoundsAuction,EstimateSimpleBounds(auction, W[g], bidderassignment, prices, rhovec[i], m))
        end
        eval(Meta.parse(string("@save \"BoundsTestMon",auction,".dat\" BoundsAuction")))
    end
end

# x to be adapted, x in [1:1:length(rhovec);]!
# requires WValuesTestMon.dat and BoundsTestMon[auction].dat
x=1
rho=rhovec[x]
for R in [50,100,200]
    for auction in [1:1:39;]
        global T=TestIncreasingDiff(auction,activebidderindeces[auction],R,1,rho)
        eval(Meta.parse(string("@save \"TestMonA",auction,"R",R,"rho",x,".dat\" T")))
    end
end
