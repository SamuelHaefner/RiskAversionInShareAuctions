include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

n = AvgNoBidders(bidderassignment)

m = 5
P = 500
W = []
for i in [1:1:length(group);]
    auctionset = group[i]
    prices = PriceBids(auctionset)
    @time push!(W, Wgamma(prices, auctionset, bidderassignment, n, m, P))
end

# indeces to be adapted!
W1 = W
@save "WValues1.dat" W1

rhovec = sort!(append!([exp(x) for x in [-10:0.5:0;]], 0))

Theta = []
for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    ThetaGroup = []
    @time for auction in auctionset
        ThetaAuction = []
        for i in [1:1:length(rhovec);]
            bounds = EstimateSimpleBounds(
                auction,
                W[g],
                bidderassignment,
                prices,
                rhovec[i],
                m,
            )
            push!(
                ThetaAuction,
                EstTheta(
                    auction,
                    W[g],
                    bidderassignment,
                    prices,
                    bounds,
                    rhovec[i],
                    m,
                ),
            )
        end
        push!(ThetaGroup, ThetaAuction)
    end
    push!(Theta, ThetaGroup)
end

# indeces to be adapted!
Theta1 = Theta
@save "Theta1.dat" Theta1