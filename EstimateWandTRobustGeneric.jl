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
    @time push!(W, Wlnorm(prices, auctionset, bidderassignment, n, m, P))
end

# indeces to be adapted!
W1 = W
@save "WValuesRobust1.dat" W1

rhovec = [[exp(x),exp(x),exp(x)] for x in sort!(append!([-10:0.5:0;],-Inf))]

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
                    rhovec[i][1],
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
@save "ThetaRobust1.dat" Theta1
