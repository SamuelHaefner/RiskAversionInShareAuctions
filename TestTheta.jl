include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

n = AvgNoBidders(bidderassignment)

m = 2
P = 500

@load "WValues1.dat" W1

W=W1

rhovec = [[exp(x),exp(x),exp(x)] for x in sort!(append!([-10:0.5:0;],-Inf))]

Theta = []
for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    ThetaGroup = []
    @time for auction in auctionset
        ThetaAuction = []
        for i in [1:1:length(rhovec);]
            bounds = EstimateSimpleBoundsNew(
                auction,
                W[g],
                bidderassignment,
                prices,
                rhovec[i],
                m,
            )
            push!(
                ThetaAuction,
                EstThetaNew(
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
@save "Theta1.dat" Theta1
