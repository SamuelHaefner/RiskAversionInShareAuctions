include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

n = AvgNoBidders(bidderassignment)

m = 5
P = 500

# load estimates of W obtained through EstimateWandTGeneric.jl
# indeces to be adapted!
@load "WValues1.dat" W1
W = W1

rhovec=[[0,0,0],[exp(-5),exp(-7),exp(-8)]]

Bounds = []
for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    BoundsGroup = []
    @time for auction in auctionset
        BoundsAuction = []
        for i in [1:1:length(rhovec);]
            simplebounds = EstimateSimpleBounds(
                auction,
                W[g],
                bidderassignment,
                prices,
                rhovec[i],
                m,
            )
            if i == 1
                push!(BoundsAuction, simplebounds)
            else
                push!(
                    BoundsAuction,
                    EstTighterBounds(
                        auction,
                        W[g],
                        bidderassignment,
                        prices,
                        simplebounds,
                        rhovec[i],
                        m,
                        10,
                        0.0001,
                    ),
                )
            end
        end
        push!(BoundsGroup, BoundsAuction)
    end
    push!(Bounds, BoundsGroup)

    # indeces to be adapted!
    Bounds1 = Bounds
    @save "Bounds1.dat" Bounds1
end
