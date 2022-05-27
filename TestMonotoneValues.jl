include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

@load "WValues1.dat" W1


#rhovec = [0,0,0]
rhovec = [exp(-5),exp(-7),exp(-8)]

auctionset = group[2]
auction = auctionset[15]
prices = PriceBids(auctionset)

c = CheckSimpleBoundsDecreasing(auction, W1[2], bidderassignment, prices, rhovec)
sum(c)/length(c)
