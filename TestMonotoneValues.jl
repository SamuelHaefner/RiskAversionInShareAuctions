include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

@load "WValues1.dat" W1


rhovec = [0,0,0]
#rhovec = [exp(-5),exp(-7),exp(-8)]

check = zeros(length(quotas))

for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    for auction in auctionset
        c = CheckSimpleBoundsDecreasing(auction, W1[g], bidderassignment, prices, rhovec)
        check[auction]=sum(c)/length(c)
    end
end

checklower = zeros(length(quotas))

for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    for auction in auctionset
        c = CheckSimpleBoundsDecreasingLower(auction, W1[g], bidderassignment, prices, rhovec)
        checklower[auction]=sum(c)/length(c)
    end
end

checkupper = zeros(length(quotas))

for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    for auction in auctionset
        c = CheckSimpleBoundsDecreasingUpper(auction, W1[g], bidderassignment, prices, rhovec)
        checkupper[auction]=sum(c)/length(c)
    end
end
