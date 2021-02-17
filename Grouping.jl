# This file determines the bidder groups
# and the auction groups used for the estimation

# Cf. Readme.md for more information.

# group auctions
group = []
push!(group, findall(x -> x > 360000, quotas))
push!(group, findall(x -> (x <= 360000) && (x >= 230000), quotas))
push!(group, findall(x -> x < 230000, quotas))

# assign bidders to three different groups based on average quantity bid
bidderassignment = []
for i in [1:1:123;]
    if AvgBid(i) < 15000
        push!(bidderassignment, 1)
    elseif AvgBid(i) > 50000
        push!(bidderassignment, 3)
    else
        push!(bidderassignment, 2)
    end
end
