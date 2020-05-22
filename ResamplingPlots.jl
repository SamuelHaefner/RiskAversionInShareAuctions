include("Auxiliary.jl")
include("Grouping.jl")

using Plots
using LaTeXStrings
using Latexify

## construct bidderassignments plot
###########################################################################
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
ab=[AvgBid(x) for x in [1:1:123;]]
bidderplot = plot(
    [1:1:123;],
    ab ./ 1000,
    seriestype = :bar,
    linecolor = convert(Array{Int64,1}, bidderassignment),
    fillcolor = convert(Array{Int64,1}, bidderassignment),
    legend = false,
    xlabel = "Bidder",
    ylabel = "Average Bid [t]",
    size=(600,250),
    #title = "Bidder Groups",
)
savefig("ResamplingPlot2a.pdf")
# construct auctiongroups plot
############################################################################
auctiongroup = []
for i in [1:1:39;]
    if quotas[i] < 230000
        push!(auctiongroup, 3)
    elseif quotas[i] > 360000
        push!(auctiongroup, 1)
    else
        push!(auctiongroup, 2)
    end
end
auctionplot = plot(
    [1:1:39;],
    quotas ./ 1000,
    seriestype = :bar,
    fillcolor = convert(Array{Int64,1}, auctiongroup),
    linecolor = convert(Array{Int64,1}, auctiongroup),
    legend = false,
    xlabel = "Auction",
    ylabel = "Quota [t]",
    size=(600,250),
    #title = "Auction Groups",
)
savefig("ResamplingPlot2b.pdf")


# get summary statistics for the groups
###########################################################################
group = [1:1:length(unique(auctiongroup));]
group_me = [length(quotas[findall(auctiongroup.==x)]) for x in group]
group_min = [minimum(quotas[findall(auctiongroup.==x)]) for x in group]
group_mean = [mean(quotas[findall(auctiongroup.==x)]) for x in group]
group_max = [maximum(quotas[findall(auctiongroup.==x)]) for x in group]

groupdescr = DataFrame([group,group_me,group_min,group_mean,group_max])
rename!(groupdescr, Symbol.(["Group","Size","Min","Mean","Max"]))
latexify(groupdescr,env=:tabular,fmt=x -> round(x, sigdigits=5))

group = [1:1:length(unique(bidderassignment));]
group_me = [length(findall(bidderassignment.==x)) for x in group]
group_min = [minimum(ab[findall(bidderassignment.==x)]) for x in group]
group_mean = [mean(ab[findall(bidderassignment.==x)]) for x in group]
group_max = [maximum(ab[findall(bidderassignment.==x)]) for x in group]

groupdescr = DataFrame([group,group_me,group_min,group_mean,group_max])
rename!(groupdescr, Symbol.(["Group","Size","Min","Mean","Max"]))
latexify(groupdescr,env=:tabular,fmt=x -> round(x, sigdigits=5))


# plot resampling algorithm
############################################################################

#for every resampling round, construct residuals
auction=20
auctionset=findall(auctiongroup.==2)
P =100
prices=PriceBids(auctionset)
n=72

bidset = []
for j in auctionset
    for i in activebidderindeces[j]
        push!(bidset, qpbid(i, j))
    end
end

#construct P resampling indeces
resampling = []
for i in [1:1:P;]
    push!(resampling, rand([1:1:length(bidset);], n))
end

#for every resampling round, construct residual supply


#bootstrap round, for every price, retrive P opponent demand and
#fit gamma distribution to it
q = []
for r in [1:1:P;]
    qvals = []
    for p in prices
        x = 0
        for i in [1:1:n;]
            bid = bidset[resampling[r]][i]
            x += sum(bid[bid.pb.>=p, :qb])
        end
        append!(qvals,x)
    end
    push!(q,qvals)
end


default(lab="")
resamplingplot = plot(
    [(-q[x] .+ quotas[auction]) / 1000 for x in [1:1:P;]],
    prices,
    linetype = :steppre,
    color = :lightblue,
    title = "Resampled Residual Supply",
    xlabel = "Quantity [t]",
    ylabel = "Price [CHF]",
    size=(400,400)
)
xlims!(0, quotas[auction] / 1000)
ylims!(0,15)

plot!(
    (-q[1] .+ quotas[auction]) / 1000,
    prices,
    linetype = :steppre,
    color = :lightblue,
    label = "Residual Supply Functions",
    legend = :bottomright,
)

vline!([0], linecolor = :black)
plot(resamplingplot)
savefig("ResamplingPlot1.pdf")

## could add some bid functions
bid1 = qpbid(17, 33)
plot!(
    pushfirst!(copy(bid1.cumqb), 0) ./ 1000,
    pushfirst!(copy(bid1.pb), bid1.pb[1]),
    linetype = :steppre,
    color = :darkblue,
    label = "Bid Functions",
)

bid2 = qpbid(9, 28)
plot!(
    pushfirst!(copy(bid2.cumqb), 0) ./ 1000,
    pushfirst!(copy(bid2.pb), bid2.pb[1]),
    linetype = :steppre,
    color = :darkblue,
)

bid3 = qpbid(2, 34)
plot!(
    pushfirst!(copy(bid3.cumqb), 0) ./ 1000,
    pushfirst!(copy(bid3.pb), bid3.pb[1]),
    linetype = :steppre,
    color = :darkblue,
)


## plot a bunch of bid functions
auction = 36
bid = qpbid(activebidderindeces[auction][1], auction)
bidfctplot = plot(
    pushfirst!(copy(bid.cumqb), 0) ./ 1000,
    pushfirst!(copy(bid.pb), bid.pb[1]),
    linetype = :steppre,
    color = :darkblue,
    label = "Bid Functions",
    legend=:bottomright,title = "Bid Functions in an Auction",
    xlabel = "Quantity [t]",
    ylabel = "Price [CHF]",
    ylims = (0,15),
    size=(400,400),
)
for i in activebidderindeces[auction]
    bid = qpbid(i, auction)
    plot!(
        pushfirst!(copy(bid.cumqb), 0) ./ 1000,
        pushfirst!(copy(bid.pb), bid.pb[1]),
        linetype = :steppre,
        color = :darkblue,
    )
end
plot(bidfctplot)
savefig("ResamplingPlot1b.pdf")
