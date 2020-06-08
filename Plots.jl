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
ab = [AvgBid(x) for x in [1:1:123;]]
bidderplot = plot(
    [1:1:123;],
    ab ./ 1000,
    seriestype = :bar,
    linecolor = convert(Array{Int64,1}, bidderassignment),
    fillcolor = convert(Array{Int64,1}, bidderassignment),
    legend = false,
    xlabel = "Bidder",
    ylabel = "Average Bid [t]",
    size = (600, 250),
    # title = "Bidder Groups",
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
    size = (600, 250),
    # title = "Auction Groups",
)
savefig("ResamplingPlot2b.pdf")


# get summary statistics for the groups
###########################################################################
group = [1:1:length(unique(auctiongroup));]
group_me = [length(quotas[findall(auctiongroup .== x)]) for x in group]
group_min = [minimum(quotas[findall(auctiongroup .== x)]) for x in group]
group_mean = [mean(quotas[findall(auctiongroup .== x)]) for x in group]
group_max = [maximum(quotas[findall(auctiongroup .== x)]) for x in group]

groupdescr = DataFrame([group,group_me,group_min,group_mean,group_max])
rename!(groupdescr, Symbol.(["Group","Size","Min","Mean","Max"]))
latexify(groupdescr,env = :tabular,fmt = x->round(x, sigdigits = 5))

group = [1:1:length(unique(bidderassignment));]
group_me = [length(findall(bidderassignment .== x)) for x in group]
group_min = [minimum(ab[findall(bidderassignment .== x)]) for x in group]
group_mean = [mean(ab[findall(bidderassignment .== x)]) for x in group]
group_max = [maximum(ab[findall(bidderassignment .== x)]) for x in group]

groupdescr = DataFrame([group,group_me,group_min,group_mean,group_max])
rename!(groupdescr, Symbol.(["Group","Size","Min","Mean","Max"]))
latexify(groupdescr,env = :tabular,fmt = x->round(x, sigdigits = 5))


# plot resampling algorithm
############################################################################

# for every resampling round, construct residuals
auction = 20
auctionset = findall(auctiongroup .== 2)
P = 100
prices = PriceBids(auctionset)
n = 72

bidset = []
for j in auctionset
    for i in activebidderindeces[j]
        push!(bidset, qpbid(i, j))
    end
end

# construct P resampling indeces
resampling = []
for i in [1:1:P;]
    push!(resampling, rand([1:1:length(bidset);], n))
end

# for every resampling round, construct residual supply


# bootstrap round, for every price, retrive P opponent demand and
# fit gamma distribution to it
q = []
for r in [1:1:P;]
    qvals = []
    for p in prices
        x = 0
        for i in [1:1:n;]
            bid = bidset[resampling[r]][i]
            x += sum(bid[bid.pb .>= p, :qb])
        end
        append!(qvals, x)
    end
    push!(q, qvals)
end


default(lab = "")
resamplingplot = plot(
    [(-q[x] .+ quotas[auction]) / 1000 for x in [1:1:P;]],
    prices,
    linetype = :steppre,
    color = :lightblue,
    title = "Resampled Residual Supply",
    xlabel = "Quantity [t]",
    ylabel = "Price [CHF]",
    size = (400, 400)
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
    legend = :bottomright,title = "Bid Functions in an Auction",
    xlabel = "Quantity [t]",
    ylabel = "Price [CHF]",
    ylims = (0, 15),
    size = (400, 400),
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


###########################################################################
## Plot bounds bounds for [bidder] in auction [auction] from 
## [Bounds] saved in Bounds[run].dat; save plot in f
###########################################################################
function PlotTighterBounds(bidder, auction, f, m, runs)
    # get group no. of auction
    g = findall([in(auction, group[x]) for x in [1:1:length(group);]])
    auctionindexingroup = findall(x->x == auction, group[g[1]])

    # yvalues
    bdsu_y = []
    bdsl_y = []
    initvu_y = []
    initvl_y = []
    for runno in [1:1:runs;]
        eval(Meta.parse(string("@load \"Bounds", runno, ".dat\" Bounds", runno)))
        eval(Meta.parse(string("Bounds = Bounds", runno)))
        [push!(bdsu_y, copy(Bounds[g[1]][auctionindexingroup[1]][2][2][bidder][x][2].vval)) for x in [1:1:m;]]
        [push!(bdsl_y, copy(Bounds[g[1]][auctionindexingroup[1]][2][2][bidder][x][1].vval)) for x in [1:1:m;]]
        [push!(initvu_y, copy(Bounds[g[1]][auctionindexingroup[1]][2][1][bidder][x][2].vval)) for x in [1:1:m;]]
        [push!(initvl_y, copy(Bounds[g[1]][auctionindexingroup[1]][2][1][bidder][x][1].vval)) for x in [1:1:m;]]
    end

    [pushfirst!(bdsu_y[x], copy(bdsu_y[x][1])) for x in [1:1:m * runs;]]
    [pushfirst!(bdsl_y[x], copy(bdsl_y[x][1])) for x in [1:1:m * runs;]]
    [pushfirst!(initvu_y[x], copy(initvu_y[x][1])) for x in [1:1:m * runs;]]
    [pushfirst!(initvl_y[x], copy(initvl_y[x][1])) for x in [1:1:m * runs;]]

    # xvalues
    bds_x = Bounds[g[1]][auctionindexingroup[1]][2][2][bidder][1][2].qval
    pushfirst!(bds_x, 0)
    initv_x = Bounds[g[1]][auctionindexingroup[1]][2][1][bidder][1][2].qval
    pushfirst!(initv_x, 0)

    # bagged estimate 
    ebdsu_y = mean(bdsu_y)
    ebdsl_y = mean(bdsl_y)
    einitvu_y = mean(initvu_y)
    einitvl_y = mean(initvl_y)

    # errors bounds
    sebdsu_y = 1.95 * std(bdsu_y)
    sebdsl_y = 1.95 * std(bdsl_y)
    seinitvu_y = 1.95 * std(initvu_y)
    seinitvl_y = 1.95 * std(initvl_y)

    sebdsu_y_u = [quantile([bdsu_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(bdsu_y[1]);]] - ebdsu_y
    sebdsu_y_l = ebdsu_y - [quantile([bdsu_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(bdsu_y[1]);]]
    sebdsl_y_u = [quantile([bdsl_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(bdsl_y[1]);]] - ebdsl_y
    sebdsl_y_l = ebdsl_y - [quantile([bdsl_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(bdsl_y[1]);]]
    seinitvu_y_u = [quantile([initvu_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(initvu_y[1]);]] - einitvu_y
    seinitvu_y_l = einitvu_y - [quantile([initvu_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(initvu_y[1]);]]
    seinitvl_y_u = [quantile([initvl_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(initvl_y[1]);]] - einitvl_y
    seinitvl_y_l = einitvl_y - [quantile([initvl_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(initvl_y[1]);]]


    plot(
        initv_x,
        [einitvu_y,einitvl_y],
        ribbon = ([seinitvu_y_l,seinitvl_y_l], [seinitvu_y_u,seinitvl_y_u]),
        fillalpha = 0.1,
        # einitvu_y,
        # ribbon=(seinitvu_y_l,seinitvu_y_u),
        linetype = :steppre,
        color = :orange,
        label = ["Standard Bounds" ""],
    )
    plot!(
        bds_x,
        [ebdsu_y,ebdsl_y],
        ribbon = ([sebdsu_y_l,sebdsl_y_l], [sebdsu_y_u,sebdsl_y_u]),
        fillalpha = 0.1,
        # ribbon=[sebdsu_y,sebdsl_y],
        linetype = :steppre,
        color = :red,
        label = ["Tigher Bounds" ""],
        title = join(["Bidder ", bidder, " in Auction ", auction]),
        xlabel = "q",
        ylabel = "p",
    )
    bid = qpBid(activebidderindeces[auction][bidder], auction)
    plot!(
        pushfirst!(copy(bid.cumqb), 0),
        pushfirst!(copy(bid.pb), bid.pb[1]),
        linetype = :steppre,
        color = :black,
        label = "Bid Schedule",
    )
    savefig(f)
end
