include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

using Plots
using LaTeXStrings
using Latexify

## Construct bidderassignments plot (left panel in Figure 3 in the Main Text)
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

# Construct auctiongroups plot (left panel in Figure 6 in the Main Text)
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


# Get summary statistics for the groups (used for right panels in Figure 6 and 3, respectively)
# Note: this does not automatically produce an *.txt output; rather
# the latexify()-output below needs to be copied by hand.
###########################################################################
group_a = [1:1:length(unique(auctiongroup));]
group_me = [length(quotas[findall(auctiongroup .== x)]) for x in group_a]
group_min = [minimum(quotas[findall(auctiongroup .== x)]) for x in group_a]
group_mean = [mean(quotas[findall(auctiongroup .== x)]) for x in group_a]
group_max = [maximum(quotas[findall(auctiongroup .== x)]) for x in group_a]

groupdescr = DataFrame([group_a,group_me,group_min,group_mean,group_max],:auto)
rename!(groupdescr, Symbol.(["Group","Size","Min","Mean","Max"]))
latexify(groupdescr,env = :tabular,fmt = x->round(x, sigdigits = 5))

group_b = [1:1:length(unique(bidderassignment));]
group_me = [length(findall(bidderassignment .== x)) for x in group_b]
group_min = [minimum(ab[findall(bidderassignment .== x)]) for x in group_b]
group_mean = [mean(ab[findall(bidderassignment .== x)]) for x in group_b]
group_max = [maximum(ab[findall(bidderassignment .== x)]) for x in group_b]

groupdescr = DataFrame([group_b,group_me,group_min,group_mean,group_max],:auto)
rename!(groupdescr, Symbol.(["Group","Size","Min","Mean","Max"]))
latexify(groupdescr,env = :tabular,fmt = x->round(x, sigdigits = 5))


# Plot resampling algorithm (right panel in Figure 7 in the Main Text)
############################################################################

auction = 20
auctionset = findall(auctiongroup .== 2)
P = 100
prices = PriceBids(auctionset)
n = 72

bidset = []
for j in auctionset
    for i in activebidderindeces[j]
        push!(bidset, qpBid(i, j))
    end
end

# construct P resampling indeces
resampling = []
for i in [1:1:P;]
    push!(resampling, rand([1:1:length(bidset);], n))
end

# for every bootstrap round and for every price, retrive P opponent demand 
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


## Plot a bunch of bid functions (left panel in Figure 7 in the Main Text)
#####################################################################################
auction = 36
bid = qpBid(activebidderindeces[auction][1], auction)
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
    local bid = qpBid(i, auction)
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
## [m=5] in data, corresponds to number of bootstrap rounds per run
## [runs=40] corresponds to the number of rounds
##
## !! the function will return an err if number of points in intial bounds 
## disagree between bootstrap rounds; this happens when W or w is missing 
## for a point; which in turn happens when the price bid was highest among 
## all submitted prices.
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
    bds_x = copy(Bounds[g[1]][auctionindexingroup[1]][2][2][bidder][1][2].qval)
    pushfirst!(bds_x, 0)
    initv_x = copy(Bounds[g[1]][auctionindexingroup[1]][2][1][bidder][1][2].qval)
    pushfirst!(initv_x, 0)

    # bagged estimate 
    ebdsu_y = mean(bdsu_y)
    ebdsl_y = mean(bdsl_y)
    einitvu_y = mean(initvu_y)
    einitvl_y = mean(initvl_y)

    # confidence bands
    sebdsu_y_u = [quantile([bdsu_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(bdsu_y[1]);]]
    sebdsu_y_l = [quantile([bdsu_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(bdsu_y[1]);]]
    sebdsl_y_u = [quantile([bdsl_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(bdsl_y[1]);]]
    sebdsl_y_l = [quantile([bdsl_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(bdsl_y[1]);]]
    seinitvu_y_u = [quantile([initvu_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(initvu_y[1]);]]
    seinitvu_y_l = [quantile([initvu_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(initvu_y[1]);]]
    seinitvl_y_u = [quantile([initvl_y[x][y] for x in [1:1:m * runs;]], 0.95) for y in [1:1:length(initvl_y[1]);]]
    seinitvl_y_l = [quantile([initvl_y[x][y] for x in [1:1:m * runs;]], 0.05) for y in [1:1:length(initvl_y[1]);]]

    plot(
        initv_x,
        #[[einitvu_y einitvu_y], [einitvl_y einitvl_y]],
        [einitvu_y einitvu_y einitvl_y einitvl_y],
        #fillrange = [[einitvu_y.+seinitvu_y_l,einitvu_y.-seinitvl_y_l], [einitvl_y.+seinitvu_y_u,einitvl_y.-seinitvl_y_u]],
        fillrange = [seinitvu_y_l seinitvl_y_l seinitvu_y_u seinitvl_y_u],
        fillalpha = 0.1,
        # einitvu_y,
        # ribbon=(seinitvu_y_l,seinitvu_y_u),
        linetype = :steppre,
        color = :orange,
        label = ["Standard Bounds" "" "" ""],
    )
    plot!(
        bds_x,
        [ebdsu_y ebdsu_y ebdsl_y ebdsl_y],
        fillrange = [sebdsu_y_l sebdsu_y_u sebdsl_y_l sebdsl_y_u],
        fillalpha = 0.1,
        # ribbon=[sebdsu_y,sebdsl_y],
        linetype = :steppre,
        color = :red,
        label = ["Tighter Bounds" "" "" ""], 
        title = join(["Bidder ", bidder, " in Auction ", auction, " from Group No. ", bidderassignment[activebidderindeces[auction][bidder]]]),
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

## Produce figure for Figure 5 in the Main Text
########################################################################################
PlotTighterBounds(1,20,"BoundsBidder1Auction20.pdf",5,40)
PlotTighterBounds(8,25,"BoundsBidder8Auction25.pdf",5,40)

###############################################################
# TighterBoundsIterations(bid, initvl, initvu, W, g, prices, rho, Q, bootstraprun, maxiter, tolerance)
###############################################################
#### Description
# Computes tighter upper and lower bounds from initial conditions ```initvl``` 
# and ```initivu``` by runing Algorithm 2. To be used below.
#### Arguments
# ```bid``` -- bid function  
# ```initvl``` -- marginal profit function, the initial condition for the lower bound  
# ```initvu``` -- marginal profit function, the initial condition for the upper bound  
# ```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
# ```g``` -- natural number, indicating the group number of the auction  
# ```prices``` -- vector of prices used for the estimation of W  
# ```rho``` -- positive real number, the risk preference   
# ```bootstraprun``` -- index of boostratp run we look at  
# ```maxiter``` -- maximum number of fixed point iterations   
# ```tolerance``` -- tolerance level (used in iteration)  
#### Return value
# A list of lists of dataframes, [[vlb,vlb,...,vlb],[vub,...,vub]], where each instance of vub is 
# the data frame containing the upper bound and vlb is the data frame containing lower bound
# in one iteration of the algorithm.
################################################################

function TighterBoundsIterations(
    bid,
    initvl,
    initvu,
    W,
    g,
    prices,
    rho,
    Q,
    bootstraprun,
    maxiter,
    tolerance,
    )

    # initial conditions, convert initvl and initvu
    vl = copy(initvl)
    vu = copy(initvu)
    vl[!, :vval] = convert(Array{Float64,1}, vl[!, :vval])
    vl[!, :qval] = convert(Array{Float64,1}, vl[!, :qval])
    vu[!, :vval] = convert(Array{Float64,1}, vu[!, :vval])
    vu[!, :qval] = convert(Array{Float64,1}, vu[!, :qval])

    vuiterations=[]
    vliterations=[]

    push!(vuiterations,vu)
    push!(vliterations,vl)

    ## go through steps in bid function in reversing order
    for bidstep in [length(bid.pb):-1:1;]
        # initialize bounds to be iterated 
        vuint = []
        vlint = []
        # partition step into 5 subintervals
        d = bid.qb[bidstep] / 5
        if bidstep == 1
            qvals = [d:d:bid.cumqb[1];]
            push!(vuint, IntV(0, bid.cumqb[1], vu))
            push!(vlint, IntV(0, bid.cumqb[1], vl))
        else
            qvals = [bid.cumqb[bidstep-1]+d:d:bid.cumqb[bidstep];]
            push!(vuint, IntV(bid.cumqb[bidstep-1], bid.cumqb[bidstep], vu))
            push!(vlint, IntV(bid.cumqb[bidstep-1], bid.cumqb[bidstep], vl))
        end

        # iteration of bounds, at most [maxiter] iterations
        for repetition in [2:1:maxiter;]
            # obtain new values for vu
            vvals = Array{Float64}(undef, 0)
            if repetition == 2
                n = 100
            else
                n = 10
            end
            for q in qvals
                # function for which root is to be determined
                fu(x) = FOC(
                    bidstep,
                    bid,
                    VarPhiU(q, x, vl),
                    W,
                    g,
                    prices,
                    rho,
                    Q,
                    bootstraprun,
                    n,
                )
                # check if we are looking at the last step 
                # (range when calling fzero differs)
                if bidstep == length(bid.pb)
                    ## check values at upper and lower end of interval
                    if fu(StepV(q, initvu)) < 0
                        push!(vvals, StepV(q, initvu))
                    elseif fu(StepV(q, initvl)) > 0
                        push!(vvals, StepV(q, initvl))
                    else
                        push!(
                            vvals,
                            try
                                fzero(
                                    fu,
                                    (StepV(q, initvl), StepV(q, initvu)),
                                    xatol = 0.0001,
                                )
                            catch
                                StepV(q, initvu)
                            end,
                        )
                    end
                else
                    ## check values at upper and lower end of interval
                    if fu(StepV(q, initvu)) < 0
                        push!(vvals, StepV(q, initvu))
                    elseif fu(StepV(bid.cumqb[bidstep] + 1, vu)) > 0
                        push!(vvals, StepV(bid.cumqb[bidstep] + 1, vu))
                    else
                        push!(
                            vvals,
                            try
                                fzero(
                                    fu,
                                    (
                                        StepV(bid.cumqb[bidstep] + 1, vu),
                                        StepV(q, initvu),
                                    ),
                                    xatol = 0.0001,
                                )
                            catch
                                StepV(q, initvu)
                            end,
                        )
                    end
                end
            end
            #replace values in vu with those in [vvals,qvals]
            vu = append!(
                vu[findall(.!in.(vu.qval, (qvals,))), :],
                DataFrame(
                    vval = pushfirst!(vvals, StepV(qvals[1], vu))[1:end-1],
                    qval = qvals,
                ),
            )
            sort!(vu, :qval)

            push!(vuiterations,vu)

            # obtain new values for vl
            vvals = Array{Float64}(undef, 0)
            if repetition == 2
                n = 100
            else
                n = 10
            end
            for q in qvals
                fl(x) = FOC(
                    bidstep,
                    bid,
                    VarPhiL(q, x, vu),
                    W,
                    g,
                    prices,
                    rho,
                    Q,
                    bootstraprun,
                    n,
                )
                if bidstep == length(bid.pb)
                    ## check values at upper and lower end of interval
                    if fl(StepV(q, initvu)) < 0 
                        push!(vvals,StepV(q, initvu))
                    elseif fl(maximum([StepV(bid.cumqb[bidstep], vl), 
                            StepV(bid.cumqb[bidstep] + 1, vl),
                            ])) > 0
                        push!(
                            vvals,
                            maximum([
                                StepV(bid.cumqb[bidstep], vl),
                                StepV(bid.cumqb[bidstep] + 1, vl),
                            ]),
                        )
                    else
                        push!(
                            vvals,
                            try
                                fzero(
                                    fl,
                                    (StepV(q, initvl), StepV(q, initvu)),
                                    xatol = 0.0001,
                                )
                            catch
                                StepV(q, initvl)
                            end,
                        )
                    end
                else
                    ## check values at upper and lower end of interval
                    if fl(StepV(q, initvu)) < 0 
                        push!(vvals,StepV(q, initvu))
                    elseif fl(maximum([StepV(bid.cumqb[bidstep], vl), 
                            StepV(bid.cumqb[bidstep] + 1, vl),
                            ])) > 0
                        push!(
                            vvals,
                            maximum([
                                StepV(bid.cumqb[bidstep], vl),
                                StepV(bid.cumqb[bidstep] + 1, vl),
                            ]),
                        )
                    else
                        push!(
                            vvals,
                            try
                                fzero(
                                    fl,
                                    (
                                        maximum([
                                            StepV(bid.cumqb[bidstep], vl),
                                            StepV(bid.cumqb[bidstep] + 1, vl),
                                        ]),
                                        StepV(q, initvu),
                                    ),
                                    xatol = 0.0001,
                                )
                            catch
                                maximum([
                                    StepV(bid.cumqb[bidstep], vl),
                                    StepV(bid.cumqb[bidstep] + 1, vl),
                                ])
                            end,
                        )
                    end
                end
            end
            #replace values in vl with those in [vvals,qvals]
            vl = append!(
                vl[findall(.!in.(vl.qval, (qvals,))), :],
                DataFrame(vval = vvals, qval = qvals),
            )
            sort!(vl, :qval)

            push!(vliterations,vl)
            
            if bidstep == 1
                push!(vuint, IntV(0, bid.cumqb[1], vu))
                push!(vlint, IntV(0, bid.cumqb[1], vl))
            else
                push!(vuint, IntV(bid.cumqb[bidstep-1], bid.cumqb[bidstep], vu))
                push!(vlint, IntV(bid.cumqb[bidstep-1], bid.cumqb[bidstep], vl))
            end

            if abs(vuint[repetition] - vuint[repetition-1]) +
               abs(vlint[repetition] - vlint[repetition-1]) < tolerance
                break
            end
        end
    end
    return [vliterations, vuiterations]
end

## Construct plot showing the iterations of the fixed point algorithm
## (Figure 8 in the Main Text)
##########################################################################

auction=25
auctiongroup=2
rho=[exp(-5),exp(-7),exp(-8)]
m=1
prices = PriceBids(group[auctiongroup])

@load "WValues1.dat" W1
W=W1

simplebounds = EstimateSimpleBoundsRobust(
                auction,
                W[auctiongroup],
                bidderassignment,
                prices,
                rho,
                m,
)

bidder=8

g = bidderassignment[findall(in.(bidderindeces,activebidderindeces[auction][bidder]))][1]
bid = qpBid(activebidderindeces[auction][bidder], auction)

Q=quotas[auction]
maxiter=10
tolerance=0.001

initvl = simplebounds[bidder][1][1]
initvu = simplebounds[bidder][1][2]

bootstraprun=1

test = TighterBoundsIterations(
    bid,
    initvl,
    initvu,
    W[auctiongroup],
    g,
    prices,
    rho[g],
    Q,
    bootstraprun,
    maxiter,
    tolerance,
)

xpoints=[]
ypointsu=[]
ypointsl=[]
for i in [1:1:length(test[1]);]
    push!(xpoints,pushfirst!(copy(test[1][i].qval),0))
    push!(ypointsl,pushfirst!(copy(test[1][i].vval),test[1][i].vval[1]))
    push!(ypointsu,pushfirst!(copy(test[2][i].vval),test[2][i].vval[1]))
end

p=plot(
    xpoints[end],
    [ypointsl[end],ypointsu[end]],
    linetype = :steppre,
    color=RGBA(0.0, 0.0, 0, 1), legend=false,
    title = join(["Bidder ", bidder, " in Auction ", auction]),
    xlabel = "q",
    ylabel = "p",
)
d=0.9/(length(test[1])-1)
for i in [1:1:length(test[1]);]
    eval(Meta.parse(string("p=plot!(xpoints[",i,"],[ypointsl[",i,"],ypointsu[",i,"]],linetype = :steppre,color=RGBA(",1-i*d,",",1-i*d,",",1-i*d,", 1))")))
end
plot(p)
savefig("BoundsAlgorithmConverge.pdf")