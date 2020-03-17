# change working directory
#cd("Dropbox/Working Papers/Share Auctions/Risk Aversion/Replication/JuliaScripts")

# packages required
using CSV
using DataFrames
using Distributions
using BSON: @save, @load  ## @save and @load
using Roots


# read in data
bids = CSV.read(
    "setofbids.csv",
    delim = ',',
    missingstring = "NA",
    decimal = '.',
    types = [
        Int64,
        Int64,
        Int64,
        Int64,
        Float64,
        Float64,
        Float64,
        Float64,
        Float64,
    ],
)

# max number of price-quantity pairs allowed
K = 5

# unit of account: 100 corresponds to CHF
un = 100

# upper bound on type space (in CHF)
vupperbar = 23

# retreive the auction indeces of all auctions
# just use the frist 39 auctions (the last auction has only one participant)
auctionindeces = unique(bids[!, :auction])[1:39]

# retreive the bidder indeces of all registered bidders
bidderindeces = unique(bids[!, :bidder])

# retrieve, for every auction, the quota, the market clearing prices,
# the number of active bidders, and the indices of the active bidders
quotas = []
clearingprices = []
activebidders = []
activebidderindeces = []
for i in [1:1:length(auctionindeces);]
    push!(quotas, bids[bids.auction.==auctionindeces[i], :quotatot][1])
    push!(
        clearingprices,
        minimum(bids[(bids.auction.==auctionindeces[i]).&(bids.qr.>0), :pb]),
    )
    push!(
        activebidders,
        length(unique(bids[bids.auction.==auctionindeces[i], :bidder])),
    )
    push!(
        activebidderindeces,
        unique(bids[bids.auction.==auctionindeces[i], :bidder]),
    )
end


#############################################################################
# function that returns bid of [bidder] in [auction]
# in a dataframe with columns [qb,pb,cumqb]
#############################################################################
function qpbid(bidder, auction)
    # retreive price-quantity pairs
    bid = bids[
        (bids.bidder.==bidder).&(bids.auction.==auctionindeces[auction]),
        [:qb, :pb],
    ]

    # make sure it is sorted such that price points are decreasing
    sort!(bid, :pb, rev = true)

    # price points need to be divided by unit of account
    bid[!, 2] /= un

    # compute cumulated quantities
    cumqb = []
    for i in [1:1:size(bid, 1);]
        push!(cumqb, sum(bid[[1:1:i;], 1]))
    end

    # add to data frame and rename new column as cumqb
    cumqb = convert(Array{Float64,1}, cumqb)
    bid = hcat(bid, cumqb)
    rename!(bid, (:x1 => :cumqb))
end

#############################################################################
## returns average number of bidders in each cluster
#############################################################################
function AvgNoBidders(bidderassignment)
    avgno=[]
    for g in sort(unique(bidderassignment))
        n=[]
        for auction in [1:1:39;]
            push!(n,length(findall(in.(activebidderindeces[auction],(findall(in.(bidderassignment,g)),)))))
        end
        push!(avgno,ceil(mean(n)))
    end
    return convert(Array{Int64,1},avgno)
end

#############################################################################
# returns the bid-to-cover ration in an auction
#############################################################################
function BidToCover(auction)
    q=0
    for bidder in activebidderindeces[auction]
        q+=qpbid(bidder,auction).cumqb[end]
    end
    return q/quotas[auction]
end
# describe([BidToCover(x) for x in [1:1:39;]])

#############################################################################
# returns the revenue in an auction
#############################################################################
function Revenue(auction)
    revenue = sum(bids[(bids.auction.==auctionindeces[auction]),:qr,])
end
# describe([Revenue(x) for x in [1:1:39;]])

##############################################################################
# returns the received quantity of a bidder in a given auction
##############################################################################
function qrec(bidder,auction)
    return sum(bids[
        (bids.bidder.==bidder).&(bids.auction.==auctionindeces[auction]),
        :qr,
    ])
end

##############################################################################
# returns the received quantity share of a bidder in a given auction
##############################################################################
function qsharerec(bidder,auction)
    return sum(bids[
        (bids.bidder.==bidder).&(bids.auction.==auctionindeces[auction]),
        :qperc,
    ])
end


##############################################################################
# returns the indeces in which a bidder is active
###############################################################################
function ActiveAuctions(bidder)
    activeauctions=[]
    for i in [1:1:39;]
        if in(bidder,activebidderindeces[i]) push!(activeauctions,i) end
    end
    return activeauctions
end
# describe([length(ActiveAuctions(x)) for x in [1:1:123;]])
# findall([length(ActiveAuctions(x)) for x in [1:1:123;]].>=35)


##############################################################################
# return the average bid of a bidder
##############################################################################
function AvgBid(bidder)
    qbid = []
    for i in ActiveAuctions(bidder)
        push!(qbid,qpbid(bidder,i).cumqb[end])
    end
    return mean(qbid)
end
# describe([AvgBid(x) for x in [1:1:123;]])

########################################################################
# return share of succesfull bidders in an auction
#########################################################################
function ShareSuccBidders(auction)
    succ=0
    for i in activebidderindeces[auction]
        if qrec(i,auction)>0 succ+=1 end
    end
    return succ/activebidders[auction]
end
# describe([ShareSuccBidders(x) for x in [1:1:123;]])

#########################################################################
# return the success rate of a bidder
#########################################################################
function SuccessRate(bidder)
    succ=0
    for i in ActiveAuctions(bidder)
        if qrec(bidder,i)>0 succ+=1 end
    end
    return succ/length(ActiveAuctions(bidder))
end
# describe([SuccessRate(x) for x in [1:1:123;]])

# group auctions
group=[]
push!(group,findall(x->x>360000, quotas))
push!(group,findall(x->(x<=360000) && (x>=230000), quotas))
push!(group,findall(x->x<230000, quotas))

# assign bidders to three different groups based on average quantity bid
bidderassignment=[]
for i in [1:1:123;]
    if AvgBid(i) < 15000 push!(bidderassignment,1)
    elseif AvgBid(i) > 50000 push!(bidderassignment,3)
    else push!(bidderassignment,2) end
end


#########################################################################
# function returning [m] bootstrap estimates of the parameters [.,.]
# of the gamma distribution of opponent demand, Q-S(p),
# when [n] opponent bidders are in the auction,
# for all price points in [prices].
#
# the resampling algorithm redraws from the set of bids in [auctionset]
# [P] times.
#########################################################################
function Wgamma(prices, auctionset, bidderassignment, n, m, P)

    # construct array of relevant bids, one array for every bidder group
    bidset = []
    for g in sort(unique(bidderassignment))
        bidsetgroup=[]
        for j in auctionset
            for i in activebidderindeces[j][findall(in.(activebidderindeces[j],(findall(in.(bidderassignment,g)),)))]
                push!(bidsetgroup, qpbid(i, j))
            end
        end
        push!(bidset,bidsetgroup)
    end

    #construct m bootstrap indeces
    bootstrapped = []
    for i in [1:1:m;]
        push!(bootstrapped, [rand([1:1:length(bidset[x]);], length(bidset[x])) for x in sort(unique(bidderassignment))])
    end

    #construct P resampling indeces, drawing [n(1),...,n(g)] bids, to be used on every
    #bootstrap sample
    resampling = []
    for g in sort(unique(bidderassignment))
        resamplinggroup = []
        nn=[]
        for h in sort(unique(bidderassignment))
            if h==g push!(nn,n[h]-1) else push!(nn,n[h]) end
        end
        for i in [1:1:P;]
            push!(resamplinggroup, [rand([1:1:length(bidset[x]);], nn[x]) for x in sort(unique(bidderassignment))])
        end
        push!(resampling,resamplinggroup)
    end


    #for every group g, every bootstrap round, for every price,  retrive P opponent demand and
    #fit gamma distribution to it
    W=[]
    for g in sort(unique(bidderassignment))
        dist = []
        for b in [1:1:m;]
            #retrieve bids
            relbids=[]
            for r in [1:1:P;]
                samplebid=[]
                for h in [1:1:length(unique(bidderassignment));]
                    append!(samplebid,bidset[h][bootstrapped[b][h]][resampling[g][r][h]])
                end
                push!(relbids,samplebid)
            end
            #for every price, return estimate of demand
            z = []
            for p in prices
                y=[]
                for r in [1:1:P;]
                    push!(y,sum([sum(relbids[r][x][relbids[r][x].pb.>=p, :qb]) for x in [1:1:length(relbids[r]);]]))
                end
                push!(z, try
                    fit(Gamma, convert(Array{Float64,1}, y))
                catch
                    missing
                end)
            end
            push!(dist,z)
        end
        push!(W,dist)
    end
    return(W)
end

#########################################################################
# function returning [m] bootstrap estimates of the parameters [.,.]
# of the log normal distribution of opponent demand, Q-S(p),
# when [n] opponent bidders are in the auction, f
# or all price points in [prices].
#
# the resampling algorithm redraws from the set of bids in [auctionset]
# [P] times.
#########################################################################
function Wlnorm(prices, auctionset, bidderassignment, n, m, P)
    # construct array of relevant bids, one array for every bidder group
    bidset = []
    for g in sort(unique(bidderassignment))
        bidsetgroup=[]
        for j in auctionset
            for i in activebidderindeces[j][findall(in.(activebidderindeces[j],(findall(in.(bidderassignment,g)),)))]
                push!(bidsetgroup, qpbid(i, j))
            end
        end
        push!(bidset,bidsetgroup)
    end

    #construct m bootstrap indeces
    bootstrapped = []
    for i in [1:1:m;]
        push!(bootstrapped, [rand([1:1:length(bidset[x]);], length(bidset[x])) for x in sort(unique(bidderassignment))])
    end

    #construct P resampling indeces, drawing [n(1),...,n(g)] bids, to be used on every
    #bootstrap sample
    resampling = []
    for g in sort(unique(bidderassignment))
        resamplinggroup = []
        nn=[]
        for h in sort(unique(bidderassignment))
            if h==g push!(nn,n[h]-1) else push!(nn,n[h]) end
        end
        for i in [1:1:P;]
            push!(resamplinggroup, [rand([1:1:length(bidset[x]);], nn[x]) for x in sort(unique(bidderassignment))])
        end
        push!(resampling,resamplinggroup)
    end


    #for every group g, every bootstrap round, for every price,  retrive P opponent demand and
    #fit gamma distribution to it
    W=[]
    for g in sort(unique(bidderassignment))
        dist = []
        for b in [1:1:m;]
            #retrieve bids
            relbids=[]
            for r in [1:1:P;]
                samplebid=[]
                for h in [1:1:length(unique(bidderassignment));]
                    append!(samplebid,bidset[h][bootstrapped[b][h]][resampling[g][r][h]])
                end
                push!(relbids,samplebid)
            end
            #for every price, return estimate of demand
            z = []
            for p in prices
                y=[]
                for r in [1:1:P;]
                    push!(y,sum([sum(relbids[r][x][relbids[r][x].pb.>=p, :qb]) for x in [1:1:length(relbids[r]);]]))
                end
                push!(z, try
                    fit(LogNormal, convert(Array{Float64,1}, y))
                catch
                    missing
                end)
            end
            push!(dist,z)
        end
        push!(W,dist)
    end
    return(W)
end

#########################################################################
## return value of \beta_bid(q)
#########################################################################
function StepBid(q, bid)
    if q <= minimum(bid[!, :cumqb])
        return bid[1, :pb]
    end
    if q > maximum(bid[!, :cumqb])
        return 0
    end
    return bid[length(bid[bid.cumqb.<q, :cumqb])+1, :pb]
end


#########################################################################
## return number of step in \beta_bid at q
#########################################################################
function NoStepBid(q, bid)
    if q <= minimum(bid[!, :cumqb])
        return 1
    end
    if q > maximum(bid[!, :cumqb])
        return length(bid.cumqb)
    end
    return length(bid[bid.cumqb.<q, :cumqb]) + 1
end

#########################################################################
## return value of v(q)
#########################################################################
function StepV(q, v)
    if q <= minimum(v[!, :qval])
        return v[1, :vval]
    end
    if q > maximum(v[!, :qval])
        return 0
    end
    return v[length(v[v.qval.<q, :qval])+1, :vval]
end

#########################################################################
## return value of \int_a^b\beta_bid(q)dq
#########################################################################
function IntBid(a, b, bid)
    x = bid[(bid.cumqb.>a).&(bid.cumqb.<b), :cumqb]
    push!(x, b)
    pushfirst!(x, a)
    y = []
    y = [StepBid(z, bid) for z in x[2:end]]
    x = x[2:end] - x[1:end-1]
    return round(x' * y; digits = 4)
end

#########################################################################
## return value of \int_a^b v(q)dq
#########################################################################
function IntV(a, b, v)
    x = v[(v.qval.>a).&(v.qval.<b), :qval]
    push!(x, b)
    pushfirst!(x, a)
    y = []
    y = [StepV(z, v) for z in x[2:end]]
    x = x[2:end] - x[1:end-1]
    return round(x' * y; digits = 4)
end

#########################################################################
## return value of \overline{\Pi}^j_i(bid,v) when the quota is Q
## WPar is an array, containing the estimated parameters for steps j=1,...,k
#########################################################################
function PiOverline(bidstep, bid, v, WPar, rho, Q, n)
    if rho == 0
        return 0
    end
    if bidstep >= length(bid.pb)
        return 0
    end

    # terms under integrand
    f(x) = exp(-rho * (IntV(bid.cumqb[bidstep], x, v) -
         IntBid(bid.cumqb[bidstep], x, bid))) *
    (StepV(x, v) - StepBid(x, bid)) * cdf(WPar[NoStepBid(x, bid)], Q - x)

    ### integrate f from bid.cumqb[j] to bid.cumqb[end]
    ### use intervals given with intpts
    intpts = v[v.qval.>bid.cumqb[bidstep], :qval]
    pushfirst!(intpts, bid.cumqb[bidstep])
    val = []
    for i in [1:1:length(intpts)-1;]
        d = (intpts[i+1] - intpts[i]) / n
        push!(val, d * sum([f(x) for x in [intpts[i]+d:d:intpts[i+1];]]))
    end

    return try
        sum(val)
    catch
        0
    end
end

##########################################################################
## returns dataframe with simple bounds, v=[vub,vlb,q]
## WPar is an array, containing the estimated parameters for steps j=1,...,k
##########################################################################
function SimpleBound(bid, WPar, rho, Q)
    ## start with last step
    vlb = DataFrame(vval = bid.pb[end], qval = bid.cumqb[end])

    ## no information in case of just one step
    if length(bid.pb) == 1
        return [vlb,DataFrame(vval = vupperbar, qval = bid.cumqb[end])]
    end

    ## using vlb on last step gives \PiOverline_i^j = 0
    if (ismissing(WPar[length(bid.pb)-1]) || ismissing(WPar[length(bid.pb)]))
        return [DataFrame(vval = bid.pb, qval = bid.cumqb),DataFrame(vval = fill(vupperbar,length(bid.cumqb)), qval = bid.cumqb)]
    else
        vub = DataFrame(
            vval = max(
                bid.pb[length(bid.pb)-1],
                min(
                    vupperbar,
                    bid.pb[length(bid.pb)-1] +
                    (bid.pb[length(bid.pb)-1] - bid.pb[length(bid.pb)]) *
                    cdf(WPar[length(bid.pb)], Q - bid.cumqb[length(bid.pb)-1]) /
                    (cdf(
                        WPar[length(bid.pb)-1],
                        Q - bid.cumqb[length(bid.pb)-1],
                    ) - cdf(
                        WPar[length(bid.pb)],
                        Q - bid.cumqb[length(bid.pb)-1],
                    )),
                ),
            ),
            qval = bid.cumqb[end],
        )
    end

    ## move backwards through steps, if more than two steps
    if length(bid.pb) >= 3
        for i in [length(bid.pb)-1:-1:2;]
            # check whether distribution estimate produced missing,
            # if so return early
            if (ismissing(WPar[i]) || ismissing(WPar[i+1]))
                vlbnew = bid.pb[1:i]
                qvalnew = bid.cumqb[1:i]
                # DataFrame(vval=fill(vupperbar,2), qval=[100, 200])
                append!(vlb, DataFrame(vval=vlbnew,qval=qvalnew))
                sort!(vlb, :qval)
                vubnew = fill(vupperbar,i)
                append!(vub, DataFrame(vval=vubnew,qval=qvalnew))
                sort!(vub, :qval)
                # construct decreasing bounds \overline{v}_i and \underline{v}_i;
                # take most conservative bound from estimated vub and vul
                if length(vlb.qval) >= 2
                    for i in [1:1:length(vlb.qval);]
                        vub[i, :vval] = maximum(vub[i:end, :vval])
                        vlb[i, :vval] = max(bid.pb[i],minimum(vlb[1:i, :vval]))
                    end
                end
                return [vlb, vub]
            else
                vlbnew = min(
                    vupperbar,
                    max(
                        bid.pb[i],
                        bid.pb[i] +
                        (bid.pb[i] - bid.pb[i+1]) *
                        (cdf(WPar[i+1], Q - bid.cumqb[i]) -
                         rho * PiOverline(i, bid, vub, WPar, rho, Q, 100)[1]) /
                        (cdf(WPar[i], Q - bid.cumqb[i]) -
                         cdf(WPar[i+1], Q - bid.cumqb[i])),
                    ),
                )
            end
            push!(vlb, [vlbnew bid.cumqb[i]])
            sort!(vlb, :qval)
            # check whether distribution estimate produced missing,
            # if so, return early
            if (ismissing(WPar[i-1]) || ismissing(WPar[i]))
                qvalnew = bid.cumqb[1:i]
                vubnew = fill(vupperbar,i)
                append!(vub, DataFrame(vval=vubnew,qval=qvalnew))
                sort!(vub, :qval)
                # construct decreasing bounds \overline{v}_i and \underline{v}_i;
                # take most conservative bound from estimated vub and vul
                if length(vlb.qval) >= 2
                    for i in [1:1:length(vlb.qval);]
                        vub[i, :vval] = maximum(vub[i:end, :vval])
                        vlb[i, :vval] = max(bid.pb[i],minimum(vlb[1:i, :vval]))
                    end
                end
                return [vlb, vub]
            else
                vubnew = max(
                    bid.pb[i-1],
                    min(
                        vupperbar,
                        bid.pb[i-1] +
                        (bid.pb[i-1] - bid.pb[i]) *
                        (cdf(WPar[i], Q - bid.cumqb[i-1]) -
                         rho * PiOverline(i - 1, bid, vlb, WPar, rho, Q,100)[1]) /
                        (cdf(WPar[i-1], Q - bid.cumqb[i-1]) -
                         cdf(WPar[i], Q - bid.cumqb[i-1])),
                    ),
                )
            end
            push!(vub, [vubnew bid.cumqb[i]])
            sort!(vub, :qval)
        end
    end

    ## first step
    if (ismissing(WPar[1]) || ismissing(WPar[2]))
        vlbnew = bid.pb[1]
    else
        vlbnew = min(
            vupperbar,
            max(
                bid.pb[1],
                bid.pb[1] +
                (bid.pb[1] - bid.pb[2]) * (cdf(WPar[2], Q - bid.cumqb[1]) -
                 rho * PiOverline(1, bid, vub, WPar, rho, Q,100)[1]) /
                (cdf(WPar[1], Q - bid.cumqb[1]) -
                 cdf(WPar[2], Q - bid.cumqb[1])),
            ),
        )
    end
    push!(vlb, [vlbnew bid.cumqb[1]])
    sort!(vlb, :qval)
    push!(vub, [vupperbar bid.cumqb[1]])
    sort!(vub, :qval)

    # construct decreasing bounds \overline{v}_i and \underline{v}_i;
    # take most conservative bound from estimated vub and vul
    if length(vlb.qval) >= 2
        for i in [1:1:length(vlb.qval);]
            vub[i, :vval] = maximum(vub[i:end, :vval])
            vlb[i, :vval] = max(bid.pb[i],minimum(vlb[1:i, :vval]))
        end
    end

    return [vlb, vub]
end

############################################################################
## returns w(p,q) where WPar contains the parameters of W(p+dp,q) and W(p,q)
############################################################################
function w(q, WPar, dp)
    return (cdf(WPar[1], q) - cdf(WPar[2], q)) / dp
end

############################################################################
## returns all submitted prices in [auctionset] in ascending order
############################################################################
function PriceBids(auctionset)
    pricebidsarr = [bids[(bids.auction.==auctionindeces[x]), :pb] for x in auctionset]
    pricebids = []
    for i in [1:1:length(pricebidsarr);]
        append!(pricebids, pricebidsarr[i])
    end
    return sort!(unique(pricebids / un))
end


###########################################################################
## Estimates simple bounds for all bidders in [auction] using bids from
## [auctionset] assuming a risk preference rho. The resampling algorithm
## uses [P] redraws on [m] bootstrap estimates. The function returns the
## price bids, the corresponding W, as well as the bounds for every bidder
## every bootstrap run.
###########################################################################
function EstimateSimpleBounds(auction, W, bidderassignments,  prices, rho, m)
    bounds = []
    for bidder in activebidderindeces[auction]
        boundsbidder = []
        for bootstraprun in [1:1:m;]
            bid = qpbid(bidder, auction)
            g = bidderassignment[bidder]
            WPar = W[g][bootstraprun][sort(
                findall(in.(prices, (bid.pb,))),
                rev = true,
            )]
            boundsbidderbootstraprun = SimpleBound(
                bid,
                WPar,
                rho,
                quotas[auction],
            )
            push!(boundsbidder, boundsbidderbootstraprun)
        end
        push!(bounds, boundsbidder)
    end
    return bounds
end

##########################################################################
## estimates share of violations
##########################################################################

function EstTheta(auction, W, bidderassignment, prices, bounds, rho, m)
    tested = zeros(length(unique(bidderassignment)),m)  #define matrix, one colum for every bootrap round, one row for every bidder group
    violated = zeros(length(unique(bidderassignment)),m) #define matrix, one column for every bootrap round, one row for every bidder group
    for bidder in [1:1:activebidders[auction];]
        g = bidderassignment[bidder]
        bid = qpbid(activebidderindeces[auction][bidder], auction)
        for bootstraprun in [1:1:m;]
            for bidstep in [length(bid.pb):-1:1;]
                uc = FOC(
                    bidstep,
                    bid,
                    bounds[bidder][bootstraprun][1],
                    W,
                    g,
                    prices,
                    rho,
                    quotas[auction],
                    bootstraprun,
                    100
                )
                lc = FOC(
                    bidstep,
                    bid,
                    bounds[bidder][bootstraprun][2],
                    W,
                    g,
                    prices,
                    rho,
                    quotas[auction],
                    bootstraprun,
                    100
                )
                if (uc > 0 || lc < 0)
                    violated[g,bootstraprun] += 1
                end
                tested[g,bootstraprun] += 1
            end
        end
    end
    return [tested, violated]
end


##########################################################################
## returns the value of the FOC
##########################################################################
function FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
    WPar = W[group][bootstraprun][sort(findall(in.(prices, (bid.pb,))), rev = true)]
    if maximum(findall(in.(prices, (bid.pb,)))) < length(prices)
        WParPlus = W[group][bootstraprun][sort(
            findall(in.(prices, (bid.pb,))) .+ 1,
            rev = true,
        )]
    else
        WParPlus = W[group][bootstraprun][sort(
            findall(in.(prices, (bid.pb,))) .+ 1,
            rev = true,
        )[2:end]]
        pushfirst!(WParPlus, missing)
    end

    dp = try
        prices[findall(in.(prices, (bid.pb,))).+1][bidstep] -
        prices[findall(in.(prices, (bid.pb,)))][bidstep]
    catch
        return 0
    end

    if ismissing(WParPlus[bidstep]) || ismissing(WPar[bidstep])
        return 0
    end
    if bidstep >= 2
        f(x) =
            exp(rho * (IntV(x, bid.cumqb[bidstep], v) -
                 IntBid(x, bid.cumqb[bidstep], bid))) *
            ((StepV(x, v) - StepBid(x, bid)) *
             (w(Q - x, [WParPlus[bidstep], WPar[bidstep]], dp) +
              rho * (x - bid.cumqb[bidstep-1]) *
              cdf(WPar[NoStepBid(x, bid)], Q - x)) -
             cdf(WPar[NoStepBid(x, bid)], Q - x))

        ### integrate f from bid.cumqb[j] to bid.cumqb[end]
        ### use intervals given with intpts
        intpts = v[
            (v.qval.>bid.cumqb[bidstep-1]).&(v.qval.<bid.cumqb[bidstep]),
            :qval,
        ]
        pushfirst!(intpts, bid.cumqb[bidstep-1])
        push!(intpts, bid.cumqb[bidstep])
        val = []
        for i in [1:1:length(intpts)-1;]
            d = (intpts[i+1] - intpts[i]) / n
            push!(val, d * sum([f(x) for x in [intpts[i]+d:d:intpts[i+1];]]))
        end
        return sum(val) + rho * (bid.cumqb[bidstep] - bid.cumqb[bidstep-1]) *
               PiOverline(bidstep, bid, v, WPar, rho, Q, n)
    else
        g(x) =
            exp(rho * (IntV(x, bid.cumqb[bidstep], v) -
                 IntBid(x, bid.cumqb[bidstep], bid))) *
            ((StepV(x, v) - StepBid(x, bid)) *
             (w(Q - x, [WParPlus[bidstep], WPar[bidstep]], dp) +
              rho * x * cdf(WPar[NoStepBid(x, bid)], Q - x)) -
             cdf(WPar[NoStepBid(x, bid)], Q - x))

        ### integrate f from bid.cumqb[j] to bid.cumqb[end]
        ### use intervals given with intpts
        intpts = v[(v.qval.>0).&(v.qval.<bid.cumqb[bidstep]), :qval]
        pushfirst!(intpts, 0)
        push!(intpts, bid.cumqb[bidstep])
        val = []
        for i in [1:1:length(intpts)-1;]
            d = (intpts[i+1] - intpts[i]) / n
            push!(val, d * sum([g(x) for x in [intpts[i]+d:d:intpts[i+1];]]))
        end
        return sum(val) + rho * bid.cumqb[bidstep] *
               PiOverline(bidstep, bid, v, WPar, rho, Q, n)
    end
end



###########################################################################
## functions to estimate tighter bounds
###########################################################################

function VarPhiU(q, v, vl)
    vnew = DataFrame(
        vval = max.(
            fill(convert(Float64, v), length(vl[vl.qval.<=q, :qval])),
            vl[vl.qval.<=q, :vval],
        ),
        qval = vl[vl.qval.<=q, :qval],
    )
    if !in(q, vl[:, :qval])
        push!(vnew, [max.(v,vl[vl.qval.>q, :vval][1]), q])
    end
    append!(vnew, vl[vl.qval.>q, :])
end

function VarPhiL(q, v, vu)
    vnew = DataFrame(
        vval = min.(
            fill(convert(Float64, v), length(vu[vu.qval.>=q, :qval])),
            vu[vu.qval.>=q, :vval],
        ),
        qval = vu[vu.qval.>=q, :qval],
    )
    if !in(q, vu[:, :qval])
        push!(vnew, [vu[vu.qval.<q, :vval][end], q])
    end
    append!(vnew, vu[vu.qval.<q, :])
    sort!(vnew,:qval)
end

#############################################################################
## Runs Algorithm 2, finding tighter bounds
#############################################################################

function TighterBounds(
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
    vl = copy(initvl)
    vu = copy(initvu)
    vl[!,:vval] = convert(Array{Float64,1}, vl[!,:vval])
    vl[!,:qval] = convert(Array{Float64,1}, vl[!,:qval])
    vu[!,:vval] = convert(Array{Float64,1}, vu[!,:vval])
    vu[!,:qval] = convert(Array{Float64,1}, vu[!,:qval])


    # partition step into 10 subintervals
    for bidstep in [length(bid.pb):-1:1;]
        print(string("b",bidstep))
        vuint = []
        vlint = []
        d = bid.qb[bidstep] / 10
        # obtain q values
        if bidstep == 1
            qvals = [d:d:bid.cumqb[1];]
            push!(vuint, IntV(0, bid.cumqb[1], vu))
            push!(vlint, IntV(0, bid.cumqb[1], vl))
        else
            qvals = [bid.cumqb[bidstep-1]+d:d:bid.cumqb[bidstep];]
            push!(vuint, IntV(bid.cumqb[bidstep-1], bid.cumqb[bidstep], vu))
            push!(vlint, IntV(bid.cumqb[bidstep-1], bid.cumqb[bidstep], vl))
        end

        for repetition in [2:1:maxiter;]
            print(string("r",repetition))
            # obtain new values for vu
            vvals = Array{Float64}(undef,0)
            if repetition == 2
                n = 100
            else
                n = 10
            end
            for q in qvals
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
                                    xatol=0.0001
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
                                    xatol=0.0001
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
                DataFrame(vval = pushfirst!(vvals,StepV(qvals[1],vu))[1:end-1], qval = qvals),
            )
            sort!(vu, :qval)

            # obtain new values for vl
            vvals = Array{Float64}(undef,0)
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
                    if fl(StepV(q, initvu)) < 0 || fl(StepV(q, initvl)) > 0
                        push!(vvals, StepV(q, initvl))
                    else
                        push!(
                            vvals,
                            try
                                fzero(
                                    fl,
                                    (StepV(q, initvl), StepV(q, initvu)),
                                    xatol=0.0001
                                )
                            catch
                                StepV(q, initvl)
                            end,
                        )
                    end
                else
                    ## check values at upper and lower end of interval
                    if fl(StepV(q, initvu)) < 0 ||
                       fl(maximum([
                        StepV(bid.cumqb[bidstep], vl),
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
                                    xatol=0.0001
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
    return [vl, vu]
end

function EstTighterBounds(auction, W, bidderassignment, prices, bounds, rho, m, maxiter,tolerance)
    tighterbounds = []
    for bidder in [1:1:activebidders[auction];]
    #for bidder in [1:1:2;]
        print(bidder)
        boundsbidder = []
        #only compute tighter bounds if allocated with positive amount
        if qrec(activebidderindeces[auction][bidder],auction) > 0
            bid = qpbid(activebidderindeces[auction][bidder], auction)
            g = bidderassignment[bidder]
            for bootstraprun in [1:1:m;]
                initvl=bounds[bidder][bootstraprun][1]
                initvu=bounds[bidder][bootstraprun][2]
                boundsbidderbootstraprun = TighterBounds(bid,initvl,initvu,W,g,prices,rho,quotas[auction],bootstraprun,maxiter,tolerance)
                push!(boundsbidder, boundsbidderbootstraprun)
            end
        else
            for bootstraprun in [1:1:m;]
                initvl=bounds[bidder][bootstraprun][1]
                initvu=bounds[bidder][bootstraprun][2]
                push!(boundsbidder, [initvl,initvu])
            end
        end
        push!(tighterbounds, boundsbidder)
    end
    return [bounds,tighterbounds]
end


###########################################################################
## Plot bounds bds obtained from Algorithm (TighterBounds())
###########################################################################
function PlotTighterBounds(bds,initvl,initvu,bidder,auction,f)
    plot(
    pushfirst!(copy(initvl.qval),0),
    [
     pushfirst!(copy(initvl.vval),initvl.vval[1]),
     pushfirst!(copy(initvu.vval),initvu.vval[1])
    ],
    linetype = :steppre,
    color = :orange,
    label = ["Initial Conditions" ""],
    )
    plot!(
    pushfirst!(copy(bds[1].qval), 0),
    #pushfirst!(bds[1].vval, bds[1].vval[1]),
    [
     pushfirst!(copy(bds[1].vval), bds[1].vval[1]),
     pushfirst!(copy(bds[2].vval), bds[2].vval[1])
    ],
    linetype = :steppre,
    color = :red,
    label = ["Estimated Bounds" ""],
    title = join(["Bidder ", bidder, " in Auction ", auction]),
    xlabel = "q",
    ylabel = "p"
    )
    bid=qpbid(bidder,auction)
    plot!(
    pushfirst!(copy(bid.cumqb),0),
    pushfirst!(copy(bid.pb),bid.pb[1]),
    linetype = :steppre,
    color = :black,
    label = "Bid Schedule",
    )
    savefig(f)
end


#########################################################################
## compute AvgP_\ell^{pre} and AvgP_u^{pre}
#########################################################################

function AvgPpre(auction,bounds)
    AvgPpreUInd = []
    AvgPpreLInd = []
    for i in [1:1:length(bounds);]
        push!(AvgPpreUInd,IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][2]))
        push!(AvgPpreLInd,IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][1]))
    end
    return([sum(AvgPpreLInd)/quotas[auction],sum(AvgPpreUInd)/quotas[auction]])
end

function AvgPpost(auction,bounds)
    AvgPpostUInd = []
    AvgPpostLInd = []
    for i in [1:1:length(bounds);]
        push!(AvgPpostLInd,IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][1])-IntBid(0,qrec(activebidderindeces[auction][i],auction),qpbid(activebidderindeces[auction][i],auction)))
        push!(AvgPpostUInd,IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][2])-IntBid(0,qrec(activebidderindeces[auction][i],auction),qpbid(activebidderindeces[auction][i],auction)))
    end
    return([sum(AvgPpostLInd)/quotas[auction],sum(AvgPpostUInd)/quotas[auction]])
end

function AvgPrat(auction,bounds)
    n=length(findall([qrec(activebidderindeces[auction][i],auction) for i in [1:1:length(bounds);]].>0))
    AvgPratUInd = []
    AvgPratLInd = []
    for i in [1:1:length(bounds);]
        push!(AvgPratLInd,(IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][1])-IntBid(0,qrec(activebidderindeces[auction][i],auction),qpbid(activebidderindeces[auction][i],auction)))/IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][1]))
        push!(AvgPratUInd,(IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][2])-IntBid(0,qrec(activebidderindeces[auction][i],auction),qpbid(activebidderindeces[auction][i],auction)))/IntV(0,qrec(activebidderindeces[auction][i],auction),bounds[i][1][2]))
    end
    return([sum(AvgPratLInd[findall(convert(BitArray{1},-isnan.(AvgPratUInd).+1))])/n,sum(AvgPratUInd[findall(convert(BitArray{1},-isnan.(AvgPratLInd).+1))])/n])
end


function ComputeBounds(auctiongroup,Bounds)
    AvgPpreL0=[]
    AvgPpreL1=[]
    AvgPpreL2=[]
    AvgPpreU0=[]
    AvgPpreU1=[]
    AvgPpreU2=[]

    AvgPpostL0=[]
    AvgPpostL1=[]
    AvgPpostL2=[]
    AvgPpostU0=[]
    AvgPpostU1=[]
    AvgPpostU2=[]

    AvgPratioL0=[]
    AvgPratioL1=[]
    AvgPratioL2=[]
    AvgPratioU0=[]
    AvgPratioU1=[]
    AvgPratioU2=[]

    for auction in auctiongroup
        # get group no. of auction
        g = findall([in(auction,group[x]) for x in [1:1:length(group);]])
        auctionindexingroup = findall(x->x==auction,group[g[1]])
        t0 = Bounds[g[1]][auctionindexingroup[1]][1]
        t1 = Bounds[g[1]][auctionindexingroup[1]][2][1]
        t2 = Bounds[g[1]][auctionindexingroup[1]][2][2]

        Pre0=AvgPpre(auction,t0)
        Pre1=AvgPpre(auction,t1)
        Pre2=AvgPpre(auction,t2)

        Post0=AvgPpost(auction,t0)
        Post1=AvgPpost(auction,t1)
        Post2=AvgPpost(auction,t2)

        Ratio0=AvgPrat(auction,t0)
        Ratio1=AvgPrat(auction,t1)
        Ratio2=AvgPrat(auction,t2)

        push!(AvgPpreL0,Pre0[1])
        push!(AvgPpreL1,Pre1[1])
        push!(AvgPpreL2,Pre2[1])

        push!(AvgPpreU0,Pre0[2])
        push!(AvgPpreU1,Pre1[2])
        push!(AvgPpreU2,Pre2[2])

        push!(AvgPpostL0,Post0[1])
        push!(AvgPpostL1,Post1[1])
        push!(AvgPpostL2,Post2[1])

        push!(AvgPpostU0,Post0[2])
        push!(AvgPpostU1,Post1[2])
        push!(AvgPpostU2,Post2[2])

        push!(AvgPratioL0,Ratio0[1])
        push!(AvgPratioL1,Ratio1[1])
        push!(AvgPratioL2,Ratio2[1])

        push!(AvgPratioU0,Ratio0[2])
        push!(AvgPratioU1,Ratio1[2])
        push!(AvgPratioU2,Ratio2[2])
    end
    return(DataFrame([AvgPpreL0,AvgPpreL1,AvgPpreL2,AvgPpreU0,AvgPpreU1,AvgPpreU2,AvgPpostL0,AvgPpostL1,AvgPpostL2,AvgPpostU0,AvgPpostU1,AvgPpostU2,AvgPratioL0,AvgPratioL1,AvgPratioL2,AvgPratioU0,AvgPratioU1,AvgPratioU2]))
end


function ComputeBoundsMeanStd(Runs)
    Estimates = []
    for i in Runs
        eval(Meta.parse(string("@load \"Bounds",i,".dat\" Bounds",i)))
        eval(Meta.parse(string("E = ComputeBounds([1:1:39;],Bounds",i,")")))
        push!(Estimates,E)
    end

    Means = Matrix(undef,39,18)
    Std = Matrix(undef,39,18)

    for a in [1:1:39;]
        for e in [1:1:18;]
            Means[a,e] = mean([Estimates[x][a,e] for x in [1:1:length(Runs);]])
            Std[a,e] =  std([Estimates[x][a,e] for x in [1:1:length(Runs);]])
        end
    end
    return([DataFrame(Means),DataFrame(Std)])
end
