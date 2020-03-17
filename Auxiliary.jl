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
function qpBid(bidder, auction)
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
    avgno = []
    for g in sort(unique(bidderassignment))
        n = []
        for auction in [1:1:39;]
            push!(
                n,
                length(findall(in.(
                    activebidderindeces[auction],
                    (findall(in.(bidderassignment, g)),),
                ))),
            )
        end
        push!(avgno, ceil(mean(n)))
    end
    return convert(Array{Int64,1}, avgno)
end

#############################################################################
# returns the bid-to-cover ration in an auction
#############################################################################
function BidToCover(auction)
    q = 0
    for bidder in activebidderindeces[auction]
        q += qpBid(bidder, auction).cumqb[end]
    end
    return q / quotas[auction]
end
# describe([BidToCover(x) for x in [1:1:39;]])

#############################################################################
# returns the revenue in an auction
#############################################################################
function Revenue(auction)
    revenue = sum(bids[(bids.auction.==auctionindeces[auction]), :qr])
end
# describe([Revenue(x) for x in [1:1:39;]])

##############################################################################
# returns the received quantity of a bidder in a given auction
##############################################################################
function qRec(bidder, auction)
    return sum(bids[
        (bids.bidder.==bidder).&(bids.auction.==auctionindeces[auction]),
        :qr,
    ])
end

##############################################################################
# returns the received quantity share of a bidder in a given auction
##############################################################################
function qShareRec(bidder, auction)
    return sum(bids[
        (bids.bidder.==bidder).&(bids.auction.==auctionindeces[auction]),
        :qperc,
    ])
end


##############################################################################
# returns the indeces in which a bidder is active
###############################################################################
function ActiveAuctions(bidder)
    activeauctions = []
    for i in [1:1:39;]
        if in(bidder, activebidderindeces[i])
            push!(activeauctions, i)
        end
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
        push!(qbid, qpBid(bidder, i).cumqb[end])
    end
    return mean(qbid)
end
# describe([AvgBid(x) for x in [1:1:123;]])

########################################################################
# return share of succesfull bidders in an auction
#########################################################################
function ShareSuccBidders(auction)
    succ = 0
    for i in activebidderindeces[auction]
        if qRec(i, auction) > 0
            succ += 1
        end
    end
    return succ / activebidders[auction]
end
# describe([ShareSuccBidders(x) for x in [1:1:123;]])

#########################################################################
# return the success rate of a bidder
#########################################################################
function SuccessRate(bidder)
    succ = 0
    for i in ActiveAuctions(bidder)
        if qRec(bidder, i) > 0
            succ += 1
        end
    end
    return succ / length(ActiveAuctions(bidder))
end
# describe([SuccessRate(x) for x in [1:1:123;]])

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
    f(x) =
        exp(
            -rho * (
                IntV(bid.cumqb[bidstep], x, v) -
                IntBid(bid.cumqb[bidstep], x, bid)
            ),
        ) *
        (StepV(x, v) - StepBid(x, bid)) *
        cdf(WPar[NoStepBid(x, bid)], Q - x)

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
    pricebidsarr =
        [bids[(bids.auction.==auctionindeces[x]), :pb] for x in auctionset]
    pricebids = []
    for i in [1:1:length(pricebidsarr);]
        append!(pricebids, pricebidsarr[i])
    end
    return sort!(unique(pricebids / un))
end



##########################################################################
## returns the value of the FOC
##########################################################################
function FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
    WPar = W[group][bootstraprun][sort(
        findall(in.(prices, (bid.pb,))),
        rev = true,
    )]
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
            exp(
                rho * (
                    IntV(x, bid.cumqb[bidstep], v) -
                    IntBid(x, bid.cumqb[bidstep], bid)
                ),
            ) * (
                (StepV(x, v) - StepBid(x, bid)) * (
                    w(Q - x, [WParPlus[bidstep], WPar[bidstep]], dp) +
                    rho *
                    (x - bid.cumqb[bidstep-1]) *
                    cdf(WPar[NoStepBid(x, bid)], Q - x)
                ) - cdf(WPar[NoStepBid(x, bid)], Q - x)
            )

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
        return sum(val) +
               rho *
        (bid.cumqb[bidstep] - bid.cumqb[bidstep-1]) *
        PiOverline(bidstep, bid, v, WPar, rho, Q, n)
    else
        g(x) =
            exp(
                rho * (
                    IntV(x, bid.cumqb[bidstep], v) -
                    IntBid(x, bid.cumqb[bidstep], bid)
                ),
            ) * (
                (StepV(x, v) - StepBid(x, bid)) * (
                    w(Q - x, [WParPlus[bidstep], WPar[bidstep]], dp) +
                    rho * x * cdf(WPar[NoStepBid(x, bid)], Q - x)
                ) - cdf(WPar[NoStepBid(x, bid)], Q - x)
            )

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
        return sum(val) +
               rho *
        bid.cumqb[bidstep] *
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
        push!(vnew, [max.(v, vl[vl.qval.>q, :vval][1]), q])
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
    sort!(vnew, :qval)
end
