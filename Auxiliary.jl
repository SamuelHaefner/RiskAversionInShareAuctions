# This file contains the auxiliary functions 
# and global variable used throughout the estimation

# Readme.md contains additional information.

# packages required
using CSV
using DataFrames
using Distributions
using BSON: @save, @load  
using Roots

### determine the global variables
#################################################

# read in data
bids = CSV.read(
    "setofbids.csv",
    DataFrame,
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

# upper bound on type space (number in cents, as are the bids)
vupperbar = 2053/un

# retreive the auction indeces of all auctions
# only use the first 39 auctions (the last auction has only one participant)
auctionindeces = unique(bids[!, :auction])[1:39]

# retreive the bidder indeces of all registered bidders
bidderindeces = unique(bids[!, :bidder])

# retrieve, for every auction, the quota, the market clearing price,
# the number of active bidders, and the indices of the active bidders
quotas = []
clearingprices = []
activebidders = []
activebidderindeces = []
for i in [1:1:length(auctionindeces);]
    push!(quotas, bids[bids.auction .== auctionindeces[i], :quotatot][1])
    push!(
        clearingprices,
        minimum(bids[(bids.auction .== auctionindeces[i]) .& (bids.qr .> 0), :pb]),
    )
    push!(
        activebidders,
        length(unique(bids[bids.auction .== auctionindeces[i], :bidder])),
    )
    push!(
        activebidderindeces,
        unique(bids[bids.auction .== auctionindeces[i], :bidder]),
    )
end

### define the auxiliary functions
#################################################

#################################################
# qpBid(bidder, auction)
#################################################
#### Description
# Retrieve the bid function of a bidder in an auction.
#### Arguments
# ```bidder``` -- bidderindex  
# ```auction``` -- auctionindex  
#### Return value
# A data frame: [qb (quantity points), pb (price points), cumqb (cumulated quantity points)].
#################################################

function qpBid(bidder, auction)
    # retreive price-quantity pairs
    bid = bids[
        (bids.bidder .== bidder) .& (bids.auction .== auctionindeces[auction]),
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

##################################################
# AvgNoBidders(bidderassignment)
##################################################
#### Description
# Determine the average number of bidders in each bidder cluster
#### Arguments
# ```bidderassignment``` -- vector, each element corresponding to a bidder, 
#                           denoting the cluster number of the respective bidder  
#### Return value
# A list of real numbers; one for each cluster number; returning the average no. of bidders.
##################################################

function AvgNoBidders(bidderassignment)
    avgno = []
    for g in sort(unique(bidderassignment))
        n = []
        for auction in [1:1:39;]
            push!(
                n,
                length(findall(in.(activebidderindeces[auction],
                    (findall(in.(bidderassignment, g)),),))),
            )
        end
        push!(avgno, ceil(mean(n)))
    end
    return convert(Array{Int64,1}, avgno)
end

##################################################
# BidToCover(auction)
##################################################
#### Description
# Determines the bid-to-cover ratio in an auction.
#### Arguments
# ```auction``` -- auctionindex
#### Return value
# Real number; corresponding to the bid-to-cover ratio.
##################################################

function BidToCover(auction)
    q = 0
    for bidder in activebidderindeces[auction]
        q += qpBid(bidder, auction).cumqb[end]
    end
    return q / quotas[auction]
end

##################################################
# Revenue(auction)
##################################################
#### Description
# Determines the revenue from an auction.
#### Arguments
# ```auction``` -- auctionindex
#### Return value
# Real number; corresponding to the auction revenue.
##################################################

function Revenue(auction)
    revenue = sum(bids[(bids.auction .== auctionindeces[auction]), :pr])
end

##################################################
# qRec(bidder, auction)
##################################################
#### Description
# Determines the allocated quantity of a bidder in a given auction.
#### Arguments
# ```bidder``` -- bidderindex   
# ```auction``` -- auctionindex
#### Return value
# Real number; corresponding to the allocated quantity.
##################################################

function qRec(bidder, auction)
    return sum(bids[
        (bids.bidder .== bidder) .& (bids.auction .== auctionindeces[auction]),
        :qr,
    ])
end

##################################################
# qShareRec(bidder, auction)
##################################################
#### Description
# Determines the share of quota that is allocated to a bidder in a given auction.
#### Arguments
# ```bidder``` -- bidderindex  
# ```auction``` -- auctionindex
#### Return value
# Real number; corresponding to the received share
##################################################

function qShareRec(bidder, auction)
    return sum(bids[
        (bids.bidder .== bidder) .& (bids.auction .== auctionindeces[auction]),
        :qperc,
    ])
end

##################################################
# ActiveAuctions(bidder)
##################################################
#### Description
# Determines the indeces of the auctions in which a bidder was active.
#### Arguments
# ```bidder``` -- bidderindex 
#### Return value
# A list of integers; corresponding to the indices of the auctions in 
# which the bidder was active.
##################################################

function ActiveAuctions(bidder)
    activeauctions = []
    for i in [1:1:39;]
        if in(bidder, activebidderindeces[i])
            push!(activeauctions, i)
        end
    end
    return activeauctions
end

##################################################
# AvgBid(bidder)
##################################################
#### Description
# Determines the average bid of a bidder across the auctions in which that bidder was active.
#### Arguments
# ```bidder``` -- bidderindex
#### Return value
# Real number; corresponding to the average bid of the bidder.
##################################################

function AvgBid(bidder)
    qbid = []
    for i in ActiveAuctions(bidder)
        push!(qbid, qpBid(bidder, i).cumqb[end])
    end
    return mean(qbid)
end

##################################################
# ShareSuccBidders(auction)
##################################################
#### Description
# Determines the  share of succesfull bidders (with non-zero allocated quantity) in an auction.
#### Arguments
# ```auction``` -- auctionindex
#### Return value
# Real number; corresponding to the share of successful bidders.
##################################################

function ShareSuccBidders(auction)
    succ = 0
    for i in activebidderindeces[auction]
        if qRec(i, auction) > 0
            succ += 1
        end
    end
    return succ / activebidders[auction]
end

##################################################
# SuccessRate(bidder)
##################################################
#### Description
# Returns the success rate of a bidder (succesfull participation/total participation).
#### Arguments
# ```bidder``` -- bidderindex
#### Return value
# Real number; corresponding to the bidder's success rate.
##################################################

function SuccessRate(bidder)
    succ = 0
    for i in ActiveAuctions(bidder)
        if qRec(bidder, i) > 0
            succ += 1
        end
    end
    return succ / length(ActiveAuctions(bidder))
end

##################################################
# StepBid(q, bid)
##################################################
#### Description
# Determines the value of $\beta_{bid}(q)$ (cf. the manuscript for a definition).
#### Arguments
# ```q``` -- positive real number
# ```bid``` -- bid function as returned from qpBid(bidder, auction)
#### Return value
# Real number; corresponding to the value of $\beta_{bid}(q)$.
###################################################

function StepBid(q, bid)
    if q <= minimum(bid[!, :cumqb])
        return bid[1, :pb]
    end
    if q > maximum(bid[!, :cumqb])
        return 0
    end
    return bid[length(bid[bid.cumqb .< q, :cumqb]) + 1, :pb]
end

###################################################
# NoStepBid(q, bid)
###################################################
#### Description
# Determines the step of the step function $\beta_b$ at ```q```. 
# That is, returns the number of downward jumps that have occured 
# strictly before ```q```, plus one.
#### Arguments
# ```q``` positive real number
# ```bid``` bid function as returned from qpBid(bidder, auction)
#### Return value
# Integer; see description.
###################################################

function NoStepBid(q, bid)
    if q <= minimum(bid[!, :cumqb])
        return 1
    end
    if q > maximum(bid[!, :cumqb])
        return length(bid.cumqb)
    end
    return length(bid[bid.cumqb .< q, :cumqb]) + 1
end

###################################################
# StepV(q, v)
###################################################
#### Description
# Returns the value of the profit function v(q).
#### Arguments
# ```q``` -- positive real number  
# ```v``` -- marginal profit function, 
#            which is a data frame with columns ```qval``` 
#            (quantities, need to be increasing) and 
#            ```vval``` (the values of v)
#### Return value
# Real number; corresponding to the value of the profit function v(q).
###################################################

function StepV(q, v)
    if q <= minimum(v[!, :qval])
        return v[1, :vval]
    end
    if q > maximum(v[!, :qval])
        return 0
    end
    return v[length(v[v.qval .< q, :qval]) + 1, :vval]
end

###################################################
# IntBid(a, b, bid)
###################################################
#### Description
# Returns the value of $\int_a^b\beta_{bid}(q)dq$.
#### Arguments
# ```a```, ```b``` -- positive real numbers
# ```bid``` -- bid function as returned from qpBid(bidder, auction)
#### Return value
# Real number; value of $\int_a^b\beta_{bid}(q)dq$.
###################################################

function IntBid(a, b, bid)
    x = bid[(bid.cumqb .> a) .& (bid.cumqb .< b), :cumqb]
    push!(x, b)
    pushfirst!(x, a)
    y = []
    y = [StepBid(z, bid) for z in x[2:end]]
    x = x[2:end] - x[1:end - 1]
    return round(x' * y; digits = 4)
end

###################################################
# IntV(a, b, v)
###################################################
#### Description
# Returns the value of $\int_a^b v(q)dq$.
#### Arguments
# ```a```, ```b``` -- positive real numbers
# ```v``` -- marginal profit function
#### Return value
# Real number; value of $\int_a^b v(q)dq$.
###################################################

function IntV(a, b, v)
    x = v[(v.qval .> a) .& (v.qval .< b), :qval]
    push!(x, b)
    pushfirst!(x, a)
    y = []
    y = [StepV(z, v) for z in x[2:end]]
    x = x[2:end] - x[1:end - 1]
    return round(x' * y; digits = 4)
end

###################################################
# PiOverline(bidstep, bid, v, WPar, rho, Q, n)
###################################################
#### Description
# Returns the value of $\overline{\Pi}_i^j(b,v)$ (cf. the manuscript for a definition).
#### Arguments
# ```bidstep``` -- positive natural number, corresponding to the number of the step under consideration, j  
# ```bid``` -- bid function  
#```v``` -- marginal profit function  
# ```WPar``` -- array, containing the estimated parameters of the distribution of D(p_i^j) for steps j=1,...,k  
# ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
# ```Q``` -- positive real number, corresponding to the quota $Q$  
# ```n``` -- positive natural number, indicating the number of support points used for integration  
#### Return value
# Real number; value of $\overline{\Pi}_i^j(b,v)$. 
###################################################

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
    ### in general, when function is called, use n=100. this gives a good approximation.
    intpts = v[v.qval .> bid.cumqb[bidstep], :qval]
    pushfirst!(intpts, bid.cumqb[bidstep])
    val = []
    for i in [1:1:length(intpts) - 1;]
        d = (intpts[i + 1] - intpts[i]) / n
        push!(val, d * sum([f(x) for x in [intpts[i] + d:d:intpts[i + 1];]]))
    end

    return try
        sum(val)
    catch
        0
    end
end

###################################################
# w(q, WPar, dp)
###################################################
#### Description
# Computes $w_i(p,q)$.
#### Arguments
# ```q``` -- positive real number  
# ```WPar``` -- array, containing two sets of estimated parameters 
#               of the distribution of D(p); one for p_i^j and one for p_i^j + dp.   
# ```dp``` -- positive real number  
#### Return value
# Real number; value of $w_i(p,q)$.
###################################################

function w(q, WPar, dp)
    return (cdf(WPar[1], q) - cdf(WPar[2], q)) / dp
end

####################################################
# PriceBids(auctionset)
####################################################
#### Description
# Returns all submitted prices in ```auctionset``` in ascending order.
#### Arguments
# ```auctionset``` -- vector, containing auction indeces
#### Return value
# List of real numbers; prices submitted in the auctions in ```auctionset```.
####################################################

function PriceBids(auctionset)
    pricebidsarr =
        [bids[(bids.auction .== auctionindeces[x]), :pb] for x in auctionset]
    pricebids = []
    for i in [1:1:length(pricebidsarr);]
        append!(pricebids, pricebidsarr[i])
    end
    return sort!(unique(pricebids / un))
end

####################################################
# FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
####################################################
#### Description
# ```FOC()``` returns the value of $F^j$ (cf. the manuscript).
#### Arguments
# ```bidstep``` -- positive natural number, corresponding to the number of the step under consideration, j  
# ```bid``` -- bid function  
# ```v``` -- marginal profit function  
# ```W``` -- estimates of W as returned from ```Wgamma()``` or ```Wlnorm()```  
# ```group``` -- natural number, indicating the group number of the auction  
# ```prices``` -- vector of submitted prices used for the estimation of W  
# ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
# ```Q``` -- positive real number, corresponding to the quota $Q$  
# ```bootstraprun``` -- natural number, indicating the boostrap run number   
# ```n``` -- positive natural number, indicating the number of support points used for integration  
#### Return value
# Real number; value of $F^j$.
###################################################

function FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
    # retrive estimates for W(p,q) at price points p (WPar) and at the respective 
    # next higher prices among the submitted bids (WParPlus)
    # in general, when function is called: n=100. this gives a good approximation.
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

    # dp is difference between price bid and next higher price bid among submitted bids 
    # for the current bidstep
    dp = try
        prices[findall(in.(prices, (bid.pb,))) .+ 1][bidstep] -
        prices[findall(in.(prices, (bid.pb,)))][bidstep]
    catch
        return 0
    end

    # make sure WPar and WParPlus is well defined
    if ismissing(WParPlus[bidstep]) || ismissing(WPar[bidstep])
        return 0
    end

    # need to distinguish between first and other steps (because in first step, 
    # integration starts at 0)
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
                    (x - bid.cumqb[bidstep - 1]) *
                    cdf(WPar[NoStepBid(x, bid)], Q - x)
                ) - cdf(WPar[NoStepBid(x, bid)], Q - x)
            )

        ### integrate f from bid.cumqb[bidstep-1] to bid.cumqb[bidstep]
        ### use intervals given with intpts
        intpts = v[
            (v.qval .> bid.cumqb[bidstep - 1]) .& (v.qval .< bid.cumqb[bidstep]),
            :qval,
        ]
        pushfirst!(intpts, bid.cumqb[bidstep - 1])
        push!(intpts, bid.cumqb[bidstep])
        val = []
        for i in [1:1:length(intpts) - 1;]
            d = (intpts[i + 1] - intpts[i]) / n
            push!(val, d * sum([f(x) for x in [intpts[i] + d:d:intpts[i + 1];]]))
        end

        ## return integral plus value of normalized interim utility beyond bidstep
        return sum(val) +
               rho *
        (bid.cumqb[bidstep] - bid.cumqb[bidstep - 1]) *
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

        ### integrate g from 0 to bid.cumqb[bidstep-1]
        ### use intervals given with intpts
        intpts = v[(v.qval .> 0) .& (v.qval .< bid.cumqb[bidstep]), :qval]
        pushfirst!(intpts, 0)
        push!(intpts, bid.cumqb[bidstep])
        val = []
        for i in [1:1:length(intpts) - 1;]
            d = (intpts[i + 1] - intpts[i]) / n
            push!(val, d * sum([g(x) for x in [intpts[i] + d:d:intpts[i + 1];]]))
        end

        ## return integral plus value of normalized interim utility beyond bidstep
        return sum(val) +
               rho *
        bid.cumqb[bidstep] *
        PiOverline(bidstep, bid, v, WPar, rho, Q, n)
    end
end

#################################################
# VarPhiU(q, v, vl)
#################################################
#### Description
# Determines $\varphi_u(q,v,v_l)$ (cf. the manuscript).
#### Arguments
# ```q``` -- positive real number  
# ```v``` -- marginal profit function  
# ```vl``` -- positive real number, corresponding to $v_l$  
#### Return value
# Real number; value of $\varphi_u(q,v,v_l)$.
#################################################

function VarPhiU(q, v, vl)
    vnew = DataFrame(
        vval = max.(fill(convert(Float64, v), length(vl[vl.qval .<= q, :qval])),
            vl[vl.qval .<= q, :vval],),
        qval = vl[vl.qval .<= q, :qval],
    )
    if !in(q, vl[:, :qval])
        push!(vnew, [max.(v, vl[vl.qval .> q, :vval][1]), q])
    end
    append!(vnew, vl[vl.qval .> q, :])
end

#################################################
# VarPhiL(q, v, vu)
#################################################
#### Description
# Determines $\varphi_l(q,v,v_u)$ (cf. the manuscript).
#### Arguments
# ```q``` -- positive real number  
# ```v``` -- marginal profit function  
# ```vu``` -- positive real number, corresponding to $v_u$  
#### Return value
# Real number; value of $\varphi_l(q,v,v_u)$.
#################################################


function VarPhiL(q, v, vu)
    vnew = DataFrame(
        vval = min.(fill(convert(Float64, v), length(vu[vu.qval .>= q, :qval])),
            vu[vu.qval .>= q, :vval],),
        qval = vu[vu.qval .>= q, :qval],
    )
    if !in(q, vu[:, :qval])
        push!(vnew, [vu[vu.qval .< q, :vval][end], q])
    end
    append!(vnew, vu[vu.qval .< q, :])
    sort!(vnew, :qval)
end
