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
        bidsetgroup = []
        for j in auctionset
            for i in activebidderindeces[j][findall(in.(
                activebidderindeces[j],
                (findall(in.(bidderassignment, g)),),
            ))]
                push!(bidsetgroup, qpBid(i, j))
            end
        end
        push!(bidset, bidsetgroup)
    end

    #construct m bootstrap indeces
    bootstrapped = []
    for i in [1:1:m;]
        push!(
            bootstrapped,
            [
                rand([1:1:length(bidset[x]);], length(bidset[x]))
                for x in sort(unique(bidderassignment))
            ],
        )
    end

    #construct P resampling indeces, drawing [n(1),...,n(g)] bids, to be used on every
    #bootstrap sample
    resampling = []
    for g in sort(unique(bidderassignment))
        resamplinggroup = []
        nn = []
        for h in sort(unique(bidderassignment))
            if h == g
                push!(nn, n[h] - 1)
            else
                push!(nn, n[h])
            end
        end
        for i in [1:1:P;]
            push!(
                resamplinggroup,
                [
                    rand([1:1:length(bidset[x]);], nn[x])
                    for x in sort(unique(bidderassignment))
                ],
            )
        end
        push!(resampling, resamplinggroup)
    end


    #for every group g, every bootstrap round, for every price,  retrive P opponent demand and
    #fit gamma distribution to it
    W = []
    for g in sort(unique(bidderassignment))
        dist = []
        for b in [1:1:m;]
            #retrieve bids
            relbids = []
            for r in [1:1:P;]
                samplebid = []
                for h in [1:1:length(unique(bidderassignment));]
                    append!(
                        samplebid,
                        bidset[h][bootstrapped[b][h]][resampling[g][r][h]],
                    )
                end
                push!(relbids, samplebid)
            end
            #for every price, return estimate of demand
            z = []
            for p in prices
                y = []
                for r in [1:1:P;]
                    push!(
                        y,
                        sum([
                            sum(relbids[r][x][relbids[r][x].pb.>=p, :qb])
                            for x in [1:1:length(relbids[r]);]
                        ]),
                    )
                end
                push!(z, try
                    fit(Gamma, convert(Array{Float64,1}, y))
                catch
                    missing
                end)
            end
            push!(dist, z)
        end
        push!(W, dist)
    end
    return (W)
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
        bidsetgroup = []
        for j in auctionset
            for i in activebidderindeces[j][findall(in.(
                activebidderindeces[j],
                (findall(in.(bidderassignment, g)),),
            ))]
                push!(bidsetgroup, qpBid(i, j))
            end
        end
        push!(bidset, bidsetgroup)
    end

    #construct m bootstrap indeces
    bootstrapped = []
    for i in [1:1:m;]
        push!(
            bootstrapped,
            [
                rand([1:1:length(bidset[x]);], length(bidset[x]))
                for x in sort(unique(bidderassignment))
            ],
        )
    end

    #construct P resampling indeces, drawing [n(1),...,n(g)] bids, to be used on every
    #bootstrap sample
    resampling = []
    for g in sort(unique(bidderassignment))
        resamplinggroup = []
        nn = []
        for h in sort(unique(bidderassignment))
            if h == g
                push!(nn, n[h] - 1)
            else
                push!(nn, n[h])
            end
        end
        for i in [1:1:P;]
            push!(
                resamplinggroup,
                [
                    rand([1:1:length(bidset[x]);], nn[x])
                    for x in sort(unique(bidderassignment))
                ],
            )
        end
        push!(resampling, resamplinggroup)
    end


    #for every group g, every bootstrap round, for every price,  retrive P opponent demand and
    #fit gamma distribution to it
    W = []
    for g in sort(unique(bidderassignment))
        dist = []
        for b in [1:1:m;]
            #retrieve bids
            relbids = []
            for r in [1:1:P;]
                samplebid = []
                for h in [1:1:length(unique(bidderassignment));]
                    append!(
                        samplebid,
                        bidset[h][bootstrapped[b][h]][resampling[g][r][h]],
                    )
                end
                push!(relbids, samplebid)
            end
            #for every price, return estimate of demand
            z = []
            for p in prices
                y = []
                for r in [1:1:P;]
                    push!(
                        y,
                        sum([
                            sum(relbids[r][x][relbids[r][x].pb.>=p, :qb])
                            for x in [1:1:length(relbids[r]);]
                        ]),
                    )
                end
                push!(z, try
                    fit(LogNormal, convert(Array{Float64,1}, y))
                catch
                    missing
                end)
            end
            push!(dist, z)
        end
        push!(W, dist)
    end
    return (W)
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
        return [vlb, DataFrame(vval = vupperbar, qval = bid.cumqb[end])]
    end

    ## using vlb on last step gives \PiOverline_i^j = 0
    if (ismissing(WPar[length(bid.pb)-1]) || ismissing(WPar[length(bid.pb)]))
        return [
            DataFrame(vval = bid.pb, qval = bid.cumqb),
            DataFrame(
                vval = fill(vupperbar, length(bid.cumqb)),
                qval = bid.cumqb,
            ),
        ]
    else
        vub = DataFrame(
            vval = max(
                bid.pb[length(bid.pb)-1],
                min(
                    vupperbar,
                    bid.pb[length(bid.pb)-1] +
                    (bid.pb[length(bid.pb)-1] - bid.pb[length(bid.pb)]) *
                    cdf(WPar[length(bid.pb)], Q - bid.cumqb[length(bid.pb)-1]) / (
                        cdf(
                            WPar[length(bid.pb)-1],
                            Q - bid.cumqb[length(bid.pb)-1],
                        ) - cdf(
                            WPar[length(bid.pb)],
                            Q - bid.cumqb[length(bid.pb)-1],
                        )
                    ),
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
                append!(vlb, DataFrame(vval = vlbnew, qval = qvalnew))
                sort!(vlb, :qval)
                vubnew = fill(vupperbar, i)
                append!(vub, DataFrame(vval = vubnew, qval = qvalnew))
                sort!(vub, :qval)
                # construct decreasing bounds \overline{v}_i and \underline{v}_i;
                # take most conservative bound from estimated vub and vul
                if length(vlb.qval) >= 2
                    for i in [1:1:length(vlb.qval);]
                        vub[i, :vval] = maximum(vub[i:end, :vval])
                        vlb[i, :vval] = max(bid.pb[i], minimum(vlb[1:i, :vval]))
                    end
                end
                return [vlb, vub]
            else
                vlbnew = min(
                    vupperbar,
                    max(
                        bid.pb[i],
                        bid.pb[i] +
                        (bid.pb[i] - bid.pb[i+1]) * (
                            cdf(WPar[i+1], Q - bid.cumqb[i]) -
                            rho * PiOverline(i, bid, vub, WPar, rho, Q, 100)[1]
                        ) / (
                            cdf(WPar[i], Q - bid.cumqb[i]) -
                            cdf(WPar[i+1], Q - bid.cumqb[i])
                        ),
                    ),
                )
            end
            push!(vlb, [vlbnew bid.cumqb[i]])
            sort!(vlb, :qval)
            # check whether distribution estimate produced missing,
            # if so, return early
            if (ismissing(WPar[i-1]) || ismissing(WPar[i]))
                qvalnew = bid.cumqb[1:i]
                vubnew = fill(vupperbar, i)
                append!(vub, DataFrame(vval = vubnew, qval = qvalnew))
                sort!(vub, :qval)
                # construct decreasing bounds \overline{v}_i and \underline{v}_i;
                # take most conservative bound from estimated vub and vul
                if length(vlb.qval) >= 2
                    for i in [1:1:length(vlb.qval);]
                        vub[i, :vval] = maximum(vub[i:end, :vval])
                        vlb[i, :vval] = max(bid.pb[i], minimum(vlb[1:i, :vval]))
                    end
                end
                return [vlb, vub]
            else
                vubnew = max(
                    bid.pb[i-1],
                    min(
                        vupperbar,
                        bid.pb[i-1] +
                        (bid.pb[i-1] - bid.pb[i]) * (
                            cdf(WPar[i], Q - bid.cumqb[i-1]) -
                            rho *
                            PiOverline(i - 1, bid, vlb, WPar, rho, Q, 100)[1]
                        ) / (
                            cdf(WPar[i-1], Q - bid.cumqb[i-1]) -
                            cdf(WPar[i], Q - bid.cumqb[i-1])
                        ),
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
                (bid.pb[1] - bid.pb[2]) * (
                    cdf(WPar[2], Q - bid.cumqb[1]) -
                    rho * PiOverline(1, bid, vub, WPar, rho, Q, 100)[1]
                ) / (
                    cdf(WPar[1], Q - bid.cumqb[1]) -
                    cdf(WPar[2], Q - bid.cumqb[1])
                ),
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
            vlb[i, :vval] = max(bid.pb[i], minimum(vlb[1:i, :vval]))
        end
    end

    return [vlb, vub]
end


###########################################################################
## Estimates simple bounds for all bidders in [auction] using bids from
## [auctionset] assuming a risk preference rho. The resampling algorithm
## uses [P] redraws on [m] bootstrap estimates. The function returns the
## price bids, the corresponding W, as well as the bounds for every bidder
## every bootstrap run.
###########################################################################
function EstimateSimpleBounds(auction, W, bidderassignments, prices, rho, m)
    bounds = []
    for bidder in activebidderindeces[auction]
        boundsbidder = []
        for bootstraprun in [1:1:m;]
            bid = qpBid(bidder, auction)
            g = bidderassignment[bidder]
            WPar = W[g][bootstraprun][sort(
                findall(in.(prices, (bid.pb,))),
                rev = true,
            )]
            boundsbidderbootstraprun =
                SimpleBound(bid, WPar, rho, quotas[auction])
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
    tested = zeros(length(unique(bidderassignment)), m)  #define matrix, one colum for every bootrap round, one row for every bidder group
    violated = zeros(length(unique(bidderassignment)), m) #define matrix, one column for every bootrap round, one row for every bidder group
    for bidder in [1:1:activebidders[auction];]
        g = bidderassignment[bidder]
        bid = qpBid(activebidderindeces[auction][bidder], auction)
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
                    100,
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
                    100,
                )
                if (uc > 0 || lc < 0)
                    violated[g, bootstraprun] += 1
                end
                tested[g, bootstraprun] += 1
            end
        end
    end
    return [tested, violated]
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
    vl[!, :vval] = convert(Array{Float64,1}, vl[!, :vval])
    vl[!, :qval] = convert(Array{Float64,1}, vl[!, :qval])
    vu[!, :vval] = convert(Array{Float64,1}, vu[!, :vval])
    vu[!, :qval] = convert(Array{Float64,1}, vu[!, :qval])


    # partition step into 10 subintervals
    for bidstep in [length(bid.pb):-1:1;]
        print(string("b", bidstep))
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
            print(string("r", repetition))
            # obtain new values for vu
            vvals = Array{Float64}(undef, 0)
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
                    if fl(StepV(q, initvu)) < 0 || fl(StepV(q, initvl)) > 0
                        push!(vvals, StepV(q, initvl))
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

function EstTighterBounds(
    auction,
    W,
    bidderassignment,
    prices,
    bounds,
    rho,
    m,
    maxiter,
    tolerance,
)
    tighterbounds = []
    for bidder in [1:1:activebidders[auction];]
        #for bidder in [1:1:2;]
        print(bidder)
        boundsbidder = []
        #only compute tighter bounds if allocated with positive amount
        if qRec(activebidderindeces[auction][bidder], auction) > 0
            bid = qpBid(activebidderindeces[auction][bidder], auction)
            g = bidderassignment[bidder]
            for bootstraprun in [1:1:m;]
                initvl = bounds[bidder][bootstraprun][1]
                initvu = bounds[bidder][bootstraprun][2]
                boundsbidderbootstraprun = TighterBounds(
                    bid,
                    initvl,
                    initvu,
                    W,
                    g,
                    prices,
                    rho,
                    quotas[auction],
                    bootstraprun,
                    maxiter,
                    tolerance,
                )
                push!(boundsbidder, boundsbidderbootstraprun)
            end
        else
            for bootstraprun in [1:1:m;]
                initvl = bounds[bidder][bootstraprun][1]
                initvu = bounds[bidder][bootstraprun][2]
                push!(boundsbidder, [initvl, initvu])
            end
        end
        push!(tighterbounds, boundsbidder)
    end
    return [bounds, tighterbounds]
end
