# This file contains the main functions for the estimation 
# of W, $\Theta$, and the different bounds.  

# Cf. Readme.md for more information about the respective functions.

function Wgamma(prices, auctionset, bidderassignment, n, m, P)

    # construct array of relevant bids, one array for every bidder group
    bidset = []
    for g in sort(unique(bidderassignment))
        bidsetgroup = []
        for j in auctionset
            for i in activebidderindeces[j][findall(in.(
                activebidderindeces[j],
                (bidderindeces[findall(in.(bidderassignment, g))],),
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

    #construct P resampling indeces, drawing [n(1),...,n(g)] bids, to be 
    #used on every bootstrap sample
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


    #for every group g, every bootstrap round, for every price, retrieve 
    #opponent demand P times and fit gamma distribution to it
    W = []
    for g in sort(unique(bidderassignment))
        dist = []
        for b in [1:1:m;]
            #retrieve relevant bids
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
            #for every price, return estimate of demand distribution
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

function Wlnorm(prices, auctionset, bidderassignment, n, m, P)
    # construct array of relevant bids, one array for every bidder group
    bidset = []
    for g in sort(unique(bidderassignment))
        bidsetgroup = []
        for j in auctionset
            for i in activebidderindeces[j][findall(in.(
                activebidderindeces[j],
                (bidderindeces[findall(in.(bidderassignment, g))],),
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

    #construct P resampling indeces, drawing [n(1),...,n(g)] bids, to be 
    #used on every bootstrap sample
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


    #for every group g, every bootstrap round, for every price, retrive 
    #opponent demand P times and fit lognormal distribution to it
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

function SimpleBound(bid, WPar, rho, Q)
    
    ## start with last quantity point
    vlb = DataFrame(vval = bid.pb[end], qval = bid.cumqb[end])
    ## no information in case of just one (p,q)-pair
    if length(bid.pb) == 1
        return [vlb, DataFrame(vval = vupperbar, qval = bid.cumqb[end])]
    end
    ## abort if WPar is not properly defined
    if (ismissing(WPar[length(bid.pb)-1]) || ismissing(WPar[length(bid.pb)]))
        return [
            DataFrame(vval = bid.pb, qval = bid.cumqb),
            DataFrame(
                vval = fill(vupperbar, length(bid.cumqb)),
                qval = bid.cumqb,
            ),
        ]
    else  
        ##using vlb on last step gives \PiOverline_i^j = 0
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

    ## move backwards through quantity points, if more than two in bid
    if length(bid.pb) >= 3
        for i in [length(bid.pb)-1:-1:2;]
            # abort early if WPar is not properly defined
            if (ismissing(WPar[i]) || ismissing(WPar[i+1]))
                vlbnew = bid.pb[1:i]
                qvalnew = bid.cumqb[1:i]
                append!(vlb, DataFrame(vval = vlbnew, qval = qvalnew))
                sort!(vlb, :qval)
                vubnew = fill(vupperbar, i)
                append!(vub, DataFrame(vval = vubnew, qval = qvalnew))
                sort!(vub, :qval)
                # construct decreasing bounds \overline{v}_i and \underline{v}_i;
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
            # abort early if WPar is not properly defined
            if (ismissing(WPar[i-1]) || ismissing(WPar[i]))
                qvalnew = bid.cumqb[1:i]
                vubnew = fill(vupperbar, i)
                append!(vub, DataFrame(vval = vubnew, qval = qvalnew))
                sort!(vub, :qval)
                # construct decreasing bounds \overline{v}_i and \underline{v}_i;
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

    ## lower bound at first quantity point
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
    ## upper bound at first quantity point is vupperbar
    push!(vub, [vupperbar bid.cumqb[1]])
    sort!(vub, :qval)

    # construct decreasing bounds \overline{v}_i and \underline{v}_i;
    if length(vlb.qval) >= 2
        for i in [1:1:length(vlb.qval);]
            vub[i, :vval] = maximum(vub[i:end, :vval])
            vlb[i, :vval] = max(bid.pb[i], minimum(vlb[1:i, :vval]))
        end
    end

    return [vlb, vub]
end

function EstimateSimpleBounds(auction, W, bidderassignment, prices, rhovec, m)
    bounds = []
    for bidder in activebidderindeces[auction]
        rho=rhovec[bidderassignment[findall(in.(bidderindeces,bidder))][1]]
        boundsbidder = []
        for bootstraprun in [1:1:m;]
            bid = qpBid(bidder, auction)
            g = bidderassignment[findall(in.(bidderindeces,bidder))][1]
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

function EstTheta(auction, W, bidderassignment, prices, bounds, rho, m)
    #define matrices, one colum for every bootrap round, one row for every bidder group
    tested = zeros(length(unique(bidderassignment)), m)  
    violated = zeros(length(unique(bidderassignment)), m) 


    for bidder in [1:1:activebidders[auction];]
        g = bidderassignment[findall(in.(bidderindeces,activebidderindeces[auction][bidder]))][1]
        bid = qpBid(activebidderindeces[auction][bidder], auction)
        for bootstraprun in [1:1:m;]
            for bidstep in [length(bid.pb):-1:1;]
                ## compute values of F^j for upper and lower bounds
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
                ## test for violation
                if (uc > 0 || lc < 0)
                    violated[g, bootstraprun] += 1
                end
                tested[g, bootstraprun] += 1
            end
        end
    end
    return [tested, violated]
end

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

    # initial conditions, convert initvl and initvu
    vl = copy(initvl)
    vu = copy(initvu)
    vl[!, :vval] = convert(Array{Float64,1}, vl[!, :vval])
    vl[!, :qval] = convert(Array{Float64,1}, vl[!, :qval])
    vu[!, :vval] = convert(Array{Float64,1}, vu[!, :vval])
    vu[!, :qval] = convert(Array{Float64,1}, vu[!, :qval])


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
    rhovec,
    m,
    maxiter,
    tolerance,
    )
    tighterbounds = []
    for bidder in [1:1:activebidders[auction];]
        rho=rhovec[bidderassignment[findall(in.(bidderindeces,activebidderindeces[auction][bidder]))][1]]
        boundsbidder = []
        #only compute tighter bounds if allocated with positive amount
        if qRec(activebidderindeces[auction][bidder], auction) > 0
            bid = qpBid(activebidderindeces[auction][bidder], auction)
            g = bidderassignment[findall(in.(bidderindeces,activebidderindeces[auction][bidder]))][1]
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
