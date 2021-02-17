# This file  contains the main functions to test for the monotonicity of F in v.

# Cf. Readme.md for more information on the respective functions.

function TestIncreasingDiff(auction, bidderset, BoundsAuction, W, R, rhoindexset, rhovec)
    T = []
    for rhoindex in rhoindexset
        Trho = []
        for bidder in bidderset
            try
                push!(
                    Trho,
                    TestIncreasingDiffBidder(
                        auction,
                        bidder,
                        BoundsAuction[rhoindex],
                        W,
                        R,
                        rhovec[rhoindex],
                    ),
                )
            catch
                push!(Trho, [-1, -1])
            end
        end
        push!(T, Trho)
    end
    return T
end

function GetQVals(g, bidstep)
    # random number of steps
    steps = rand([1:1:20;], 1)[1]
    # random q values between steps
    if bidstep == 1
        qval = [rand([1:1:g.qval[bidstep]-1;], 1)[1], g.qval[bidstep]]
    else
        qval = [g.qval[bidstep-1], g.qval[bidstep]]
    end
    append!(qval, rand([qval[1]:1:qval[2];], steps))
    qval = sort(unique(qval))
end

function GetHigherF(g, h, qval)
    vval = []
    push!(vval, rand([StepV(qval[1], g):0.01:StepV(qval[1], h);], 1)[1])
    for i in [2:1:length(qval);]
        push!(
            vval,
            rand([
                StepV(qval[i], g):0.01:minimum([StepV(qval[i], h), vval[i-1]]);
            ])[1],
        )
    end
    f1 = DataFrame([vval, qval])

    rename!(f1, (:x1 => :vval))
    rename!(f1, (:x2 => :qval))
end


function CheckFOCMonotone(
    bidstep,
    bid,
    vub,
    vlb,
    W,
    g,
    prices,
    rho,
    Q,
    bootstraprun,
    R,
    )
    
    FOCdiff = []
    for i in [1:1:R;]
        qval = GetQVals(vlb, bidstep)
        fun1 = GetHigherF(vlb, vub, qval)
        fun2 = GetHigherF(fun1, vub, qval)

        fun1 = [fun1; vlb[vlb.qval.>vlb.qval[bidstep], [:vval, :qval]]]
        fun2 = [fun2; vlb[vlb.qval.>vlb.qval[bidstep], [:vval, :qval]]]

        push!(
            FOCdiff,
            FOC(bidstep, bid, fun2, W, g, prices, rho, Q, bootstraprun, 10) -
            FOC(bidstep, bid, fun1, W, g, prices, rho, Q, bootstraprun, 10),
        )
    end
    return FOCdiff
end

function LoadWandBounds(auction)
    eval(Meta.parse(string(
        "@load \"BoundsTestMon",
        auction,
        ".dat\" BoundsAuction",
    )))
    @load "WValuesTestMon.dat" W
    return [BoundsAuction, W]
end

function TestIncreasingDiffBidder(auction, bidder, Bounds, W, R, rho)

    g = findall([in(auction, group[x]) for x in [1:1:length(group);]])[1]
    
    auctionset = group[g]
    prices = PriceBids(auctionset)

    test = 0
    stepviol = 0
    
    bid = qpBid(activebidderindeces[auction][bidder], auction)
    
    biddergroup = bidderassignment[activebidderindeces[auction][bidder]]
    vubr = Bounds[bidder][1][2]
    vlbr = Bounds[bidder][1][1]
   
    for bidstep in [1:1:length(bid.pb);]
        teststep = copy(test)
        test +=
            length(findall(
                CheckFOCMonotone(
                    bidstep,
                    bid,
                    vubr,
                    vlbr,
                    W[g],
                    biddergroup,
                    prices,
                    rho,
                    quotas[auction],
                    1,
                    R,
                ) .< 0,
            )) / R
        if test > teststep
            stepviol += 1
        end
    end
    return ([test, stepviol])
end

