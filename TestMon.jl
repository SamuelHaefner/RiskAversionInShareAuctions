# This file  contains the main functions to test for the monotonicity of F in v.

# For a description of the procedure, see the manuscript.
# For more information about how the functions below relate to each other,
# see the Readme.md (especially Comment 1).


#######################################################
# TestIncreasingDiff(auction, bidderset, BoundsAuction, W, R, rhoindexset, rhovec)
#######################################################  
#### Description
# Determines the number of price-quantity pairs for which 
# the monotonicity assumption on $F^j$ is violated for a given auction.
#### Arguments
# ```auction``` -- auction index  
# ```bidderset``` -- bidder index set (set of bidders to be looked at)  
# ```BoundsAuction``` -- first element in list returned by ```LoadWandBounds()```  
# ```W``` -- second element in list returned by ```LoadWandBounds()```  
# ```R``` -- number of tests per bid step  
# ```rhoindexset``` -- vector of indexes, referring to rhovec, for which the test is run  
# ```rhovec``` -- vector of rho values  
#### Return value
# A list of length of the ```rhoindexset```, each entry is again a list 
# of length of the ```bidderset```, and each entry of that list consists of 
# two numbers: (1) the total number of price-quantity pairs tested for 
# that bidder (2) the total number of price-quantity pairs for which a 
# violation of monotonicity was found.
#######################################################

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


########################################################
# GetQVals(g, bidstep)
########################################################  
#### Description
# Draws a random number (uniformly between 1 and 20) of equidistant 
# x-values that lie between the x-values of the ```bidstep```-th 
# price-quantity pair and the (```bidstep```-1)-th price-quantity pair in ```g```.
#### Arguments
# ```g``` -- decreasing step function function with a number of steps of at least *bidstep*  
# ```bidstep``` -- natural number, indicating the number of the step in the bid function  
#### Return value
# List of x-values.
########################################################

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


########################################################
# GetHigherF(g, h, qval)
########################################################
#### Description
# Compute a random step function between ```g``` and ```h``` with x values given in ```qval```
#### Arguments
# ```g```, ```h``` -- decreasing step function, satisfying ```g``` > ```h```    
# ```qvals``` -- x-values returned by ```GetQVals()```
#### Return value
# Data frame containing the step function, [vval,qval].
########################################################

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


#########################################################
# CheckFOCMonotone(bidstep, bid, vub, vlb, W, g, prices, rho, Q, bootstraprun, R)
#########################################################
#### Description
# Evaluates $F^j$ at the two different, ordered profit functions ```R``` 
# times and reports the differences.
#### Arguments
# ```bidstep``` -- natural number, the bid step to be looked at  
# ```bid``` -- bid function  
# ```vub``` -- marginal profit function, (simple) upper bound  
# ```vlb``` -- marginal profit function, (simple) lower bound  
# ```W``` -- estimate of W as obtained by ```Wgamma()``` or ```Wlnorm()```
# ```g``` -- natural number, indicating the group number of the auction  
# ```prices``` -- vector of prices used for the estimation of W  
# ```rho``` -- real number, risk preference $\rho$  
# ```Q``` -- real number, the quota $Q$  
# ```bootstraprun``` -- the number of the bootrap round  
# ```R``` -- the number of tests per bid step  
#### Return value
# List of length ```R```, containing the differences.
##########################################################

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


##########################################################
# LoadWandBounds(auction)
##########################################################
#### Description
# Loads estimates of W and of simple bounds for ```auction``` from corresponding .dat files.
#### Arguments
# ```auction``` -- auction index  
#### Return value
# List containing W and simple bounds for ```auction```
##########################################################

function LoadWandBounds(auction)
    eval(Meta.parse(string(
        "@load \"BoundsTestMon",
        auction,
        ".dat\" BoundsAuction",
    )))
    @load "WValuesTestMon.dat" W
    return [BoundsAuction, W]
end


##########################################################
# TestIncreasingDiffBidder(auction, bidder, Bounds, W, R, rho)
##########################################################
#### Description
# Determines for a bidder in a given auction whether monotonicity 
# of $F^j$ is violated for the submitted price-quantity pairs. 
#### Arguments
# ```auction``` -- auction index  
# ```bidder``` -- bidder index  
# ```Bounds``` -- first element in list returned by ```LoadWandBounds(```)  
# ```W``` -- second element in list returned by ```LoadWandBounds()```  
# ```R``` -- number of tests per bid step  
# ```rho``` -- real number, risk preference $\rho$  
#### Return value
# Two-dimensional list, 
# [number of tests performed (corresponds to the number of bid steps), number of violations]
##########################################################

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

