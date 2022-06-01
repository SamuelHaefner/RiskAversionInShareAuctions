############################################################
# CheckDecreasing(bid, WPar, rho, Q)
############################################################
#### Description
# Computes upper and lower bounds on the rationalizable profit 
# functions at the submitted quantity points and checks whether there
# is a decreasing function rationalizing these bounds.
#### Arguments
# ```bid``` -- bid function  
# ```WPar``` -- array, containing the estimated parameters of W(p,q) at (p_i^j) for steps j=1,...,k  
# ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
# ```Q``` -- positive real number, corresponding to the quota $Q$
#### Return value
#  zero (if there is none) or one (if there is one)
############################################################

function CheckDecreasing(bid, WPar, rho, Q)
    
    ## start with last quantity point
    vlb = DataFrame(vval = bid.pb[end], qval = bid.cumqb[end])
    ## no information in case of just one (p,q)-pair
    if length(bid.pb) == 1
        return 1
    end
    ## abort if WPar is not properly defined
    if (ismissing(WPar[length(bid.pb)-1]) || ismissing(WPar[length(bid.pb)]))
        return 1
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
                # check condition: vul[i] \geq vub[i] for all i = 2,...,\ell_i
                for i in [1:1:length(vlb.qval);]
                    if vub.vval[i] < vlb.vval[i] return 0 end
                end
                return 1
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
                # check condition: vul[i] \geq vub[i] for all i = 2,...,\ell_i
                for i in [1:1:length(vlb.qval);]
                    if vub.vval[i] < vlb.vval[i] return 0 end
                end
                return 1
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

    # check condition: vul[i] \geq vub[i] for all i = 2,...,\ell_i
    for i in [1:1:length(vlb.qval);]
        if vub.vval[i] < vlb.vval[i] return 0 end
    end
    return 1
end

############################################################
# CheckDecreasingNew(bid, WPar, rho, Q)
############################################################
#### Description
# Computes upper and lower bounds on the rationalizable profit 
# functions at the submitted quantity points and checks whether there
# is a decreasing function rationalizing these bounds.
#### Arguments
# ```bid``` -- bid function  
# ```WPar``` -- array, containing the estimated parameters of W(p,q) at (p_i^j) for steps j=1,...,k  
# ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
# ```Q``` -- positive real number, corresponding to the quota $Q$
#### Return value
#  zero (if there is none) or one (if there is one)
############################################################

function CheckDecreasingNew(bid, WPar, rho, Q)
    
    ## start with last quantity point
    vlb = DataFrame(vval = bid.pb[end], qval = bid.cumqb[end])
    ## no information in case of just one (p,q)-pair
    if length(bid.pb) == 1
        return 1
    end
    ## abort if WPar is not properly defined
    if (ismissing(WPar[length(bid.pb)-1]) || ismissing(WPar[length(bid.pb)]))
        return 1
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
                        vub[i, :vval] = max(bid.pb[i], minimum(vub[1:i, :vval]))
                        vlb[i, :vval] = maximum(vlb[i:end, :vval])
                    end
                end
                # check condition: vul[i] \geq vub[i] for all i = 1,...,\ell_i
                for i in [1:1:length(vlb.qval);]
                    if vub.vval[i] < vlb.vval[i] return 0 end
                end
                return 1
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
                        vub[i, :vval] = max(bid.pb[i], minimum(vub[1:i, :vval]))
                        vlb[i, :vval] = maximum(vlb[i:end, :vval])
                    end
                end
                # check condition: vul[i] \geq vub[i] for all i = 1,...,\ell_i
                for i in [1:1:length(vlb.qval);]
                    if vub.vval[i] < vlb.vval[i] return 0 end
                end
                return 1
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
            vub[i, :vval] = max(bid.pb[i], minimum(vub[1:i, :vval]))
            vlb[i, :vval] = maximum(vlb[i:end, :vval])
        end
    end
    # check condition: vul[i] \geq vub[i] for all i = 1,...,\ell_i
    for i in [1:1:length(vlb.qval);]
        if vub.vval[i] < vlb.vval[i] return 0 end
    end
    return 1
end

###########################################################
# CheckDecreasing{Upper,Lower}(bid, WPar, rho, Q)
############################################################
#### Description
# Computes upper and lower bounds on the rationalizable profit 
# functions at the submitted quantity points and checks whether the 
# upper/lower bounds are decreasing
#### Arguments
# ```bid``` -- bid function  
# ```WPar``` -- array, containing the estimated parameters of W(p,q) at (p_i^j) for steps j=1,...,k  
# ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
# ```Q``` -- positive real number, corresponding to the quota $Q$
#### Return value
#  zero (if there is none) or one (if there is one)
############################################################

function CheckDecreasingUpper(bid, WPar, rho, Q)
    
    ## start with last quantity point
    vlb = DataFrame(vval = bid.pb[end], qval = bid.cumqb[end])
    ## no information in case of just one (p,q)-pair
    if length(bid.pb) == 1
        return 1
    end
    ## abort if WPar is not properly defined
    if (ismissing(WPar[length(bid.pb)-1]) || ismissing(WPar[length(bid.pb)]))
        return 1
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
                # check if decreasing: vub[i-1] \geq vub[i] for all i = 2,...,\ell_i
                for i in [1:1:length(vub.qval)-1;]
                    if vub.vval[i] < vub.vval[i+1] return 0 end
                end
                return 1
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
                # check if decreasing: vub[i-1] \geq vub[i] for all i = 2,...,\ell_i
                for i in [1:1:length(vub.qval)-1;]
                    if vub.vval[i] < vub.vval[i+1] return 0 end
                end
                return 1
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

    # check if decreasing: vub[i-1] \geq vub[i] for all i = 2,...,\ell_i
    for i in [1:1:length(vub.qval)-1;]
        if vub.vval[i] < vub.vval[i+1] return 0 end
    end
    return 1
end

function CheckDecreasingLower(bid, WPar, rho, Q)
    
    ## start with last quantity point
    vlb = DataFrame(vval = bid.pb[end], qval = bid.cumqb[end])
    ## no information in case of just one (p,q)-pair
    if length(bid.pb) == 1
        return 1
    end
    ## abort if WPar is not properly defined
    if (ismissing(WPar[length(bid.pb)-1]) || ismissing(WPar[length(bid.pb)]))
        return 1
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
                # check if decreasing: vlb[i-1] \geq vlb[i] for all i = 2,...,\ell_i
                for i in [1:1:length(vlb.qval)-1;]
                    if vlb.vval[i] < vlb.vval[i+1] return 0 end
                end
                return 1
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
                # check if decreasing: vlb[i-1] \geq vlb[i] for all i = 2,...,\ell_i
                for i in [1:1:length(vlb.qval)-1;]
                    if vlb.vval[i] < vlb.vval[i+1] return 0 end
                end
                return 1
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

    # check if decreasing: vlb[i-1] \geq vlb[i] for all i = 2,...,\ell_i
    for i in [1:1:length(vlb.qval)-1;]
        if vlb.vval[i] < vlb.vval[i+1] return 0 end
    end
    return 1
end


############################################################
# CheckSimpleBoundsDecreasing(auction, W, bidderassignment,  prices, rhovec, m)
############################################################
#### Description
# Computes ```CheckDecreasing()``` for all bidders in ```auction``` (for first bootstrap estimate).
#### Arguments
# ```auction``` -- auction index  
# ```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
# ```bidderassignment``` -- bidder assignment vector  
# ```prices``` -- vector of submitted prices used for the estimation of W  
# ```rhovec``` -- vector of positive real number, each number corresponding to the risk preference $\rho_g$ in bidder group $g$  
# ```m``` -- number of bootstrap rounds to be estimated
#### Return value
# A vector of zeros or ones, one element per bidder in the auction
#############################################################

function CheckSimpleBoundsDecreasing(auction, W, bidderassignment, prices, rhovec)
    check = zeros(length(activebidderindeces[auction]))  
    i=1
    for bidder in activebidderindeces[auction]
        bid = qpBid(bidder, auction)
        rho=rhovec[bidderassignment[findall(in.(bidderindeces,bidder))][1]]
        bootstraprun=1
        g = bidderassignment[findall(in.(bidderindeces,bidder))][1]
        WPar = W[g][bootstraprun][sort(
                findall(in.(prices, (bid.pb,))),
                rev = true,
            )]
        check[i] = CheckDecreasing(bid, WPar, rho, quotas[auction])
        i+=1
    end
    return check
end

function CheckSimpleBoundsAltDecreasing(auction, W, bidderassignment, prices, rhovec)
    check = zeros(length(activebidderindeces[auction]))  
    i=1
    for bidder in activebidderindeces[auction]
        bid = qpBid(bidder, auction)
        rho=rhovec[bidderassignment[findall(in.(bidderindeces,bidder))][1]]
        bootstraprun=1
        g = bidderassignment[findall(in.(bidderindeces,bidder))][1]
        WPar = W[g][bootstraprun][sort(
                findall(in.(prices, (bid.pb,))),
                rev = true,
            )]
        check[i] = CheckDecreasingAlt(bid, WPar, rho, quotas[auction])
        i+=1
    end
    return check
end

############################################################
# CheckSimpleBoundsDecreasing{Upper,Lower}(auction, W, bidderassignment,  prices, rhovec, m)
############################################################
#### Description
# Computes ```CheckDecreasing{Upper,Lower}()``` for all bidders in ```auction``` 
# (for first bootstrap estimate).
#### Arguments
# ```auction``` -- auction index  
# ```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
# ```bidderassignment``` -- bidder assignment vector  
# ```prices``` -- vector of submitted prices used for the estimation of W  
# ```rhovec``` -- vector of positive real number, each number corresponding to the risk preference $\rho_g$ in bidder group $g$  
# ```m``` -- number of bootstrap rounds to be estimated
#### Return value
# A vector of zeros or ones, one element per bidder in the auction
#############################################################

function CheckSimpleBoundsDecreasingUpper(auction, W, bidderassignment, prices, rhovec)
    check = zeros(length(activebidderindeces[auction]))  
    i=1
    for bidder in activebidderindeces[auction]
        bid = qpBid(bidder, auction)
        rho=rhovec[bidderassignment[findall(in.(bidderindeces,bidder))][1]]
        bootstraprun=1
        g = bidderassignment[findall(in.(bidderindeces,bidder))][1]
        WPar = W[g][bootstraprun][sort(
                findall(in.(prices, (bid.pb,))),
                rev = true,
            )]
        check[i] = CheckDecreasingUpper(bid, WPar, rho, quotas[auction])
        i+=1
    end
    return check
end

function CheckSimpleBoundsDecreasingLower(auction, W, bidderassignment, prices, rhovec)
    check = zeros(length(activebidderindeces[auction]))  
    i=1
    for bidder in activebidderindeces[auction]
        bid = qpBid(bidder, auction)
        rho=rhovec[bidderassignment[findall(in.(bidderindeces,bidder))][1]]
        bootstraprun=1
        g = bidderassignment[findall(in.(bidderindeces,bidder))][1]
        WPar = W[g][bootstraprun][sort(
                findall(in.(prices, (bid.pb,))),
                rev = true,
            )]
        check[i] = CheckDecreasingLower(bid, WPar, rho, quotas[auction])
        i+=1
    end
    return check
end
