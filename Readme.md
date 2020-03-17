This files contains information about the Julia scripts and functions used in:

#### **Risk Aversion in Share Auctions: Estimating Rents from TRQs in Switzerland**

Author: Samuel Häfner, University of St. Gallen, [samuel.haefner@unisg.ch](mailto:samuel.haefner@unisg.ch) 

------

The julia files contain both the functions for estimation as well as the scripts to read and process the obtained estimates.

The main files are the following:

- **Auxiliary.jl** --- reads in the required packages, the data set and defines the required auxiliary functions and  global variables.
- **Grouping.jl** --- determines the bidder groups and auction groups used for the estimation.
- **Estimation.jl** --- contains the main functions for the estimation of W and the bounds.
- **TestMon.jl** --- contains the main functions to test for the monotonicity of F in v.

The following files contain the scripts used for estimation:

- **EstimateWandTGeneric.jl** --- generic script to compute and save estimates of W and Theta.
- **EstimateTighterBoundsGeneric.jl** --- generic script to compute standard and tighter bounds, using W from EstimateWandTGeneric.jl. 
- **TestMonGeneric.jl** --- generic script to test for monotonicity of F, determining W and the bounds for one bootstrap round, computing violations of monotonicity.

The following files contain the scripts to read in and process the obtained estimates:

- **ReadEstimates.jl** --- contains the functions to read in the estimates.
- **ReadT.jl** --- script to read in estimates of T and produce the plots. 
- **ReadTestMon.jl** --- script to read in the monotonicity violations obtained and saved in TestMonGeneric.jl.

The last file contains scripts used for the plots:

- **ResamplingPlots.jl** --- scripts to generate the group plots and the resampling plots. 

------

**Auxiliary.jl** --- reads in the required packages, the data set and defines a few important global variables.

*Global variables:*

- bids: a data frame containing the bidding data.
- K: max number of price-quantity pairs.
- un: unit of account (100 corresponds to CHF/kg).
- vupperbar: upper bound on type space (in CHF).
- auctionindeces: vector containing the auctionindeces corresponding to those in the data.
- bidderindeces: vector containing the bidderindeces corresponding to those in the data.
- quotas: vector with the quotas.
- clearingprices: vector with the market clearing prices.
- activebidders: vector with the number of active bidders.
- activebidderindeces: list of vectors with the indeces of the active bidders.

*Auxiliary functions:*

- qpBid(bidder,auction)
- AvgNoBidders(bidderassignment)
- BidToCover(auction)
- Revenue(auction)
- qRec(bidder,auction)
- qShareRec(bidder,auction)
- ActiveAuctions(bidder)
- AvgBid(bidder)
- ShareSuccBidders(auction)
- SuccessRate(bidder)
- StepBid(q, bid)
- NoStepBid(q, bid)
- StepV(q, v)
- IntBid(a, b, bid)
- IntV(a, b, v)
- PiOverline(bidstep, bid, v, WPar, rho, Q, n)
- w(q, WPar, dp)
- PriceBids(auctionset)
- FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
- VarPhiU(q, v, vl)
- VarPhiL(q, v, vu)


------

**Grouping.jl** --- determines the bidder groups and auction groups used for the estimation

- group: list containing, for each group, a vector of the respective auction indeces 
- bidderassignment: a vector determining for every bidder the group assignment {1,2,3}

------

**Estimation.jl** --- contains the main functions for estimation 

*Estimation of W:*

- Wgamma(prices, auctionset, bidderassignment, n, m, P)
- Wlnorm(prices, auctionset, bidderassignment, n, m, P)

Estimation of Simple Bounds:

- SimpleBound(bid, WPar, rho, Q)
- EstimateSimpleBounds(auction, W, bidderassignments,  prices, rho, m)

Estimation of Theta:

- EstTheta(auction, W, bidderassignment, prices, bounds, rho, m)

Estimation of Tighter Bounds:

- TighterBounds(bid, initvl, initvu, W, g, prices, rho, Q, bootstraprun, maxiter, tolerance)

- EstTighterBounds(auction, W, bidderassignment, prices, bounds, rho, m, maxiter, tolerance)

  ​

------

**TestMon.jl** --- contains the main functions to test for the monotonicity of F in v



------

