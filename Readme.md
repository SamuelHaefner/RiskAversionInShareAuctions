This file contains information about the Julia code used in:

# **[Risk Aversion in Share Auctions: Estimating Rents from TRQs in Switzerland](https://dx.doi.org/10.2139/ssrn.3397027)**

Author: [Samuel Häfner](https://samuelhaefner.github.io), University of St. Gallen, [samuel.haefner@unisg.ch](mailto:samuel.haefner@unisg.ch) 

------
------
## Overview

The julia files contain both the functions for estimation as well as the scripts to read and further process the estimates.

The following main files contain all relevant functions and global variables (for more detailed information see below):

- **Auxiliary.jl** --- reads in the required packages, the data set and defines the required auxiliary functions and  global variables.
- **Grouping.jl** --- determines the bidder groups and auction groups used for the estimation.
- **Estimation.jl** --- contains the main functions for the estimation of W and the bounds.
- **TestMon.jl** --- contains the main functions to test for the monotonicity of F in v.

The following files contain the generic scripts used to produce the scripts for the actual estimation (for more information, see the files themselves). The estimation was conducted at [sciCORE](http://scicore.unibas.ch/) scientific computing core facility at the University of Basel. The Julia version used was 1.2. Computation time depended on the script, ranging from 2-3 hours (WandTGeneric.jl and TestMonGeneric.jl) up to 48 − 72 hours (EstimateBoundsGeneric.jl).

- **EstimateWandTGeneric.jl** --- generic script to compute and save estimates of W and Theta.
- **EstimateWandTRobustGeneric.jl** --- generic script to compute and save robustness checks for W and Theta, using a log-normal distribution rather than a gamma distribution when estimating W.
- **EstimateBoundsGeneric.jl** --- generic script to compute standard and tighter bounds, using W from EstimateWandTGeneric.jl. 
- **TestMonGeneric.jl** --- generic script to test for monotonicity of F, determining W and the bounds for one bootstrap round, computing violations of monotonicity.

The following files contain the scripts to produce the .jl and .sh files that are needed reto run the estimates on a SLURM workload manager. For each *Generic.jl-File above, they split up the jobs into a manageable number of bootstrap rounds. 
- **EstimateWandTSLURM.jl**
- **EstimateWandTRobustSLURM.jl**
- **EstimateBoundsSLURM.jl**
- **TestMonSLURM.jl**

The following files contain the scripts used to read in and process the estimates produced with above *Slurm.jl scripts:

- **ReadEstimates.jl** --- contains the functions to read in the estimates and to produce the plot of the bounds.
- **ReadT.jl** --- script to read in estimates of T and produce the plots. 
- **ReadTRobust.jl** --- same as above, but using the alternative estimates (robustness check).
- **ReadTestMon.jl** --- script to read in the monotonicity violations obtained and saved in TestMonGeneric.jl.

The last file contains scripts used for the plots:

- **ResamplingPlots.jl** --- scripts to generate the group plots and the resampling plots. 

------
------

## Estimation.jl

### Estimation of W:

- **Wgamma(prices, auctionset, bidderassignment, n, m, P)**  
returns estimates of W(p,q) for the prices in the vector *prices*, estimating a gamma distribution
  - **prices** -- vector of submitted prices to be used for the estimation of W.
  - **auctionset** -- vector with indeces of auctions to be used for the estimation of W.
  - **bidderassignment** -- bidder assignment vector
  - **n** -- vector, each entry corresponding to the (average) number of active bidders from a given group
  - **m** -- number of bootstrap rounds to be estimated
  - **P** -- number of rounds of resampling used for estimation

- **Wlnorm(prices, auctionset, bidderassignment, n, m, P)**  
same as above, estimating a log normal distribution

### Estimation of Simple Bounds:

- **SimpleBound(bid, WPar, rho, Q)**  
returns a dataframe with simple bounds, v=[vub,vlb,q], based on *bid* and *WPar*
  - **bid** -- bid function
  - **WPar** -- array, containing the estimated parameters of the distribution of (p_i^j) for steps j=1,...,k
  - **rho** -- positive real number, corresponding to the risk preference $\rho$
  - **Q** -- positive real number, corresponding to the quota $Q$

- **EstimateSimpleBounds(auction, W, bidderassignment,  prices, rho, m)**  
returns estimates of simple bounds for all bidders in *auction* using bids from *auctionset* assuming a risk preference *rho*
  - **auction** -- auction index
  - **W** -- estimate of W, as returned from Wgamma() or Wlnorm()
  - **bidderassignment** -- bidder assignment vector
  - **prices** -- vector of submitted prices used for the estimation of W
  - **rho** -- positive real number, corresponding to the risk preference $\rho$
  - **m** -- number of bootstrap rounds to be estimated


### Estimation of Theta:

- **EstTheta(auction, W, bidderassignment, prices, bounds, rho, m)**  
returns share of inequality violations for a given auction

  - **auction** -- auction index
  - **W** -- estimate of W, as returned from Wgamma() or Wlnorm()
  - **bidderassignment** -- bidder assignment vector
  - **prices** -- vector of prices used for the estimation of W
  - **bounds** -- estimated simple bounds from EstimateSimpleBounds()
  - **rho** -- positive real number, the risk preference 
  - **m** -- number of bootstrap rounds
  

### Estimation of Tighter Bounds:

- **TighterBounds(bid, initvl, initvu, W, g, prices, rho, Q, bootstraprun, maxiter, tolerance)**  
returns tighter bounds from initial conditions *initvl* and *initivu*, runing Algorithm 2

  - **bid** -- bid function
  - **initvl** -- marginal profit function, the initial condition for the lower bound
  - **initvu** -- marginal profit function, the initial condition for the upper bound
  - **W** -- estimate of W, as returned from Wgamma() or Wlnorm()
  - **g** -- natural number, indicating the group number of the auction
  - **prices** -- vector of prices used for the estimation of W
  - **rho** -- positive real number, the risk preference 
  - **bootstraprun** -- index of boostratp run we look at
  - **maxiter** -- maximum number of fixed point iterations 
  - **tolerance** -- tolerance level (used in iteration)
   
- **EstTighterBounds(auction, W, bidderassignment, prices, bounds, rho, m, maxiter, tolerance)**  
returns tighter bounds for all bidders in an *auction* and all *m* bootstrap rounds

  - **auction** -- auction index
  - **W** -- estimate of W, as returned from Wgamma() or Wlnorm()
  - **bidderassignment** -- bidder assignment vector
  - **prices** -- vector of prices used for the estimation of W
  - **bounds** -- estimated simple bounds from EstimateSimpleBounds()
  - **rho** -- positive real number, the risk preference 
  - **m** -- number of boostratp runs
  - **maxiter** -- maximum number of fixed point iterations 
  - **tolerance** -- tolerance level (used in iteration)

------
------

## TestMon.jl

### The main function:

- **TestIncreasingDiff(auction, bidderset, BoundsAuction, W, R, rhoindexset, rhovec)**  
returns a list of length of the *rhoindexset*, each entry is again a list of length of the *bidderset*, and each entry of that list consists of two numbers: (1) the total number of bid steps tested for that bidder (2) the total number of bidsteps for which a violation of monotonicity was found.

  - **auction** -- auction index
  - **bidderset** -- bidder index set (set of bidders to be looked at)
  - **BoundsAuction** -- first element in list returned by LoadWandBounds()
  - **W** -- second element in list returned by LoadWandBounds()
  - **R** -- number of tests per bid step
  - **rhoindexset** -- vector of indexes, referring to rhovec, for which the test is run
  - **rhovec** -- vector of rho values
 
### The auxiliary functions:
- **GetQVals(g, bidstep)**  
returns the xvalues of the initial step function *g* and a random number of random values between the x-values of the *bidstep*-th step and the *bidstep*-1-th step

  - **g** -- decreasing step function function with a number of steps of at least *bidstep*
  - **bidstep** -- natural number, indicating the number of the step in the bid function
  
- **GetHigherF(g, h, qval)**  
returns a random step function between *g* and *h* with x values given in *qval*

  - **g**, **h** -- decreasing step function, satisfying *g* > *h*
  - **qvals** -- x-values returned by GetQVals()

- **CheckFOCMonotone(bidstep,bid,vub,vlb,W,g,prices,rho,Q,bootstraprun,R)**  
returns *R* differences in F^j, evaluated at the two different, ordered profit functions

  - **bidstep** -- natural number, the bid step to be looked at
  - **bid** -- bid function
  - **vub** -- marginal profit function, (simple) upper bound
  - **vlb** -- marginal profit function, (simple) lower bound
  - **W** -- estimate of W as obtained by Wgamma() or Wlnorm()
  - **g** -- natural number, indicating the group number of the auction
  - **prices** -- vector of prices used for the estimation of W
  - **rho** -- real number, risk preference $\rho$
  - **Q** -- real number, the quota $Q$
  - **bootstraprun** -- the number of the bootrap round
  - **R** -- the number of tests per bid step

- **LoadWandBounds(auction)**  
returns list containing W and simple bounds for *auction*
  - **auction** -- auction index
  
- **TestIncreasingDiffBidder(auction, bidder, Bounds, W, R, rho)**  
returns two numbers: [number of tests performed (corresponds to the number of bid steps), number of violations]

  - **auction** -- auction index
  - **bidder** -- bidder index
  - **Bounds** -- first element in list returned by LoadWandBounds()
  - **W** -- second element in list returned by LoadWandBounds()
  - **R** -- number of tests per bid step
  - **rho** -- real number, risk preference $\rho$
 
------
------

## Auxiliary.jl

### Global variables:

- **bids**: a data frame containing the bidding data.
- **K**: max number of price-quantity pairs.
- **un**: unit of account (100 corresponds to CHF/kg).
- **vupperbar**: upper bound on type space (in CHF).
- **auctionindeces**: vector containing the auctionindeces corresponding to those in the data.
- **bidderindeces**: vector containing the bidderindeces corresponding to those in the data.
- **quotas**: vector with the quotas.
- **clearingprices**: vector with the market clearing prices.
- **activebidders**: vector with the number of active bidders.
- **activebidderindeces**: list of vectors with the indeces of the active bidders.

### Auxiliary functions:

- **qpBid(bidder, auction)**  
returns the bid function of a bidder in an auction, which consists in dataframe with columns *qb* (quantity points), *pb* (price points), *cumqb* (cumulated quantity points)

  - **bidder** -- bidderindex  
  - **auction** -- auctionindex  

- **AvgNoBidders(bidderassignment)**  
returns the average number of bidders in each bidder cluster
  
  - **bidderassignment** -- vector, each element corresponding to a bidder, denoting the cluster number of the respective bidder  
  
- **BidToCover(auction)**  
returns the bid-to-cover ratio in an auction

  - **auction** -- auctionindex

- **Revenue(auction)**  
returns the revenue in an auction

  - **auction** -- auctionindex
  
- **qRec(bidder, auction)**  
returns the received quantity of a bidder in a given auction

  - **bidder** -- bidderindex  
  - **auction** -- auctionindex

- **qShareRec(bidder, auction)**  
returns the received quantity share of a bidder in a given auction

  - **bidder** -- bidderindex  
  - **auction** -- auctionindex
  
- **ActiveAuctions(bidder)**  
returns the indeces in which a bidder is active

  - **bidder** -- bidderindex 

- **AvgBid(bidder)**  
returns the average bid of a bidder

  - **bidder** -- bidderindex
- **ShareSuccBidders(auction)**  
returns the  share of succesfull bidders in an auction

  - **auction** -- auctionindex

- **SuccessRate(bidder)**  
returns the success rate of a bidder

  - **bidder** -- bidderindex

- **StepBid(q, bid)**
returns the value of $\beta_b(q)$ (cf. the manuscript for a definition)

  - **q** -- positive real number
  - **bid** -- bid function as returned from qpBid(bidder, auction)

- **NoStepBid(q, bid)**  
returns the step number at *q* of the step function $\beta_b$

  - **q** positive real number
  - **bid** bid function as returned from qpBid(bidder, auction)

- **StepV(q, v)**  
returns the value of v(q)

  - **q** -- positive real number
  - **v** -- marginal profit function, which is a data frame with columns *qval* (quantities, need to be increasing) and *vval* (the values of v)

- **IntBid(a, b, bid)**  
returns the value of $\int_a^b\beta_b(q)dq$
 
  - **a**, **b** -- positive real numbers
  - **bid** -- bid function as returned from qpBid(bidder, auction)
  
- **IntV(a, b, v)**  
returns the value of $\int_a^b v(q)dq$

  - **a**, **b** -- positive real numbers
  - **v** -- marginal profit function

- **PiOverline(bidstep, bid, v, WPar, rho, Q, n)**  
returns the value of $\overline{\Pi}_i^j(b,v)$

  - **bidstep** -- positive natural number, corresponding to the number of the step under consideration, j
  - **bid** -- bid function
  - **v** -- marginal profit function
  - **WPar** -- array, containing the estimated parameters of the distribution of D(p_i^j) for steps j=1,...,k
  - **rho** -- positive real number, corresponding to the risk preference $\rho$
  - **Q** -- positive real number, corresponding to the quota $Q$
  - **n** -- positive natural number, indicating the number of support points used for integration

- **w(q, WPar, dp)**  
returns w(p,q)

  - **q** -- positive real number
  - **WPar** -- array, containing the estimated parameters of the distribution of D(p) both for p_i^j and p_i^j + dp, j=1,...,k, 
  - **dp** -- positive real number

- **PriceBids(auctionset)**  
returns all submitted prices in *auctionset* in ascending order

  - **auctionset** -- vector, containing auction indeces

- **FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)**  
returns the value of the FOC, called $F^j$ in the text

  - **bidstep** -- positive natural number, corresponding to the number of the step under consideration, j
  - **bid** -- bid function
  - **v** -- marginal profit function
  - **W** -- estimates of W as returned from Wgamma() or Wlnorm(), cf. below
  - **group** -- natural number, indicating the group number of the auction
  - **prices** -- vector of submitted prices used for the estimation of W
  - **rho** -- positive real number, corresponding to the risk preference $\rho$
  - **Q** -- positive real number, corresponding to the quota $Q$
  - **bootstraprun** -- natural number, indicating the boostrap run number 
  - **n** -- positive natural number, indicating the number of support points used for integration

- **VarPhiU(q, v, vl)**  
returns $\varphi_u(q,v,v_l)$ (in the text)

  - **q** -- positive real number
  - **v** -- marginal profit function
  - **vl** -- positive real number, corresponding to $v_l$

- **VarPhiL(q, v, vu)**  
returns $\varphi_l(q,v,v_u)$ (in the text)

  - **q** -- positive real number
  - **v** -- marginal profit function
  - **vu** -- positive real number, corresponding to $v_u$

------
------

## Grouping.jl

- **group**: list containing, for each auction group, a vector of the respective auction indeces 
- **bidderassignment**: a vector determining for every bidder the bidder group assignment {1,2,3}

------