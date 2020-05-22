This file contains information about the Julia code used in:

# **[Risk Aversion in Share Auctions: Estimating Rents from TRQs in Switzerland](https://dx.doi.org/10.2139/ssrn.3397027)**

Author: [Samuel Häfner](https://samuelhaefner.github.io), University of St. Gallen, [samuel.haefner@unisg.ch](mailto:samuel.haefner@unisg.ch) 


# Overview

The replication files contain both the functions for estimation as well as the scripts to read and further process the estimates. 

## Main Files
The following main files contain all relevant functions and global variables (for more detailed information on the functions and variables in these files see the sections below):

- ```Auxiliary.jl``` --- reads in the required packages, the data set and defines the required auxiliary functions and  global variables.  
- ```Grouping.jl``` --- determines the bidder groups and auction groups used for the estimation.  
- ```Estimation.jl``` --- contains the main functions for the estimation of W, $\Theta$, and the bounds.  
- ```TestMon.jl``` --- contains the main functions to test for the monotonicity of $F^j$ in v.

## Generic Scripts
The following files contain the generic scripts used to produce the scripts for the actual estimation (for more information, see the files themselves). The estimation was conducted at [sciCORE](http://scicore.unibas.ch/) scientific computing core facility at the University of Basel. The Julia version used was 1.2. Computation time depended on the script, ranging from 2-3 hours (WandTGeneric.jl and TestMonGeneric.jl) up to 48 − 72 hours (EstimateBoundsGeneric.jl).

- ```EstimateWandTGeneric.jl``` --- generic script to compute and save estimates of W and Theta.  
- ```EstimateWandTRobustGeneric.jl``` --- generic script to compute and save robustness checks for W and Theta, using a log-normal distribution rather than a gamma distribution when estimating W.  
- ```EstimateBoundsGeneric.jl``` --- generic script to compute standard and tighter bounds, using W from EstimateWandTGeneric.jl.   
- ```TestMonGeneric.jl``` --- generic script to test for monotonicity of F, determining W and the bounds for one bootstrap round, computing violations of monotonicity.

## Scripts for SLURM
The following files contain the scripts to produce the .jl and .sh files that are needed reto run the estimates on a SLURM workload manager. For each *Generic.jl-File above, they split up the jobs into a manageable number of bootstrap rounds.  

- ```EstimateWandTSLURM.jl```  
- ```EstimateWandTRobustSLURM.jl```  
- ```EstimateBoundsSLURM.jl```  
- ```TestMonSLURM.jl```  

## Scripts to Process the Estimates
The following files contain the scripts used to read in and process the estimates produced with above *Slurm.jl scripts:

- ```ReadEstimates.jl``` --- contains the functions to read in the estimates and to produce the plot of the bounds.
- ```ReadT.jl``` --- script to read in estimates of T and produce the plots. 
- ```ReadTRobust.jl``` --- same as above, but using the alternative estimates (robustness check).
- ```ReadTestMon.jl``` --- script to read in the monotonicity violations obtained and saved in TestMonGeneric.jl.

## Scripts for Plots in the Paper
The last file contains scripts used for the plots:

- ```ResamplingPlots.jl``` --- scripts to generate the group plots and the resampling plots. 

# The Data  
The data are contained in the file ```setofbids.csv```. Each of the 12401 rows corresponds to a submitted price-quantity pair. The columns are the following:

| Variable | Description |
|--- | --- |
|```auction```| auction id |    
|```quotatot``` | quota (in kg) | 
|```bid_id``` | bid id |  
|```bidder``` | bidder id | 
|```qb```| quantity point (in kg) |  
|```pb```| price point (in Swiss cents) | 
|```qr```| resulting quantity |  
|```pr```| resulting payment | 
|```qperc```| percentage of total quota |  
|


# An Example

The following code goes through the basic computations. Doing one bootstrap round, the code first estimates W(p,q), then computes $\Theta(\rho)$, and last determines the bounds for $\rho=0$ and $\rho=\rho^*$. 

```julia
###############
## load the data set. 
## define global variables and relevant functions
###############
include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")
include("TestMon.jl")

###############
## the following lines compute the estimates of W(p,q).  
###############

# obtain the average number of active bidders for each bidder group  
n = AvgNoBidders(bidderassignment)

# do m=1 bootstrap round
m = 1

# draw P=200 opponent demand functions when resampling (cf. Algorithm 1)
P = 200

# W is computed for each auction group separately
W = []
for i in [1:1:length(group);]  
    auctionset = group[i]
    prices = PriceBids(auctionset)
    push!(W, Wgamma(prices, auctionset, bidderassignment, n, m, P))
end

###############
## the following computes all the values required 
## to compute $\Theta(\rho)$ for the values in rhovec. 
###############

# define the vector of rho values
rhovec = sort!(append!([exp(x) for x in [-10:0.5:0;]], 0))

# the values are computed for each auction 
# in each auction group separately
Theta = []
for g in [1:1:length(group);]
    auctionset = group[g]
    prices = PriceBids(auctionset)
    ThetaGroup = []
    for auction in auctionset
        ThetaAuction = []
        for i in [1:1:length(rhovec);]
            bounds = EstimateSimpleBounds(
                auction,
                W[g],
                bidderassignment,
                prices,
                rhovec[i],
                m,
            )
            push!(
                ThetaAuction,
                EstTheta(
                    auction,
                    W[g],
                    bidderassignment,
                    prices,
                    bounds,
                    rhovec[i],
                    m,
                ),
            )
        end
        push!(ThetaGroup, ThetaAuction)
    end
    push!(Theta, ThetaGroup)
end

###############
## the following code computes both the simple and the tighter bounds 
## for rho=0 and rho=rho^* (only the first for rho=0, and 
## both for rho=rho*) for all active bidders in a given auction
###############

# define the two values for rho
rhovec = [0, exp(-5.5)]

# take the second auction in the first auction group, g=1
g = 1
auction = group[g][2]
BoundsAuction = []
for i in [1:1:length(rhovec);]
  simplebounds = EstimateSimpleBounds(
    auction,
    W[g],
    bidderassignment,
    prices,
    rhovec[i],
    m,
  )
  if i == 1
    push!(BoundsAuction, simplebounds)
  else
    push!(
      BoundsAuction,
      EstTighterBounds(
        auction,
        W[g],
        bidderassignment,
        prices,
        simplebounds,
        rhovec[i],
        m,
        10,
        0.0001,
      ),
    )
  end
end
```

# Estimation.jl


```
Wgamma(prices, auctionset, bidderassignment, n, m, P) 
```
#### Description
Estimates W(p,q) using a gamma distribution. 
#### Arguments
```prices``` -- vector of submitted prices to be used for the estimation of W.  
```auctionset``` -- vector with indeces of auctions to be used for the estimation of W.  
```bidderassignment``` -- bidder assignment vector  
```n``` -- vector, each entry corresponding to the (average) number of active bidders from a given group  
```m``` -- number of bootstrap rounds to be estimated  
```P``` -- number of rounds of resampling used for estimation
#### Return value
A list of parameter estimates for the distribution of W(p,q). For each group in ```bidderassignment```, each bootstrap round, and each price in the vector ```prices```.

-----
```
Wlnorm(prices, auctionset, bidderassignment, n, m, P)
```
#### Description
Same as ```Wgamma()```, yet using a log normal distribution.

------


```
SimpleBound(bid, WPar, rho, Q)
```
#### Description
Computes upper and lower bounds on the rationalizable profit functions as described in Proposition 3 and equations (9)-(10).
#### Arguments
```bid``` -- bid function  
```WPar``` -- array, containing the estimated parameters of W(p,q) at (p_i^j) for steps j=1,...,k  
```rho``` -- positive real number, corresponding to the risk preference $\rho$  
```Q``` -- positive real number, corresponding to the quota $Q$
#### Return value
A list of dataframes, [vlb,vub], where vub is the upper bound and vlb is the lower bound.

-----
```
EstimateSimpleBounds(auction, W, bidderassignment,  prices, rho, m)
```
#### Description
Computes ```SimpleBounds()``` for all bidders in ```auction```.
#### Arguments
```auction``` -- auction index  
```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
```bidderassignment``` -- bidder assignment vector  
```prices``` -- vector of submitted prices used for the estimation of W  
```rho``` -- positive real number, corresponding to the risk preference $\rho$  
```m``` -- number of bootstrap rounds to be estimated
#### Return value
A list of objects returned by ```SimpleBounds()```, one entry per bidder.

-----

```
EstTheta(auction, W, bidderassignment, prices, bounds, rho, m)
``` 
#### Description
Determines violations of the inequalities in Proposition 4 for a given auction.
#### Arguments
  ```auction``` -- auction index  
  ```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```   
  ```bidderassignment``` -- bidder assignment vector  
  ```prices``` -- vector of prices used for the estimation of W  
  ```bounds``` -- estimated simple bounds from ```EstimateSimpleBounds()```  
  ```rho``` -- positive real number, the risk preference   
  ```m``` -- number of bootstrap rounds  
#### Return value
A list of two matrices, with columns corresponding to bootrap round and rows corresponding to bidder groups. The first matrix counts the number of violations, the second matrix counts the total number of test inequalities.  

-----
```
TighterBounds(bid, initvl, initvu, W, g, prices, rho, Q, bootstraprun, maxiter, tolerance)
```
#### Description
Computes tighter upper and lower bounds from initial conditions ```initvl``` and ```initivu``` by runing Algorithm 2.
#### Arguments
 ```bid``` -- bid function  
 ```initvl``` -- marginal profit function, the initial condition for the lower bound  
 ```initvu``` -- marginal profit function, the initial condition for the upper bound  
 ```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
 ```g``` -- natural number, indicating the group number of the auction  
 ```prices``` -- vector of prices used for the estimation of W  
```rho``` -- positive real number, the risk preference   
 ```bootstraprun``` -- index of boostratp run we look at  
 ```maxiter``` -- maximum number of fixed point iterations   
 ```tolerance``` -- tolerance level (used in iteration)  
 #### Return value
A list of dataframes, [vlb,vub], where vub is the upper bound and vlb is the lower bound.

-----

```
EstTighterBounds(auction, W, bidderassignment, prices, bounds, rho, m, maxiter, tolerance)
```
#### Description
Runs ```TighterBounds()``` for all bidders in an ```auction``` and for ```m``` bootstrap rounds.
#### Arguments

```auction``` -- auction index  
```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
 ```bidderassignment``` -- bidder assignment vector  
 ```prices``` -- vector of prices used for the estimation of ```W```  
 ```bounds``` -- estimated simple bounds from ```EstimateSimpleBounds()```  
```rho``` -- positive real number, the risk preference  
```m``` -- number of boostratp runs  
```maxiter``` -- maximum number of fixed point iterations  
```tolerance``` -- tolerance level (used in iteration)  
#### Return value
List containing for each biddera and each bootstrapround a two dimensional list [```bounds```,tighterbounds], where tighterbounds is the object returned by ```TighterBounds()```.


# TestMon.jl

```
TestIncreasingDiff(auction, bidderset, BoundsAuction, W, R, rhoindexset, rhovec)
```  
#### Description
Determines the number of price-quantity pairs for which the monotonicity assumption on $F^j$ is violated for a given auction.
#### Arguments
  ```auction``` -- auction index  
  ```bidderset``` -- bidder index set (set of bidders to be looked at)  
  ```BoundsAuction``` -- first element in list returned by ```LoadWandBounds()```  
  ```W``` -- second element in list returned by ```LoadWandBounds()```  
  ```R``` -- number of tests per bid step  
  ```rhoindexset``` -- vector of indexes, referring to rhovec, for which the test is run  
  ```rhovec``` -- vector of rho values  
 #### Return value
A list of length of the ```rhoindexset```, each entry is again a list of length of the ```bidderset```, and each entry of that list consists of two numbers: (1) the total number of price-quantity pairs tested for that bidder (2) the total number of price-quantity pairs for which a violation of monotonicity was found.

----
```
GetQVals(g, bidstep)
```  
#### Description
Draws a random number (uniformly between 1 and 20) of equidistant x-values between the x-values of the ```bidstep```-th price-quantity pair in ```g``` and the (```bidstep```-1)-th price-quantity pair
#### Arguments
  ```g``` -- decreasing step function function with a number of steps of at least *bidstep*  
  ```bidstep``` -- natural number, indicating the number of the step in the bid function  
#### Return value
List of x-values.

----
```
GetHigherF(g, h, qval)
```
#### Description
Compute a random step function between ```g``` and ```h``` with x values given in ```qval```
#### Arguments
  ```g```, ```h``` -- decreasing step function, satisfying ```g``` > ```h```    
  ```qvals``` -- x-values returned by ```GetQVals()```
#### Return value
Data frame containing the step function, [vval,qval].

-----

```
CheckFOCMonotone(bidstep, bid, vub, vlb, W, g, prices, rho, Q, bootstraprun, R)
```
#### Description
Evaluates $F^j$ at the two different, ordered profit functions ```R``` times and reports the differences.
#### Arguments
  ```bidstep``` -- natural number, the bid step to be looked at  
  ```bid``` -- bid function  
  ```vub``` -- marginal profit function, (simple) upper bound  
  ```vlb``` -- marginal profit function, (simple) lower bound  
  ```W``` -- estimate of W as obtained by ```Wgamma()``` or ```Wlnorm()```
  ```g``` -- natural number, indicating the group number of the auction  
  ```prices``` -- vector of prices used for the estimation of W  
  ```rho``` -- real number, risk preference $\rho$  
  ```Q``` -- real number, the quota $Q$  
  ```bootstraprun``` -- the number of the bootrap round  
  ```R``` -- the number of tests per bid step  
#### Return value
List of length ```R```, containing the differences.

----
```
LoadWandBounds(auction)
```
#### Description
Loads estimates of W and of simple bounds for ```auction``` from corresponding .dat files.
#### Arguments
  ```auction``` -- auction index  
#### Return value
List containing W and simple bounds for ```auction```

----

```
TestIncreasingDiffBidder(auction, bidder, Bounds, W, R, rho)
```
#### Description
Determines for a bidder in a given auction whether monotonicity of $F^j$ is violated for the submitted price-quantity pairs. 
#### Arguments
  ```auction``` -- auction index  
  ```bidder``` -- bidder index  
  ```Bounds``` -- first element in list returned by ```LoadWandBounds(```)  
  ```W``` -- second element in list returned by ```LoadWandBounds()```  
  ```R``` -- number of tests per bid step  
  ```rho``` -- real number, risk preference $\rho$  
#### Return value
Two-dimensional list, [number of tests performed (corresponds to the number of bid steps), number of violations]



# Auxiliary.jl

#### Global variables:
The script defines the following global variables that are used throughout.

```bids```: a data frame containing the bidding data  
```K```: max number of price-quantity pairs  
```un```: unit of account (100 corresponds to CHF/kg)  
```vupperbar```: upper bound on type space (in CHF)  
```auctionindeces```: vector containing the auctionindeces corresponding to those in the data  
```bidderindeces```: vector containing the bidderindeces corresponding to those in the data  
```quotas```: vector with the quotas  
```clearingprices```: vector with the market clearing prices  
```activebidders```: vector with the number of active bidders  
```activebidderindeces```: list of vectors with the indeces of the active bidders  

-----

```
qpBid(bidder, auction)
```
#### Description
```qpBid(bidder, auction)``` returns the bid function of a bidder in an auction, which is a dataframe [qb (quantity points), pb (price points), cumqb (cumulated quantity points)]
#### Arguments
  ```bidder``` -- bidderindex  
  ```auction``` -- auctionindex  
#### Return value
```qpBid(bidder, auction)``` returns the bid function of a bidder in an auction, which is a dataframe [qb (quantity points), pb (price points), cumqb (cumulated quantity points)]

-----
```
AvgNoBidders(bidderassignment)
```
#### Description
```AvgNoBidders()``` returns the average number of bidders in each bidder cluster
#### Arguments
  ```bidderassignment``` -- vector, each element corresponding to a bidder, denoting the cluster number of the respective bidder  
#### Return value
```AvgNoBidders()``` returns the average number of bidders in each bidder cluster

-----
```
BidToCover(auction)
```
#### Description
```BidToCover()``` returns the bid-to-cover ratio in an auction
#### Arguments
  ```auction``` -- auctionindex
#### Return value
```BidToCover()``` returns the bid-to-cover ratio in an auction

-----
```
Revenue(auction)
```
#### Description
```Revenue()``` returns the revenue in an auction
#### Arguments
  ```auction``` -- auctionindex
#### Return value
```Revenue()``` returns the revenue in an auction

-----
```
qRec(bidder, auction)
```
#### Description
```qRec()``` returns the received quantity of a bidder in a given auction
#### Arguments
  ```bidder``` -- bidderindex   
  ```auction``` -- auctionindex
#### Return value
```qRec()``` returns the received quantity of a bidder in a given auction

----
```
qShareRec(bidder, auction)
```
#### Description
```qShareRec()```  returns the received quantity share of a bidder in a given auction
#### Arguments
  ```bidder``` -- bidderindex  
  ```auction``` -- auctionindex
#### Return value
```qShareRec()```  returns the received quantity share of a bidder in a given auction

----
```
ActiveAuctions(bidder)
```
#### Description
```ActiveAuctions()``` returns the indeces in which a bidder is active
#### Arguments
  ```bidder``` -- bidderindex 
#### Return value
```ActiveAuctions()``` returns the indeces in which a bidder is active

-----
```
AvgBid(bidder)
```
#### Description
```AvgBid()``` returns the average bid of a bidder
#### Arguments
  ```bidder``` -- bidderindex
#### Return value
```AvgBid()``` returns the average bid of a bidder

-----
```
ShareSuccBidders(auction)
```
#### Description
```ShareSuccBidders()```  returns the  share of succesfull bidders in an auction
#### Arguments
  ```auction``` -- auctionindex
#### Return value
```ShareSuccBidders()```  returns the  share of succesfull bidders in an auction

----
```
SuccessRate(bidder)
```
#### Description
```SuccessRate()``` returns the success rate of a bidder
#### Arguments
  ```bidder``` -- bidderindex
#### Return value
```SuccessRate()``` returns the success rate of a bidder

----
```
StepBid(q, bid)
```
#### Description
```StepBid()``` returns the value of $\beta_b(q)$ (cf. the manuscript for a definition)
#### Arguments
  ```q``` -- positive real number
  ```bid``` -- bid function as returned from qpBid(bidder, auction)
#### Return value
```StepBid()``` returns the value of $\beta_b(q)$ (cf. the manuscript for a definition)

----
```
NoStepBid(q, bid)
```
#### Description
```NoStepBid()```  returns the step number at ```q``` of the step function $\beta_b$
#### Arguments
  ```q``` positive real number
  ```bid``` bid function as returned from qpBid(bidder, auction)
#### Return value
```NoStepBid()```  returns the step number at ```q``` of the step function $\beta_b$

-----
```
StepV(q, v)
```
#### Description
```StepV()``` returns the value of v(q)
#### Arguments
  ```q``` -- positive real number
  ```v``` -- marginal profit function, which is a data frame with columns ```qval``` (quantities, need to be increasing) and ```vval``` (the values of v)
#### Return value
```StepV()``` returns the value of v(q)

-----
 ```
 IntBid(a, b, bid)
 ```
 #### Description
 ```IntBid()``` returns the value of $\int_a^b\beta_b(q)dq$
 #### Arguments
  ```a```, ```b``` -- positive real numbers
  ```bid``` -- bid function as returned from qpBid(bidder, auction)
#### Return value
```IntBid()``` returns the value of $\int_a^b\beta_b(q)dq$

----
```
IntV(a, b, v)
```
#### Description
```IntV()``` returns the value of $\int_a^b v(q)dq$
#### Arguments
  ```a```, ```b``` -- positive real numbers
  ```v``` -- marginal profit function
#### Return value
```IntV()``` returns the value of $\int_a^b v(q)dq$

----
```
PiOverline(bidstep, bid, v, WPar, rho, Q, n)
```
#### Description
```PiOverline()``` returns the value of $\overline{\Pi}_i^j(b,v)$
#### Arguments
  ```bidstep``` -- positive natural number, corresponding to the number of the step under consideration, j  
  ```bid``` -- bid function  
  ```v``` -- marginal profit function  
  ```WPar``` -- array, containing the estimated parameters of the distribution of D(p_i^j) for steps j=1,...,k  
  ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
  ```Q``` -- positive real number, corresponding to the quota $Q$  
  ```n``` -- positive natural number, indicating the number of support points used for integration  
#### Return value
```PiOverline()``` returns the value of $\overline{\Pi}_i^j(b,v)$

-----
```
w(q, WPar, dp)
```
#### Description
```w()``` returns w(p,q)
#### Arguments
  ```q``` -- positive real number  
  ```WPar``` -- array, containing the estimated parameters of the distribution of D(p) both for p_i^j and p_i^j + dp, j=1,...,k   
  ```dp``` -- positive real number  
#### Return value
```w()``` returns w(p,q)

------
```
PriceBids()
```
#### Description
```PriceBids(auctionset)``` returns all submitted prices in ```auctionset``` in ascending order
#### Arguments
  ```auctionset``` -- vector, containing auction indeces
#### Return value
```PriceBids(auctionset)``` returns all submitted prices in ```auctionset``` in ascending order

-----
```
FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
```
#### Description
```FOC()``` returns the value of the FOC, called $F^j$ in the text
#### Arguments
  ```bidstep``` -- positive natural number, corresponding to the number of the step under consideration, j  
  ```bid``` -- bid function  
  ```v``` -- marginal profit function  
  ```W``` -- estimates of W as returned from ```Wgamma()``` or ```Wlnorm()```  
  ```group``` -- natural number, indicating the group number of the auction  
  ```prices``` -- vector of submitted prices used for the estimation of W  
  ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
  ```Q``` -- positive real number, corresponding to the quota $Q$  
  ```bootstraprun``` -- natural number, indicating the boostrap run number   
  ```n``` -- positive natural number, indicating the number of support points used for integration  
#### Return value
```FOC()``` returns the value of the FOC, called $F^j$ in the text

-----
```
VarPhiU(q, v, vl)
```
#### Description
```VarPhiU()``` returns $\varphi_u(q,v,v_l)$ (in the text)
#### Arguments
  ```q``` -- positive real number  
  ```v``` -- marginal profit function  
  ```vl``` -- positive real number, corresponding to $v_l$  
#### Return value
```VarPhiU()``` returns $\varphi_u(q,v,v_l)$ (in the text)

-----
```
VarPhiL(q, v, vu)
```
#### Description
```VarPhiL()``` returns $\varphi_l(q,v,v_u)$ (in the text)
#### Arguments
  ```q``` -- positive real number  
  ```v``` -- marginal profit function  
  ```vu``` -- positive real number, corresponding to $v_u$  
#### Return value
```VarPhiL()``` returns $\varphi_l(q,v,v_u)$ (in the text)

# Grouping.jl

This script defines the following global variables used for the estimation of W:

```group```: list containing, for each auction group, a vector of the respective auction indeces  
```bidderassignment```: a vector determining for every bidder the bidder group assignment {1,2,3}

