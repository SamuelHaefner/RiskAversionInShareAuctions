This file contains the required information for the replication of the results in:

# **[Risk Aversion in Share Auctions: Estimating Rents from TRQs in Switzerland](https://dx.doi.org/10.2139/ssrn.3397027)**

Author: [Samuel Häfner](https://samuelhaefner.github.io), Web3 Foundation, Zug, Switzerland, [samuel@web3.foundation](mailto:samuel@web3.foundation).

# Table of Contents

1. [Overview](#Overview)  
2. [Data](#Data)
3. [Estimates](#Estimates) 
4. [Example](#Example)  
5. [Scripts](#Scripts)  
  a. [Estimation.jl](##Estimation.jl)  
  b. [Testmon.jl](##Testmon.jl)  
  c. [Auxiliary.jl](##Auxiliary.jl) 


# Overview of the replication files

This [repository](https://github.com/SamuelHaefner/RiskAversionInShareAuctions) contains the following items:
1. The *main code* used for estimation. The code comes in four main files: 
   - ```Auxiliary.jl``` - Reads in the required packages, the data set and defines the required auxiliary functions and  global variables.  
   - ```Grouping.jl``` - Determines the bidder groups and auction groups used for the estimation.  
   - ```Estimation.jl``` - Contains the main functions for the estimation of W, $\Theta$, and the bounds.  
   - ```TestMon.jl``` - Contains the main functions to test for the monotonicity of $F^j$ in v.   
2. The *scripts to conduct the estimation*. In particular, the following files contain the scripts to produce the .jl and .sh files that are needed to run the estimation on a SLURM workload manager.  
   - ```EstimateWandTSLURM.jl``` - Produces the required .jl and .sh files to compute and save estimates of W and Theta.  
   - ```EstimateWandTRobustSLURM.jl``` - Produces the required files to compute and save the robustness checks for Theta, using a log-normal distribution rather than a gamma distribution when estimating W.
   - ```EstimateBoundsSLURM.jl``` - Produces the required files to compute standard and tighter bounds, using the estimates of W(p,q) obtained with EstimateWandTSLURM.jl.
   - ```TestMonSLURM.jl``` - Produces the required files to test for monotonicity of F.   
3. *Templates* used by the scripts above. Essentially, the scripts above split up the jobs into a manageable number of bootstrap rounds. To do so, they require the following generic .jl and .sh templates. For further information, see the files themselves.
   - ```EstimateWandTGeneric.jl```   
   - ```EstimateWandTRobustGeneric.jl```    
   - ```EstimateBoundsGeneric.jl```  
   - ```TestMonGeneric.jl```  
   - ```EstimateWandTGeneric.sh```  
   - ```EstimateWandTRobustGeneric.sh```
   - ```EstimateBoundsGeneric.sh```
   - ```TestMonGeneric.sh```
3. The scripts to *read and further process the estimates* produced with above *Slurm.jl scripts. Once the estimates are saved in the respective .dat files, these scripts can be run as they are.
   - ```ReadEstimates.jl``` - Contains the functions to read in the estimates and produce the respective tables.
   - ```ReadT.jl``` - Script to read in estimates of T and produce the plots and tables. 
   - ```ReadTRobust.jl``` - Same as above, but using the alternative estimates (robustness check).
   - ```ReadTestMon.jl``` - Contains the script to read in the monotonicity violations obtained and saved with TestMonGeneric.jl and produces the respective tables.
4. Script to produce additional plots.
   - ```Plots.jl``` - Script to generate the group plots, the resampling plots, and the plot showing the estimated bounds (which is done with function ```PlotTighterBounds()```; cf. the file for more information).  

*Comment 1:* The script ```TestMonWandBounds.jl``` produces the required files with the estimates of W(p,q) and the bounds for the monotonicity check. This script needs to be run before the TestMon*.jl scripts, using the bash script ```TestMonWandBounds.sh```.

*Comment 2:* Computation was conducted at [sciCORE](http://scicore.unibas.ch/) (scientific computing core) facility at the University of Basel, using a SLURM workload manager. The Julia version used was 1.5. Computation time depended on the script, ranging from 2-3 hours (WandTGeneric.jl and TestMonGeneric.jl) up to 48 − 72 hours (EstimateBoundsGeneric.jl).


# Data  
The data are contained in the file ```setofbids.csv```, which is also in the [repository](https://github.com/SamuelHaefner/RiskAversionInShareAuctions). Each row corresponds to a submitted price-quantity pair. The columns are the following:

| Variable | Description |
|--- | --- |
|```auction```| auction id |    
|```quotatot``` | quota (in kg) | 
|```bid_id``` | bid id |  
|```bidder``` | bidder id | 
|```qb```| quantity point (not cumulative, in kg) |  
|```pb```| price point (in Swiss cents) | 
|```qr```| resulting quantity (from that price-quantity pair)|  
|```pr```| resulting payment | 
|```qperc```| percentage of total quota |   

# Estimates
The estimates that I obtained and reported in the manuscript can be read in with the ```Read*.jl``` scripts as described in the [Overview section](#Overview). The estimates are saved in ```*.dat``` files. A zip archive of all files can be downloaded [here](https://drive.google.com/file/d/1bqyIMnCVvJJlmcStgyGKClkse68GCYYS/view?usp=sharing).


# Tables and Plots
The following table explains which scripts are used to produce the tables and plots in the [main paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3397027).

| Figure/Table | Script | Comment|
|---|---|---|
|Figure 1 | none | tikz-script in .tex file |
|Figure 2 | ```Plots.jl```| Lines 8-33 (left panel). Lines 74-82 (right panel). |
|Figure 3 | ```ReadT.jl``` | Lines 1-69; 115-211.|
|Figure 4 | ```Plots.jl```| Function ```PlotTighterBounds()```; see the file for further information. |
|Figure 5  | ```Plots.jl```| Lines 35-59 (left panel). Lines 64-72 (right panel). |
|Figure 6  |  ```Plots.jl``` | Lines 152-177 (left panel). Lines 37-46; 88-149 (right panel).|
|Table 1 | ```Auxilary.jl``` | Various functions in ```Auxiliary.jl``` are used. For a comprehensive overview see the corresponding [section](#Scripts) below. |
|Table 2 | ```ReadEstimates.jl``` | See file for more information. |
|Table 3 | ```ReadEstimates.jl``` | See file for more information. |

The following table explains which scripts are used to produce the tables and plots in the [supplementary appendix](https://samuelhaefner.github.io/SupplementaryAppendix.pdf) of the paper.

| Figure/Table | Script | Line Number |
|---|---|---|
| Figure 1 |```ReadT.jl``` | Lines 1-69; 213-467.|
| Figure 2 | ```ReatT.jl``` (right figure) and ```ReadTRobust.jl``` (left figure) | See files for more information. |
| Table 1 | ```ReadT.jl``` | Lines 1-69; 213-467. |
| Table 2 | ```ReadT.jl``` | See file for more information. |
| Table 3 | ```ReatT.jl``` (right table) and ```ReadTRobust.jl``` (left table) | See files for more information. |
| Tables in Section D | ```ReadEstimates.jl``` | |



# Example

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
## to compute $\Theta_g(\rho)$ for the values in rhovec. 
###############

# define the vector of rho values for with \Theta_g is computed
rhovec = sort!(append!([exp(x) for x in [-10:0.5:0;]], 0))

# the values are computed for each auction 
# in each auction group [g] separately
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
## for rho=(0,0,0) and rho=rho^* (only the first for rho=0, and 
## both for rho=rho*) for all active bidders in a given auction
###############

# define the two vectors of values for rho within the three groups
rhovec = [[0,0,0],[exp(-5),exp(-7),exp(-8)]]

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

# Scripts
This section provides the details about the functions and global variables defined in the four main files. 

## The main global variables
The **main global variables** used throughout the replication files are defined in the file ```Auxiliary.jl``` and are the followoing.

| Variable | Description |
|---|---|
|```bids```| a data frame containing the bidding data | 
|```K```| max number of price-quantity pairs|  
|```un```| unit of account (100 corresponds to CHF/kg)|  
|```vupperbar``` | upper bound on type space (in CHF)|  
|```auctionindeces```| vector containing the auctionindeces corresponding to those in the data| 
|```bidderindeces```| vector containing the bidderindeces corresponding to those in the data | 
|```quotas```| vector with the quotas|  
|```clearingprices```| vector with the market clearing prices | 
|```activebidders```| vector with the number of active bidders | 
|```activebidderindeces```| list of vectors with the indeces of the active bidders | 

In the file ```Grouping.jl``` the following, **additional global variables** used for the estimation of W are defined:

| Variable | Description |
|---|---|
|```group```| list containing, for each auction group, a vector of the respective auction indeces | 
|```bidderassignment```| a vector determining for every bidder the bidder group assignment {1,2,3}|

In the following subsections, I discuss the functions defined in the respective files ```Estimation.jl```, ```TestMon.jl```, and ```Auxiliary.jl```.

## Estimation.jl
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
A list of data frames, [vlb,vub], where vub is the upper bound and vlb is the lower bound.

-----
```
EstimateSimpleBounds(auction, W, bidderassignment,  prices, rhovec, m)
```
#### Description
Computes ```SimpleBounds()``` for all bidders in ```auction```.
#### Arguments
```auction``` -- auction index  
```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
```bidderassignment``` -- bidder assignment vector  
```prices``` -- vector of submitted prices used for the estimation of W  
```rhovec``` -- vector of positive real number, each number corresponding to the risk preference $\rho_g$ in bidder group $g$  
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
EstTighterBounds(auction, W, bidderassignment, prices, bounds, rhovec, m, maxiter, tolerance)
```
#### Description
Runs ```TighterBounds()``` for all bidders in an ```auction``` and for ```m``` bootstrap rounds.
#### Arguments

```auction``` -- auction index  
```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
 ```bidderassignment``` -- bidder assignment vector  
 ```prices``` -- vector of prices used for the estimation of ```W```  
 ```bounds``` -- estimated simple bounds from ```EstimateSimpleBounds()```  
```rhovec``` -- vector of positive real number, each number corresponding to the risk preference $\rho_g$ in bidder group $g$  
```m``` -- number of boostratp runs  
```maxiter``` -- maximum number of fixed point iterations  
```tolerance``` -- tolerance level (used in iteration)  
#### Return value
List containing for each biddera and each bootstrapround a two dimensional list [```bounds```,tighterbounds], where tighterbounds is the object returned by ```TighterBounds()```.


## TestMon.jl

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
Draws a random number (uniformly between 1 and 20) of equidistant x-values that lie between the x-values of the ```bidstep```-th price-quantity pair and the (```bidstep```-1)-th price-quantity pair in ```g```.
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

## Auxiliary.jl
```
qpBid(bidder, auction)
```
#### Description
Retrieve the bid function of a bidder in an auction.
#### Arguments
  ```bidder``` -- bidderindex  
  ```auction``` -- auctionindex  
#### Return value
A data frame: [qb (quantity points), pb (price points), cumqb (cumulated quantity points)].

-----
```
AvgNoBidders(bidderassignment)
```
#### Description
Determine the average number of bidders in each bidder cluster
#### Arguments
  ```bidderassignment``` -- vector, each element corresponding to a bidder, denoting the cluster number of the respective bidder  
#### Return value
A list of real numbers.

-----
```
BidToCover(auction)
```
#### Description
Determines the bid-to-cover ratio in an auction.
#### Arguments
  ```auction``` -- auctionindex
#### Return value
Real number.

-----
```
Revenue(auction)
```
#### Description
Determines the revenue from an auction.
#### Arguments
  ```auction``` -- auctionindex
#### Return value
Real number.

-----
```
qRec(bidder, auction)
```
#### Description
Determines the allocated quantity of a bidder in a given auction.
#### Arguments
  ```bidder``` -- bidderindex   
  ```auction``` -- auctionindex
#### Return value
Real number.

----
```
qShareRec(bidder, auction)
```
#### Description
Determines the share of quota that is allocated to a bidder in a given auction.
#### Arguments
  ```bidder``` -- bidderindex  
  ```auction``` -- auctionindex
#### Return value
Real number.

----
```
ActiveAuctions(bidder)
```
#### Description
Determines the indeces of the auctions in which a bidder is active
#### Arguments
  ```bidder``` -- bidderindex 
#### Return value
A list of integers.

-----
```
AvgBid(bidder)
```
#### Description
Determines the average bid of a bidder across the auctions in which that bidder was active.
#### Arguments
  ```bidder``` -- bidderindex
#### Return value
Real number.

-----
```
ShareSuccBidders(auction)
```
#### Description
Determines the  share of succesfull bidders (with non-zero allocated quantity) in an auction.
#### Arguments
  ```auction``` -- auctionindex
#### Return value
Real number.

----
```
SuccessRate(bidder)
```
#### Description
Returns the success rate of a bidder (succesfull participation/total participation).
#### Arguments
  ```bidder``` -- bidderindex
#### Return value
Real number.

----
```
StepBid(q, bid)
```
#### Description
Determines the value of $\beta_b(q)$ (cf. the manuscript for a definition).
#### Arguments
  ```q``` -- positive real number
  ```bid``` -- bid function as returned from qpBid(bidder, auction)
#### Return value
Real number.

----
```
NoStepBid(q, bid)
```
#### Description
Determines the step of the step function $\beta_b$ at ```q```. That is, returns the number of downward jumps that have occured strictly before ```q```, plus one.
#### Arguments
  ```q``` positive real number
  ```bid``` bid function as returned from qpBid(bidder, auction)
#### Return value
Integer.

-----
```
StepV(q, v)
```
#### Description
Returns the value of the profit function v(q).
#### Arguments
  ```q``` -- positive real number  
  ```v``` -- marginal profit function, which is a data frame with columns ```qval``` (quantities, need to be increasing) and ```vval``` (the values of v)
#### Return value
Real number.

-----
 ```
 IntBid(a, b, bid)
 ```
 #### Description
 Returns the value of $\int_a^b\beta_b(q)dq$
 #### Arguments
  ```a```, ```b``` -- positive real numbers
  ```bid``` -- bid function as returned from qpBid(bidder, auction)
#### Return value
Real number.

----
```
IntV(a, b, v)
```
#### Description
Returns the value of $\int_a^b v(q)dq$.
#### Arguments
  ```a```, ```b``` -- positive real numbers
  ```v``` -- marginal profit function
#### Return value
Real number.

----
```
PiOverline(bidstep, bid, v, WPar, rho, Q, n)
```
#### Description
Returns the value of $\overline{\Pi}_i^j(b,v)$ (cf. the manuscript for a definition).
#### Arguments
  ```bidstep``` -- positive natural number, corresponding to the number of the step under consideration, j  
  ```bid``` -- bid function  
  ```v``` -- marginal profit function  
  ```WPar``` -- array, containing the estimated parameters of the distribution of D(p_i^j) for steps j=1,...,k  
  ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
  ```Q``` -- positive real number, corresponding to the quota $Q$  
  ```n``` -- positive natural number, indicating the number of support points used for integration  
#### Return value
Real number.

-----
```
w(q, WPar, dp)
```
#### Description
Computes $w_i(p,q)$.
#### Arguments
  ```q``` -- positive real number  
  ```WPar``` -- array, containing two sets of estimated parameters of the distribution of D(p); one for p_i^j and one for p_i^j + dp.   
  ```dp``` -- positive real number  
#### Return value
Real number.

------
```
PriceBids()
```
#### Description
Returns all submitted prices in ```auctionset``` in ascending order.
#### Arguments
  ```auctionset``` -- vector, containing auction indeces
#### Return value
List of real numbers.

-----
```
FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
```
#### Description
```FOC()``` returns the value of $F^j$ (cf. the manuscript).
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
Real number.

-----
```
VarPhiU(q, v, vl)
```
#### Description
Determines $\varphi_u(q,v,v_l)$ (cf. the manuscript).
#### Arguments
  ```q``` -- positive real number  
  ```v``` -- marginal profit function  
  ```vl``` -- positive real number, corresponding to $v_l$  
#### Return value
Real number.

-----
```
VarPhiL(q, v, vu)
```
#### Description
Determines $\varphi_l(q,v,v_u)$ (cf. the manuscript).
#### Arguments
  ```q``` -- positive real number  
  ```v``` -- marginal profit function  
  ```vu``` -- positive real number, corresponding to $v_u$  
#### Return value
Real number.

