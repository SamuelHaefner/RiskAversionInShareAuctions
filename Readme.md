This file contains information about the Julia code used in:

# **[Risk Aversion in Share Auctions: Estimating Rents from TRQs in Switzerland](https://dx.doi.org/10.2139/ssrn.3397027)**

Author: [Samuel Häfner](https://samuelhaefner.github.io), University of St. Gallen, [samuel.haefner@unisg.ch](mailto:samuel.haefner@unisg.ch) 

------
## Overview

The julia files contain both the functions for estimation as well as the scripts to read and further process the estimates. The following main files contain all relevant functions and global variables (for more detailed information see below):

- ```Auxiliary.jl``` --- reads in the required packages, the data set and defines the required auxiliary functions and  global variables.  
- ```Grouping.jl``` --- determines the bidder groups and auction groups used for the estimation.  
- ```Estimation.jl``` --- contains the main functions for the estimation of W and the bounds.  
- ```TestMon.jl``` --- contains the main functions to test for the monotonicity of F in v.

The following files contain the generic scripts used to produce the scripts for the actual estimation (for more information, see the files themselves). The estimation was conducted at [sciCORE](http://scicore.unibas.ch/) scientific computing core facility at the University of Basel. The Julia version used was 1.2. Computation time depended on the script, ranging from 2-3 hours (WandTGeneric.jl and TestMonGeneric.jl) up to 48 − 72 hours (EstimateBoundsGeneric.jl).

- ```EstimateWandTGeneric.jl``` --- generic script to compute and save estimates of W and Theta.  
- ```EstimateWandTRobustGeneric.jl``` --- generic script to compute and save robustness checks for W and Theta, using a log-normal distribution rather than a gamma distribution when estimating W.  
- ```EstimateBoundsGeneric.jl``` --- generic script to compute standard and tighter bounds, using W from EstimateWandTGeneric.jl.   
- ```TestMonGeneric.jl``` --- generic script to test for monotonicity of F, determining W and the bounds for one bootstrap round, computing violations of monotonicity.

The following files contain the scripts to produce the .jl and .sh files that are needed reto run the estimates on a SLURM workload manager. For each *Generic.jl-File above, they split up the jobs into a manageable number of bootstrap rounds.  

- ```EstimateWandTSLURM.jl```  
- ```EstimateWandTRobustSLURM.jl```  
- ```EstimateBoundsSLURM.jl```  
- ```TestMonSLURM.jl```  

The following files contain the scripts used to read in and process the estimates produced with above *Slurm.jl scripts:

- ```ReadEstimates.jl``` --- contains the functions to read in the estimates and to produce the plot of the bounds.
- ```ReadT.jl``` --- script to read in estimates of T and produce the plots. 
- ```ReadTRobust.jl``` --- same as above, but using the alternative estimates (robustness check).
- ```ReadTestMon.jl``` --- script to read in the monotonicity violations obtained and saved in TestMonGeneric.jl.

The last file contains scripts used for the plots:

- ```ResamplingPlots.jl``` --- scripts to generate the group plots and the resampling plots. 

------

## Estimation.jl

### (a) Estimation of W:
```Wgamma()``` returns estimates of W(p,q) for all p in the vector ```prices```, estimating a gamma distribution:

```
Wgamma(prices, auctionset, bidderassignment, n, m, P) 
```

```prices``` -- vector of submitted prices to be used for the estimation of W.  
```auctionset``` -- vector with indeces of auctions to be used for the estimation of W.  
```bidderassignment``` -- bidder assignment vector  
```n``` -- vector, each entry corresponding to the (average) number of active bidders from a given group  
```m``` -- number of bootstrap rounds to be estimated  
```P``` -- number of rounds of resampling used for estimation


```Wlnorm()``` is same as ```Wgamma()```, yet estimating a log normal distribution:

```
Wlnorm(prices, auctionset, bidderassignment, n, m, P)
```

### (b) Computation of Simple Bounds:

```SimpleBound()``` returns a dataframe with simple bounds, v=[vub,vlb,q], based on ```bid``` and ```WPar```
```
SimpleBound(bid, WPar, rho, Q)
```

```bid``` -- bid function  
```WPar``` -- array, containing the estimated parameters of the distribution of (p_i^j) for steps j=1,...,k  
```rho``` -- positive real number, corresponding to the risk preference $\rho$  
```Q``` -- positive real number, corresponding to the quota $Q$

```EstimateSimpleBounds()```returns estimates of simple bounds for all bidders in ```auction``` using bids from ```auctionset``` assuming a risk preference ```rho```
```
EstimateSimpleBounds(auction, W, bidderassignment,  prices, rho, m)
```

```auction``` -- auction index  
```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
```bidderassignment``` -- bidder assignment vector  
```prices``` -- vector of submitted prices used for the estimation of W  
```rho``` -- positive real number, corresponding to the risk preference $\rho$  
```m``` -- number of bootstrap rounds to be estimated


### (c) Computation of Theta:

```EstTheta()``` returns share of inequality violations for a given auction:
```
EstTheta(auction, W, bidderassignment, prices, bounds, rho, m)
``` 

  ```auction``` -- auction index  
  ```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```   
  ```bidderassignment``` -- bidder assignment vector  
  ```prices``` -- vector of prices used for the estimation of W  
  ```bounds``` -- estimated simple bounds from ```EstimateSimpleBounds()```  
  ```rho``` -- positive real number, the risk preference   
  ```m``` -- number of bootstrap rounds  
  

### (d) Computation of Tighter Bounds:

```TighterBounds()``` returns tighter bounds from initial conditions ```initvl``` and ```initivu```, runing Algorithm 2

```
TighterBounds(bid, initvl, initvu, W, g, prices, rho, Q, bootstraprun, maxiter, tolerance)
```

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

```EstTighterBounds()``` returns tighter bounds for all bidders in an ```auction``` and all ```m``` bootstrap rounds

```
EstTighterBounds(auction, W, bidderassignment, prices, bounds, rho, m, maxiter, tolerance)
```


```auction``` -- auction index  
```W``` -- estimate of W, as returned from ```Wgamma()``` or ```Wlnorm()```  
 ```bidderassignment``` -- bidder assignment vector  
 ```prices``` -- vector of prices used for the estimation of ```W```  
 ```bounds``` -- estimated simple bounds from ```EstimateSimpleBounds()```  
```rho``` -- positive real number, the risk preference  
```m``` -- number of boostratp runs  
```maxiter``` -- maximum number of fixed point iterations  
```tolerance``` -- tolerance level (used in iteration)  

------

## TestMon.jl

### (a) The main function:

```TestIncreasingDiff()``` returns a list of length of the ```rhoindexset```, each entry is again a list of length of the ```bidderset```, and each entry of that list consists of two numbers: (1) the total number of bid steps tested for that bidder (2) the total number of bidsteps for which a violation of monotonicity was found.

```
TestIncreasingDiff(auction, bidderset, BoundsAuction, W, R, rhoindexset, rhovec)
```  

  ```auction``` -- auction index  
  ```bidderset``` -- bidder index set (set of bidders to be looked at)  
  ```BoundsAuction``` -- first element in list returned by ```LoadWandBounds()```  
  ```W``` -- second element in list returned by ```LoadWandBounds()```  
  ```R``` -- number of tests per bid step  
  ```rhoindexset``` -- vector of indexes, referring to rhovec, for which the test is run  
  ```rhovec``` -- vector of rho values  
 
### (b) The auxiliary functions:
```GetQVals(g, bidstep)``` returns the x values of the initial step function ```g``` and a random number of random values between the x-values of the ```bidstep```-th step and the ```bidstep```-1-th step
```
GetQVals(g, bidstep)
```  

  ```g``` -- decreasing step function function with a number of steps of at least *bidstep*  
  ```bidstep``` -- natural number, indicating the number of the step in the bid function  
  
```GetHigherF()``` returns a random step function between ```g``` and ```h``` with x values given in ```qval```
```
GetHigherF(g, h, qval)
```

  ```g```, ```h``` -- decreasing step function, satisfying ```g``` > ```h```    
  ```qvals``` -- x-values returned by GetQVals()

```CheckFOCMonotone()``` returns ```R``` differences in $F^j$, evaluated at the two different, ordered profit functions
```
CheckFOCMonotone(bidstep, bid, vub, vlb, W, g, prices, rho, Q, bootstraprun, R)
```

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

```LoadWandBounds()```  returns list containing W and simple bounds for ```auction```
```
LoadWandBounds(auction)
```
  ```auction``` -- auction index  
  
```TestIncreasingDiffBidder()``` returns two numbers: [number of tests performed (corresponds to the number of bid steps), number of violations]

```
TestIncreasingDiffBidder(auction, bidder, Bounds, W, R, rho)
```

  ```auction``` -- auction index  
  ```bidder``` -- bidder index  
  ```Bounds``` -- first element in list returned by ```LoadWandBounds(```)  
  ```W``` -- second element in list returned by ```LoadWandBounds()```  
  ```R``` -- number of tests per bid step  
  ```rho``` -- real number, risk preference $\rho$  
 
------

## Auxiliary.jl

### (a) Global variables:
Throughout the scripts I use the following global variables:

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

### (b) Auxiliary functions:

```qpBid(bidder, auction)``` returns the bid function of a bidder in an auction, which is a dataframe [qb (quantity points), pb (price points), cumqb (cumulated quantity points)]

```
qpBid(bidder, auction)
```

  ```bidder``` -- bidderindex  
  ```auction``` -- auctionindex  

```AvgNoBidders()``` returns the average number of bidders in each bidder cluster
```
AvgNoBidders(bidderassignment)
```
  
  ```bidderassignment``` -- vector, each element corresponding to a bidder, denoting the cluster number of the respective bidder  
  
```BidToCover()``` returns the bid-to-cover ratio in an auction
```
BidToCover(auction)
```

  ```auction``` -- auctionindex

```Revenue()``` returns the revenue in an auction
```
Revenue(auction)
```

  ```auction``` -- auctionindex
  
```qRec()``` returns the received quantity of a bidder in a given auction
```
qRec(bidder, auction)
```

  ```bidder``` -- bidderindex   
  ```auction``` -- auctionindex

```qShareRec()```  returns the received quantity share of a bidder in a given auction
```
qShareRec(bidder, auction)
```
  ```bidder``` -- bidderindex  
  ```auction``` -- auctionindex
  
```ActiveAuctions()``` returns the indeces in which a bidder is active
```
ActiveAuctions(bidder)
```
  ```bidder``` -- bidderindex 

```AvgBid()``` returns the average bid of a bidder
```
AvgBid(bidder)
```
  ```bidder``` -- bidderindex

```ShareSuccBidders()```  returns the  share of succesfull bidders in an auction
```
ShareSuccBidders(auction)
```
  ```auction``` -- auctionindex

```SuccessRate()``` returns the success rate of a bidder
```
SuccessRate(bidder)
```
  ```bidder``` -- bidderindex

```StepBid()``` returns the value of $\beta_b(q)$ (cf. the manuscript for a definition)
```
StepBid(q, bid)
```
  ```q``` -- positive real number
  ```bid``` -- bid function as returned from qpBid(bidder, auction)

```NoStepBid()```  returns the step number at ```q``` of the step function $\beta_b$
```
NoStepBid(q, bid)
```
  ```q``` positive real number
  ```bid``` bid function as returned from qpBid(bidder, auction)

```StepV()``` returns the value of v(q)
```
StepV(q, v)
```
  ```q``` -- positive real number
  ```v``` -- marginal profit function, which is a data frame with columns ```qval``` (quantities, need to be increasing) and ```vval``` (the values of v)

```IntBid()``` returns the value of $\int_a^b\beta_b(q)dq$
 ```
 IntBid(a, b, bid)
 ```
  ```a```, ```b``` -- positive real numbers
  ```bid``` -- bid function as returned from qpBid(bidder, auction)
  
```IntV()``` returns the value of $\int_a^b v(q)dq$
```
IntV(a, b, v)
```
  ```a```, ```b``` -- positive real numbers
  ```v``` -- marginal profit function

```PiOverline()``` returns the value of $\overline{\Pi}_i^j(b,v)$
```
PiOverline(bidstep, bid, v, WPar, rho, Q, n)
```
  ```bidstep``` -- positive natural number, corresponding to the number of the step under consideration, j  
  ```bid``` -- bid function  
  ```v``` -- marginal profit function  
  ```WPar``` -- array, containing the estimated parameters of the distribution of D(p_i^j) for steps j=1,...,k  
  ```rho``` -- positive real number, corresponding to the risk preference $\rho$  
  ```Q``` -- positive real number, corresponding to the quota $Q$  
  ```n``` -- positive natural number, indicating the number of support points used for integration  

```w()``` returns w(p,q)
```
w(q, WPar, dp)
```
  ```q``` -- positive real number  
  ```WPar``` -- array, containing the estimated parameters of the distribution of D(p) both for p_i^j and p_i^j + dp, j=1,...,k   
  ```dp``` -- positive real number  

```PriceBids(auctionset)``` returns all submitted prices in ```auctionset``` in ascending order
```
PriceBids()
```
  ```auctionset``` -- vector, containing auction indeces

```FOC()``` returns the value of the FOC, called $F^j$ in the text
```
FOC(bidstep, bid, v, W, group, prices, rho, Q, bootstraprun, n)
```
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

```VarPhiU()``` returns $\varphi_u(q,v,v_l)$ (in the text)
```
VarPhiU(q, v, vl)
```
  ```q``` -- positive real number  
  ```v``` -- marginal profit function  
  ```vl``` -- positive real number, corresponding to $v_l$  

```VarPhiL()``` returns $\varphi_l(q,v,v_u)$ (in the text)
```
VarPhiL(q, v, vu)
```
  ```q``` -- positive real number  
  ```v``` -- marginal profit function  
  ```vu``` -- positive real number, corresponding to $v_u$  

------

## Grouping.jl

This script define the following global variables used for the estimation of W:

```group```: list containing, for each auction group, a vector of the respective auction indeces  
```bidderassignment```: a vector determining for every bidder the bidder group assignment {1,2,3}

------