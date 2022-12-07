## this script produces the data overview table (Table 1)

include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")

using Plots
using LaTeXStrings
using Latexify

#quotas: quotas (from Auxiliary.jl)

#bid-to-cover ratio
btcratio = [BidToCover(x) for x in [1:length(auctionindeces);]]

#number of bidders: activebidders (from Auxiliary.jl)

#market clearing price: clearingprices (from Auxiliary.jl)

#revenue
revenue = [Revenue(x) for x in [1:length(auctionindeces);]]


#average total quanitity bid per bidder
avbidperbidder = [AvgBid(x) for x in bidderindeces]


#share of successful bidders
succbidders = [ShareSuccBidders(x) for x in [1:length(auctionindeces);]]


#success rate per bidder
succrate = [SuccessRate(x) for x in bidderindeces]


#success rate per bidder*
succratestar = []
for i in bidderindeces
    if length(ActiveAuctions(i)) >= 35
        push!(succratestar,succrate[i])
    end
end

#share of allocated quantity (if positive)
sharealloc = []
for i in bidderindeces
    for j in [1:length(auctionindeces);]
        push!(sharealloc,qShareRec(i,j))
    end
end
sharealloc = sharealloc[sharealloc.>0]

# construct rows of table
quotas_row = [mean(quotas),minimum(quotas),quantile(quotas,0.25),quantile(quotas,0.75),maximum(quotas)]./1000
btcratio_row = [mean(btcratio),minimum(btcratio),quantile(btcratio,0.25),quantile(btcratio,0.75),maximum(btcratio)]
activebidders_row = [mean(activebidders),minimum(activebidders),quantile(activebidders,0.25),quantile(activebidders,0.75),maximum(activebidders)]
clearingprices_row = [mean(clearingprices),minimum(clearingprices),quantile(clearingprices,0.25),quantile(clearingprices,0.75),maximum(clearingprices)]./100
revenue_row = [mean(revenue),minimum(revenue),quantile(revenue,0.25),quantile(revenue,0.75),maximum(revenue)]./1000000
avbidperbidder_row = [mean(avbidperbidder),minimum(avbidperbidder),quantile(avbidperbidder,0.25),quantile(avbidperbidder,0.75),maximum(avbidperbidder)]
succbidders_row = [mean(succbidders),minimum(succbidders),quantile(succbidders,0.25),quantile(succbidders,0.75),maximum(succbidders)]
succrate_row = [mean(succrate),minimum(succrate),quantile(succrate,0.25),quantile(succrate,0.75),maximum(succrate)]
succratestar_row = [mean(succratestar),minimum(succratestar),quantile(succratestar,0.25),quantile(succratestar,0.75),maximum(succratestar)]
sharealloc_row = [mean(sharealloc),minimum(sharealloc),quantile(sharealloc,0.25),quantile(sharealloc,0.75),maximum(sharealloc)]./100


# create dataframe
TableData = DataFrame(hcat(
    quotas_row,
    btcratio_row,
    activebidders_row,
    clearingprices_row,
    revenue_row,
    avbidperbidder_row,
    succbidders_row,
    succrate_row,
    succratestar_row,
    sharealloc_row
)',:auto)
rename!(TableData,[:Mean,:Min,:Pctl25,:Pctl75,:Max])


## save as latex table
t=latexify(
    TableData,
    env = :tabular,
)
write("Table1.txt",t)