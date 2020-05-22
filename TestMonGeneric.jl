include("Auxiliary.jl")
include("Estimation.jl")
include("Grouping.jl")
include("TestMon.jl")

n = AvgNoBidders(bidderassignment)
rhovec = [exp(x) for x in [-10:0.5:0;]]

# x to be adapted, x in [1:1:length(rhovec);]!
x = 1
rho = rhovec[x]
for R in [50,100,200]
    for auction in [1:1:39;]
        # requires WValuesTestMon.dat and BoundsTestMon[auction].dat in working directory
        global WB = LoadWandBounds(auction)
        global T = TestIncreasingDiff(auction, activebidderindeces[auction], WB[1], WB[2], R, 1, rho)
        eval(Meta.parse(string("@save \"TestMonA", auction, "R", R, "rho", x, ".dat\" T")))
    end
end
