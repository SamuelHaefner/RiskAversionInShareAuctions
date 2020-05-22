## constructs SLURM files for EstimateWandTRobust

jlscript = readlines("EstimateWandTRobustGeneric.jl")
shscript = readlines("EstimateWandTRobustGeneric.sh")

for run in [1:1:40;]
    outfile = string("EstimateWandTRobust",run,".jl")
    jlscript[17] = string("W",run," = W")
    jlscript[18] = string("@save \"WValuesRobust",run,".dat\" W",run)
    jlscript[57] = string("Theta",run," = Theta")
    jlscript[58] = string("@save \"ThetaRobust",run,".dat\" Theta",run)
    open(outfile,"w") do file
        for i in [1:1:length(jlscript);]
            println(file, jlscript[i])
        end
    end

    outfile = string("EstimateWandTRobust",run,".sh")
    shscript[3]=string("#SBATCH --job-name=EWandTRob",run)
    shscript[5]=string("julia EstimateWandTRobust",run,".jl")
    open(outfile,"w") do file
        for i in [1:1:length(shscript);]
            println(file, shscript[i])
        end
    end
end

mainshscript = []
push!(mainshscript,"#!/bin/bash")
for run in [1:1:40;]
    push!(mainshscript,string("sbatch EstimateWandTRobust",run,".sh"))
end
outfile = string("EstimateWandTRobust.sh")
open(outfile,"w") do file
    for i in [1:1:length(mainshscript);]
        println(file, mainshscript[i])
    end
end

