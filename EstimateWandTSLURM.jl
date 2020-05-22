## constructs SLURM files for EstimateWandT

jlscript = readlines("EstimateWandTGeneric.jl")
shscript = readlines("EstimateWandTGeneric.sh")

for run in [1:1:40;]
    outfile = string("EstimateWandT",run,".jl")
    jlscript[17] = string("W",run," = W")
    jlscript[18] = string("@save \"WValues",run,".dat\" W",run)
    jlscript[57] = string("Theta",run," = Theta")
    jlscript[58] = string("@save \"Theta",run,".dat\" Theta",run)
    open(outfile,"w") do file
        for i in [1:1:length(jlscript);]
            println(file, jlscript[i])
        end
    end

    outfile = string("EstimateWandT",run,".sh")
    shscript[3]=string("#SBATCH --job-name=EWandT",run)
    shscript[5]=string("julia EstimateWandT",run,".jl")
    open(outfile,"w") do file
        for i in [1:1:length(shscript);]
            println(file, shscript[i])
        end
    end
end

mainshscript = []
push!(mainshscript,"#!/bin/bash")
for run in [1:1:40;]
    push!(mainshscript,string("sbatch EstimateWandT",run,".sh"))
end
outfile = string("EstimateWandT.sh")
open(outfile,"w") do file
    for i in [1:1:length(mainshscript);]
        println(file, mainshscript[i])
    end
end

