## constructs SLURM files for EstimateBounds

jlscript = readlines("EstimateBoundsGeneric.jl")
shscript = readlines("EstimateBoundsGeneric.sh")

for run in [1,2,3,4,5,6,35,36,37,38,39,40]
    outfile = string("EstimateBounds_add",run,".jl")
    jlscript[12] = string("@load \"WValues",run,".dat\" W",run)
    jlscript[13] = string("W=W",run)
    jlscript[18] = string("for g in [3]")
    jlscript[57] = string("Bounds_add",run," = Bounds")
    jlscript[58] = string("@save \"Bounds_add",run,".dat\" Bounds_add",run)
    open(outfile,"w") do file
        for i in [1:1:length(jlscript);]
            println(file, jlscript[i])
        end
    end

    outfile = string("EstimateBounds_add",run,".sh")
    shscript[4]=string("#SBATCH --job-name=ETB",run)
    shscript[6]=string("julia EstimateBounds_add",run,".jl")
    open(outfile,"w") do file
        for i in [1:1:length(shscript);]
            println(file, shscript[i])
        end
    end
end

mainshscript = []
push!(mainshscript,"#!/bin/bash")
for run in [1,2,3,4,5,6,35,36,37,38,39,40]
    push!(mainshscript,string("sbatch EstimateBounds_add",run,".sh"))
end
outfile = string("EstimateBounds_add.sh")
open(outfile,"w") do file
    for i in [1:1:length(mainshscript);]
        println(file, mainshscript[i])
    end
end

