## constructs SLURM files for EstimateBounds

jlscript = readlines("EstimateBoundsGeneric.jl")
shscript = readlines("EstimateBoundsGeneric.sh")

for run in [1:1:40;]
    outfile = string("EstimateBounds",run,".jl")
    jlscript[12] = string("@load \"WValues",run,".dat\" W",run)
    jlscript[13] = string("W=W",run)
    jlscript[57] = string("Bounds",run," = Bounds")
    jlscript[58] = string("@save \"Bounds",run,".dat\" Bounds",run)
    open(outfile,"w") do file
        for i in [1:1:length(jlscript);]
            println(file, jlscript[i])
        end
    end

    outfile = string("EstimateBounds",run,".sh")
    shscript[4]=string("#SBATCH --job-name=ETB",run)
    shscript[6]=string("julia EstimateBounds",run,".jl")
    open(outfile,"w") do file
        for i in [1:1:length(shscript);]
            println(file, shscript[i])
        end
    end
end

mainshscript = []
push!(mainshscript,"#!/bin/bash")
for run in [1:1:40;]
    push!(mainshscript,string("sbatch EstimateBounds",run,".sh"))
end
outfile = string("EstimateBounds.sh")
open(outfile,"w") do file
    for i in [1:1:length(mainshscript);]
        println(file, mainshscript[i])
    end
end

