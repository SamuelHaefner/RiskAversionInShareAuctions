## constructs SLURM files for TestMon

jlscript = readlines("TestMonGeneric.jl")
shscript = readlines("TestMonGeneric.sh")

for run in [1:1:20;]
    outfile = string("TestMon",run,".jl")
    jlscript[10] = string("x = ",run)
    open(outfile,"w") do file
        for i in [1:1:length(jlscript);]
            println(file, jlscript[i])
        end
    end

    outfile = string("TestMon",run,".sh")
    shscript[4]=string("#SBATCH --job-name=TestMon",run)
    shscript[6]=string("julia TestMon",run,".jl")
    open(outfile,"w") do file
        for i in [1:1:length(shscript);]
            println(file, shscript[i])
        end
    end
end

mainshscript = []
push!(mainshscript,"#!/bin/bash")
for run in [1:1:20;]
    push!(mainshscript,string("sbatch TestMon",run,".sh"))
end
outfile = string("TestMon.sh")
open(outfile,"w") do file
    for i in [1:1:length(mainshscript);]
        println(file, mainshscript[i])
    end
end