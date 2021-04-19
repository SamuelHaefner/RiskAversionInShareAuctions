#!/bin/bash
#SBATCH --qos=1day
#SBATCH --mem-per-cpu=48G
#SBATCH --job-name=TestMon
module load Julia
julia TestMon1.jl
