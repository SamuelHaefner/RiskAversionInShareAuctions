#!/bin/bash
#SBATCH --qos=1day
#SBATCH --job-name=TestMon
module load Julia
julia TestMon1.jl
