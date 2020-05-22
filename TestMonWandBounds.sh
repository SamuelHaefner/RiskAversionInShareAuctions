#!/bin/bash
#SBATCH --qos=1day
#SBATCH --job-name=TestMonWB
module load Julia
julia TestMonWandBounds.jl
