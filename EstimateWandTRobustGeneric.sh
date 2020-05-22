#!/bin/bash
#SBATCH --qos=1day
#SBATCH --job-name=EWandTRob1
module load Julia
julia EstimateWandTRobustGeneric.jl
