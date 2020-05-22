#!/bin/bash
#SBATCH --qos=1day
#SBATCH --job-name=EWandT1
module load Julia
julia EstimateWandTGeneric.jl
