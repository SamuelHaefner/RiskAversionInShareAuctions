#!/bin/bash
#SBATCH --qos=1week
#SBATCH --mem-per-cpu=48G
#SBATCH --job-name=ETB1
module load Julia
julia EstimateBounds1.jl
