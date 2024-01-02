#!/bin/sh
#
#SBATCH --job-name=simulate_trajectories
#SBATCH --time=3-23:59:00
#SBATCH --cpus-per-task=1
#SBATCH --qos=1week
#SBATCH --mem-per-cpu=8G

#

/scicore/home/neher/niralu56/.juliaup/bin/julia trajectories_random_beta.jl