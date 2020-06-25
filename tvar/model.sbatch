#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
module load matlab/96
matlab -nodesktop < /scratch/kais/rstarBrookings2017/tvar/MainModel1.m
