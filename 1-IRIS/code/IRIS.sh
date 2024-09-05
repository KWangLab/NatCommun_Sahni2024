#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=120g
#SBATCH --gres=lscratch:20
#SBATCH --partition=norm,ccr
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,FAIL
#SBATCH --cpus-per-task=40

module load R
Rscript /data/sahnis2/final_antipdl1_project/IRIS/IRIS_final.R
