#!/bin/bash
#SBATCH --job-name="Betweness"                                                                               
#SBATCH -D .
#SBATCH --output=Betweness%j.out                                                                                       
#SBATCH --error=Betweness%j.err
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=07:30:00
#SBATCH --qos=bsc_ls

module load R
Rscript Betweness_significance_analysis_cluster_2.R

