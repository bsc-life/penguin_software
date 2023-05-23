#!/bin/bash
#SBATCH --job-name="Betweness1"                                                                               
#SBATCH -D .
#SBATCH --output=Betweness1%j.out                                                                                       
#SBATCH --error=Betweness1%j.err
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=07:30:00
#SBATCH --qos=bsc_ls

module load R
Rscript Betweness_significance_analysis_cluster_4.R

