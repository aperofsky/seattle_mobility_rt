#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load R/4.2
Rscript /home/perofskyamc/SFS_Rt_Block_Bootstrap/CCF_H1N1_Rt_vs_mobility_epidemia_block_bootstrap_5mo_window_Biowulf_exp_decay.R > Rtry_biowulf_h1n1_block_bootstrap_exp_decay.out