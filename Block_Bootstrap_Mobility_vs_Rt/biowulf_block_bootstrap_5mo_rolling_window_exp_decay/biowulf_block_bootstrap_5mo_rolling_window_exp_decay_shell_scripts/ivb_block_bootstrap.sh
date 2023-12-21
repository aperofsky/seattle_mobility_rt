#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load R/4.2
Rscript /home/perofskyamc/SFS_Rt_Block_Bootstrap/CCF_IVB_Rt_vs_mobility_epidemia_block_bootstrap_5mo_window_Biowulf_exp_decay.R > Rtry_biowulf_ivb_block_bootstrap_exp_decay.out