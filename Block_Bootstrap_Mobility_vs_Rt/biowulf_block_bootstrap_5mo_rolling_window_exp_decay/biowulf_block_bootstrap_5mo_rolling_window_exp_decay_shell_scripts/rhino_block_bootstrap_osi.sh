#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load R/4.2
Rscript /home/perofskyamc/SFS_Rt_Block_Bootstrap/CCF_Rhino_Rt_vs_mobility_epidemia_block_bootstrap_5mo_window_Biowulf_osi_exp_decay.R > Rtry_biowulf_rhino_block_bootstrap_osi_exp_decay.out