#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load R/4.2
Rscript /home/perofskyamc/SFS_Rt_Block_Bootstrap/block_bootstrap_snowstorm_spearman.R > Rtry_block_bootstrap_snowstorm_sp.out