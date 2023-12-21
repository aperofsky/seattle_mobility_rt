#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load R/4.2
Rscript /home/perofskyamc/SFS_Rt_Epidemia/cov_HKU1_NL63_rt_biowulf_epidemia.R > Rtry_cov_HKU1_NL63_rt_biowulf_epidemia.out