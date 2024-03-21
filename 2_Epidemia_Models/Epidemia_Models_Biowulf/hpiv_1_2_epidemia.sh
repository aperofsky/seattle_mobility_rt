#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load R/4.2
Rscript /home/perofskyamc/SFS_Rt_Epidemia/hpiv_1_2_rt_biowulf_epidemia.R > Rtry_hpiv_1_2_rt_biowulf_epidemia.out