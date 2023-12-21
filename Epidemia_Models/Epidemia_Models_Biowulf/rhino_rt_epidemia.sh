#!/bin/bash
export TMPDIR=/lscratch/$SLURM_JOB_ID
module load R/4.2
Rscript /home/perofskyamc/SFS_Rt_Epidemia/rhino_rt_biowulf_epidemia.R > Rtry_rhino_rt_biowulf_epidemia.out