library(dplyr)
library(stringr)

##### create swarm file to run 5 month rolling block bootstraps on NIH Biowulf cluster
## Filenames for R scripts
file_names <- list.files("3_Block_Bootstrap_Mobility_vs_Rt/post_2019_rolling_window_block_bootstrap/BB_5mo_rolling_window_R_scripts_spearman/")
file_names
fc_files <- file_names[1:26]
fc_files
length(fc_files)

## create swarm file by looping through individual R scripts
file_list <- list()
for (i in fc_files) {
  ## pathogen name
  first <- str_split_1(i, "CCF_")
  second <- str_split_1(first[2], "_Rt")
  pathogen_name <- second[1]

  # R output file name (note: not outputs from analysis)
  ## is there masking, osi, or rebound in the name?
  masking <- stringr::str_split(i, "_")[[1]][stringr::str_split(i, "_")[[1]] == "masking"]
  osi <- stringr::str_split(i, "_")[[1]][stringr::str_split(i, "_")[[1]] == "osi"]
  rebound <- stringr::str_split(i, "_")[[1]][stringr::str_split(i, "_")[[1]] == "rebound"]
  
  ## some pathogens have different R scripts for pre- and post-pandemic dynamics
  ## some pathogens have more than one R script due to different time frames of availability for behavioral indicators (e.g., masking, osi)
  full_name <- if_else(length(masking) == 0, pathogen_name, paste(pathogen_name, masking, sep = "_"))
  full_name <- if_else(length(osi) == 0, full_name, paste(full_name, osi, sep = "_"))
  full_name <- if_else(length(rebound) == 0, full_name, paste(full_name, rebound, sep = "_"))

  ## folder in home directory on Biowulf
  temp <- paste("Rscript /home/perofskyamc/SFS_Rt_Block_Bootstrap/", i, " > ", "Rtry_", full_name, ".out", sep = "")
  file_list[[i]] <- temp
}
file_list
swarm <- do.call("rbind", file_list)
head(swarm)
write.table(swarm, row.names = F, col.names = F, file = "3_Block_Bootstrap_Mobility_vs_Rt/post_2019_rolling_window_block_bootstrap/Run_5mo_rolling_block_bootstrap_swarm.txt", quote = FALSE)
