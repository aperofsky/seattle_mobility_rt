## Estimate pathogen daily effective reproduction numbers (Rt) using semi-mechanistic epidemiological models

* `Epidemia_Models_Biowulf` subfolder: 
    *   Contains R scripts, shell scripts, and various inputs for fitting transmission models using the `Epidemia` R package. `epidemia_model_sbatch_commands.txt` includes individual batch scripts to run transmission models on the NIH Biowulf Cluster. 
        - `sbatch` commands and shell scripts (files ending with `.sh`) are used to specify the resources you need to run your jobs, such as the number of nodes you want to run your jobs on and how much memory you’ll need. Slurm then schedules your jobs based on the availability of the resources you’ve specified. 
        - These models can be run on a local machine but they are computationally intensive and can take several hours to several days, depending on the pathogen and time period. 
        - To run these models on the NIH Biowulf cluster, copy/paste all the files into one folder located in your `home` directory on the cluster (here my folder is called `SFS_Rt_Epidemia`). Login to the cluster, change your current directory to `SFS_Rt_Epidemia`, and copy/paste the individual sbatch commands into your terminal window. Rt estimates will be saved to your `data` directory on the cluster.
        - After models have finished running, transfer Rt estimates in the `data` directory on Biowulf to the `epidemia_rt_output` subfolder.
            
* `1_epidemia_Rt_figures.R`:
    *   Combine Rt estimates from all pathogens
    *   Compare time series of SFS/SCAN COVID incidence to confirmed King County cases (Figure S27)
    *   Plot pathogen Rt and incidences during winter 2018 - 2018 (Figure S5)
    *   Plot pathogen Rt and incidences for whole study period (Figure 1)

* `2_change_in_rt_nonparametric_bootstrap.R`:
    *   Compare Rt before and after 2 major events: Feb 2019 snowstorm and Feb 29, 2020 State of Emergency Declaration, using non-parametric bootstrap tests for the ratio of two means
