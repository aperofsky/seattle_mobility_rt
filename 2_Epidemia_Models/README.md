## Estimate respiratory virus daily effective reproduction numbers (Rt) using semi-mechanistic epidemiological models

* `Epidemia_Models_Biowulf/`
    *   This folder contains R scripts, shell scripts, and various inputs for fitting transmission models using the `Epidemia` R package.
    *   The `epidemia_model_sbatch_commands.txt` file includes individual batch scripts to run transmission models on the NIH Biowulf Cluster.
        -   `sbatch` commands and shell scripts (`.sh`) specify the resources you need to run your jobs, such as the number of nodes you want to run your jobs on and how much memory you’ll need. Slurm then schedules your jobs based on the availability of the resources you’ve specified. 
        -   These models can be run on a local machine but they will take several hours to several days, depending on the pathogen and time period. Although these models also take a similar amount of time to run on the cluster, they can be run simultaneously and you can use your local machine to do other things while you wait for the jobs to finish. :)
    *   To run transmission models on the NIH Biowulf cluster, transfer all the files in `Epidemia_Models_Biowulf/` to a folder named `SFS_Rt_Epidemia/` in your `home` directory on the cluster.
    *   In your terminal window, `ssh` into the cluster, change your current directory to `SFS_Rt_Epidemia/`, and copy/paste the individual sbatch commands to submit jobs to the cluster.
    *   Rt estimates will be saved to your `data` directory on the cluster.
    *   Transfer the Rt estimates saved in your `data` directory on Biowulf to `Epidemia_Models_Biowulf/epidemia_rt_output/` on your local machine.
            
* `1_epidemia_Rt_figures.R`
    *   Compile Rt estimates from all pathogens and time periods into one data frame.
    *   Figure S27: Compare the time series of Seattle Flu Study COVID-19 incidence to the time series of confirmed COVID-19 cases in King County.
    *   Figure S5: Pathogen Rt estimates and incidences during winter 2018 - 2019.
    *   Figure 1: Pathogen Rt estimates and incidences during the whole study period, late 2018 - summer 2022.

* `2_change_in_rt_nonparametric_bootstrap.R`
    *   Use non-parametric bootstrap tests for the ratio of two means to compare changes in Rt before and after 2 major events:
        -   Major snowstorm in February 2019
        -   COVID-19 State of Emergency Declaration in late February 2020