## Estimate time series cross-correlations between pathogen Rt and population mobility

**Note: Sections 1 and 2 are computationally intensive and should be run on a HPC cluster.**

`1_snowstorm_rolling_window_block_bootstrap/`
*   Estimate daily rolling cross-correlations between endemic virus Rt and mobility during winter 2018-2019, using the NIH Biowulf cluster.
*   Transfer the following files to a folder named `SFS_Rt_Block_Bootstrap/` in your `home` directory: 
    *   `block_bootstrap_snowstorm_spearman.R`
    *   `combined_rt_mobility_15day_mv_avg.rds`
    *   `utils.R`
    *   `snowstorm_sp.sh`
*   `ssh` into Biowulf and change your current directory to `SFS_Rt_Block_Bootstrap/`.
*   Run the following `sbatch` command:
```
sbatch --gres=lscratch:500 --time=1-00:00:00 --mail-type=BEGIN,FAIL,TIME_LIMIT_90,END --exclusive --mem=60g snowstorm_sp.sh
```
*   Outputs will be saved to your `data` directory on Biowulf.
*   When the job is complete, transfer outputs to `biowulf_block_bootstrap_1mo_rolling_window_output_spearman/` on your local machine.
*   Runtime: The job can take several hours to run on the cluster because the cross-correlations for each pathogen are run sequentially. The code could be reworked so that there are separate jobs for each pathogen, either via individual `sbatch` commands or a `swarm` script.

`2_post_2019_rolling_window_block_bootstrap/`
*   Estimate cross-correlations between Rt and mobility during Fall 2019 to Summer 2022, , using the NIH Biowulf cluster.
*   `create_swarm_5mo_block_boostrap.R`
    *   Create `swarm` script  `Run_5mo_rolling_block_bootstrap_swarm.txt` to submit several jobs at one time to the cluster. `Run_5mo_rolling_block_bootstrap_swarm.txt` is already provided but this R script shows how to produce a `swarm` script for submitting many similar jobs simultaneously.
*   Transfer the following files to the `SFS_Rt_Block_Bootstrap` folder in your `home` directory on the cluster. 
    *   `Run_5mo_rolling_block_bootstrap_swarm.txt`
    *   All the scripts within `BB_5mo_rolling_window_R_scripts_spearman/` (transfer individual files and not the folder itself)
    *   `combined_rt_mobility_15day_mv_avg.rds` 
*   `ssh` into Biowulf and change your current directory to `SFS_Rt_Block_Bootstrap/`.
*   Run the following `swarm` command:
```
swarm -f Run_5mo_rolling_block_bootstrap_swarm.txt -g 5 --gres=lscratch:500 --time=12:00:00 --module R/4.2
```

*   Outputs will be saved to your `data` directory on Biowulf.
*   If you'd like to run cross-correlation analyses individually, `sbatch` commands are listed in `block_bootstrap_sbatch_commands_spearman.txt`. To run `sbatch` commands, copy/paste them into your terminal window.   
*   After jobs are complete, transfer outputs to `BB_5mo_rolling_window_output_spearman/` on your local machine.
*   Runtime: Most jobs complete within 2 hours using `swarm`. Jobs submitted via individual `sbatch` commands take less time (< 30 minutes). The above `swarm` command's memory and node specifications could be optimized to make the jobs run faster.

`3_make_block_bootstrap_figures/`
*   Figures showing cross-correlations and optimal lags averaged by calendar month:
    *   `1_CCF_Rt_vs_mobility_epidemia_snowstorm_by_month.R`
        *   Figure S6: Daily cross-correlations between endemic virus Rt and mobility during winter 2018-2019, averaged by month.
    *   `2_1_CCF_Rt_vs_mobility_epidemia_2019_2020_by_month_all_indicators.R`
        *   Figure S10: Weekly cross-correlations between Rt and mobility during winter 2019-2020, averaged by month. Shows all mobility metrics considered in the study.
    *   `2_2_CCF_Rt_vs_mobility_epidemia_2019_2020_by_month_select_indicators.R`
        *   Figure 4: Same as Figure S10 but shows only a subset of mobility metrics that strongly correlate with endemic virus Rt during Fall 2019.
    *   Runtime for each script: < 45 seconds

*   Figures showing daily or weekly rolling cross-correlations:
    *   `3_Feb_2019_snowstorm_rolling_block_bootstrap_figure.R`
        *   Figure S7: Daily rolling cross-correlations between Rt and mobility during winter 2018-2019.
        *   Runtime: < 15 seconds
    *   `4_compile_5mo_block_bootstrap_results.R`
        *   Compile rolling cross-correlation estimates for weeks spanning Fall 2019 to Summer 2022.
        *   Runtime: < 1 second
    *   `5_2019_2020_block_bootstrap_figures_and_cross_correlations.R`
        *   Figure S8: Weekly rolling cross-correlations between Rt and mobility during winter 2019 - 2020.
        *   Runtime: < 30 seconds
    *   `6_COVID_wave_block_bootstrap_figures_and_cross_correlations.R`
        *   Figure 5: Weekly rolling cross-correlations between SARS-CoV-2 Rt and mobility, during Spring 2020 to Summer 2022.
        *   Runtime: < 10 seconds
    *   `7_endemic_virus_rebound_block_bootstrap_figures_and_cross_correlations.R`
        *   Weekly cross-correlations between Rt and mobility during endemic virus rebound in 2020 (Figure S13: non-enveloped viruses) and 2021 (Figure S15: enveloped viruses).
        *   Runtime: < 15 seconds
