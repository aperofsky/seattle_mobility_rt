## `4_Block_Bootstrap_Mobility_vs_Rt` folder: Estimate time series cross-correlations between mobility and Rt

* `1_snowstorm_rolling_window_block_bootstrap` folder
    *   Estimate cross-correlations between mobility and Rt during winter 2018-2019
    *   Analyses performing block bootstrap cross-correlations were run on the NIH Biowulf Cluster, with most jobs completing within a few hours. They can be run on a local machine but each individual analysis could take several hours.
    1.  Transfer `block_bootstrap_snowstorm_spearman.R`, `combined_rt_mobility_15day_mv_avg.rds`, `utils.R`, and `snowstorm_sp.sh` files to the `SFS_Rt_Block_Bootstrap` folder in your `home` directory on the cluster.
    2.  Login to Biowulf and change your current directory to the `SFS_Rt_Block_Bootstrap` folder.
    3.  Run the following sbatch command:
```
sbatch --gres=lscratch:500 --time=1-00:00:00 --mail-type=BEGIN,FAIL,TIME_LIMIT_90,END --exclusive --mem=60g snowstorm_sp.sh
```
        Outputs will be saved to your `data` directory on Biowulf. 
    4.  Transfer outputs to `biowulf_block_bootstrap_1mo_rolling_window_output_spearman` subfolder on your local machine.

* `2_post_2019_rolling_window_block_bootstrap` folder
    *   Estimate cross-correlations between mobility and Rt from Fall 2019 to Summer 2022
    1.  `create_swarm_5mo_block_boostrap.R`: Create swarm script to submit several jobs at one time (`Run_5mo_rolling_block_bootstrap_swarm.txt`). Individual `sbatch` commands are located `block_bootstrap_sbatch_commands_spearman.txt` if you'd like to run cross-correlation analyses individually.
    2.  Transfer `Run_5mo_rolling_block_bootstrap_swarm.txt`, all the scripts within the `BB_5mo_rolling_window_R_scripts_spearman` folder, and `combined_rt_mobility_15day_mv_avg.rds` to the `SFS_Rt_Block_Bootstrap` folder in your `home` directory on the cluster.
    3.  Login to Biowulf and change your current directory to the `SFS_Rt_Block_Bootstrap` folder.
    4.  Run the following `swarm` command:
```
swarm -f Run_5mo_rolling_block_bootstrap_swarm.txt -g 5 --gres=lscratch:500 --time=12:00:00 --module R/4.2
```
        Outputs will be saved to your `data` directory on Biowulf. Most jobs are completed within 2 hours using `swarm`. Jobs submitted via individual `sbatch` commands take less time (less than 30 minutes).
    5.   Transfer outputs to `BB_5mo_rolling_window_output_spearman` subfolder on your local machine.
    
* `3_make_block_bootstrap_figures` folder
    *   Figures showing cross-correlations and optimal lags averaged by calendar month:
        *   `1_CCF_Rt_vs_mobility_epidemia_snowstorm_by_month.R` (Figure S6)
        *   `2_1_CCF_Rt_vs_mobility_epidemia_2019_2020_by_month_all_indicators.R` (Figure S10)
        *   `2_2_CCF_Rt_vs_mobility_epidemia_2019_2020_by_month_select_indicators.R` (Figure 4)

    *   Figures showing daily or weekly rolling cross-correlations:
        *   `3_Feb_2019_snowstorm_rolling_block_bootstrap_figure.R`: Figure of rolling daily cross-correlations between Rt and mobility during winter 2018-2019.
        *   `4_compile_5mo_block_bootstrap_results.R`: Compile rolling weekly cross-correlation estimates for Fall 2019 to Summer 2022.
        *   `5_2019_2020_block_bootstrap_figures_and_cross_correlations.R`: Figure of rolling weekly cross-correlations between endemic virus Rt and mobility during winter 2019-2020, before the start of the COVID-19 pandemic.
        *   `6_COVID_wave_block_bootstrap_figures_and_cross_correlations.R`: Figure of rolling weekly cross-correlations between SARS-CoV-2 Rt and mobility, masking, and NPI stringency across four individual COVID-19 waves (winter 2020-2021, Alpha wave, Delta wave, Omicron BA.1 wave).
        * `7_endemic_virus_rebound_block_bootstrap_figures_and_cross_correlations.R`: Figures of rolling weekly cross-correlations between endemic virus Rt and mobility during 2020 (non-enveloped viruses) and 2021 (enveloped viruses)
