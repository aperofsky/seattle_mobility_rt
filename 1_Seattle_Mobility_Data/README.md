## Create cell phone mobility datasets and combine with data on masking and the stringency of COVID-19 non-pharmaceutical interventions (NPIs)

* `1_import_SafeGraph_Data.R`
    *   Note: Input and output data cannot be shared but you can run this code if you have access to the raw SafeGraph data.
    *   Import, organize, and clean raw SafeGraph foot traffic data.
    *   Join SafeGraph data to census data (e.g., population sizes for census block groups and King County).
    *   Create foot traffic mobility indices: daily visits to various categories of points of interest (POIs), weekly visits aggregated by visitor home census block groups (CBGs).
    
* `2_process_Safegraph_data_and_city_maps.R`
    *   Note: Input and output data cannot be shared but you can run this code if you have access to the raw SafeGraph data.
    *   Combine SafeGraph Weekly Patterns dataset with Home Panel and Visitor Panel datasets. Adjust for changes in panel size over time.
    *   Figure 3: Seattle mobility network and network degree distributions at key epidemiological time points.
    
* `3_external_datasets_behavior`
    *   Import and combine publicly accessible data sources for % devices staying home, % masking in public, and the stringency of NPIs.
    *   Figure S24: Combine data from SafeGraph and Meta Data for Good to create a custom % devices staying at home metric.

* `4_mobility_data_df.R`
    *   Note: Requires outputs from `1_import_SafeGraph_Data.R` and `2_process_Safegraph_data_and_city_maps.R` but these data cannot be shared.
    *   Combine SafeGraph Weekly Patterns mobility metrics with other data sources on mobility, masking, and NPI stringency.
    *   Estimate the percentage change in baseline activity for SafeGraph foot traffic metrics, with the baseline for each metric defined as the mean foot traffic volume during 2019 (excluding national holidays).
    *   Processed mobility metrics are saved as `mobility_metrics_for_epidemia.rds` in the `mobility_data/` subfolder.

* `5_manuscript_mobility_figures`
    *   **Note: Start here if you do not have access to the raw SafeGraph data.**
    *   Figure S4: % change in baseline for visits to individual POI categories.
    *   Figure 2: A. % change in baseline for large scale movements. B. % change in baseline for visits to various POI categories. C. % devices staying home and % individuals masking in public.
    *   Figure S17: Same as Figure 2 but focused on the Omicron BA.1 wave during winter 2021-2022.
