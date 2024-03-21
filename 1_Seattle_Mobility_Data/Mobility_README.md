##  `1_Seattle_Mobility_Data` folder: Create cell phone mobility datasets and combine with data on masking and NPI stringency

* `1_import_SafeGraph_Data.R`:
    *   Note: Input and output data cannot be shared but you can run this code if you have access to the raw SafeGraph data.
    *   Import, organize, and clean raw SafeGraph foot traffic data
    *   Join SafeGraph data to census data (e.g., population sizes for census block groups and King County)
    *   Create foot traffic mobility indices: daily visits to various categories of points of interest (POIs), weekly visits aggregated by visitor home census block group (CBG)
    
* `2_process_Safegraph_data_and_city_maps.R`:
    *   Note: Input and output data cannot be shared but you can run this code if you have access to the raw SafeGraph data.
    *   Combine SafeGraph Weekly Patterns dataset with Home Panel and Visitor Panel datasets. Adjust for changes in panel size over time.
    *   Plot Seattle mobility network and network degree distributions at different time points (Figure 3)
    
* `3_external_datasets_behavior`:
    *   Import and combine publicly accessible data sources for % staying home, % masking, and NPI stringency
    *   Combine % staying home metrics from SafeGraph and Meta Data for Good to create custom metric (Figure S24)

* `4_mobility_data_df.R`:
    *   Note: requires outputs from `1_import_SafeGraph_Data.R` and `2_process_Safegraph_data_and_city_maps.R` but these can't be shared.
    *   Combine SafeGraph Weekly Patterns mobility metrics with other data sources on mobility, masking, and NPI stringency (saved as `mobility_metrics_for_epidemia.rds` in `mobility_data` folder)

* `5_manuscript_mobility_figures`:
    *   **Note: Start here if you do not have access to the raw SafeGraph data.**
    *   Figure S4: % change in baseline for individual POI categories
    *   Figure 2: % change in baseline for large scale movements and different POI categories; % staying home and % masking
    *   Figure S17: Same as Figure 2 but focused on time period of Omicron BA.1 wave during winter 2021-2022
