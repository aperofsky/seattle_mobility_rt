## Short-term forecasts of rhinovirus, adenovirus, and SARS-CoV-2 Rt

### Model inputs
`1_Climate_data/`
*   `noaa_local_climate_data.R`
    *   Daily records of precipitation, average wet bulb temperature, and average relative humidity are publicly accessible through the National Centers for Environmental Information’s [U.S. Local Climatological Database](https://www.ncei.noaa.gov/products/land-based-station/local-climatological-data).
    *   Import, clean, and interpolate missing values for meteorological data collected in Seattle, Washington. 
    *   Runtime: < 1 second

`2_Vaccination_coverage_data/`
*   Daily COVID-19 vaccination coverage in King County, Washington. Data are publicly accessible through the [Public Health – Seattle & King County COVID-19 Vaccination dashboard](https://kingcounty.gov/en/dept/dph/health-safety/disease-illness/covid-19/data/vaccination).

`3_VOC_data/`
*   Nextstrain-curated SARS-CoV-2 sequence metadata from GISAID can be downloaded via the [nextstrain-cli tool](https://docs.nextstrain.org/projects/cli/en/stable/). We do not provide this file due to its large size.
    *   Install the [Nextstrain CLI](https://docs.nextstrain.org/en/latest/install.html).
    *   In your terminal window, change the current directory to `6_LASSO_Predictive_Models/3_VOC_data/`.
    *   Download the Nextstrain-curated sequence metadata and rename the file to `gisaid_metadata.tsv.gz`:
```
nextstrain remote download s3://nextstrain-ncov-private/metadata.tsv.gz
mv metadata.tsv.gz gisaid_metadata.tsv.gz
```        

*   `multinomial_KC_voc_prediction_for_paper.R`
    *   Import `gisaid_metadata.tsv.gz` and filter metadata to sequences collected in King County, Washington during 2020 - 2022 with good QC scores. 
        *   Due to its large size, it will take several minutes to import `gisaid_metadata.tsv.gz`.
        *   Save King County metadata file as `king_county_SC2_clades_2020_2022.csv`.
    *   Use a multinomial logistic regression (MLR) model to predict the daily frequencies of SARS-CoV-2 clades circulating in King County.
    *   Figure S22: Predicted SARS-CoV-2 clade frequencies in King County, 2020 - 2022.
    *   Runtime: < 5 seconds (when starting with King County sequence metadata file `king_county_SC2_clades_2020_2022.csv`)

### Forecasting models
`4_Forecasting_Models/`
*   Use moving window LASSO regression models to produce 7-day ahead forecasts of rhinovirus, adenovirus, and SARS-CoV-2 Rt.
    *   Model covariates include activity of the target virus during the previous 2 weeks (14 autoregressive terms), cell phone mobility trends, the co-circulation of other viruses, and meteorological data (temperature, precipitation, humidity). The SARS-CoV-2 model also includes masking rates and NPI stringency as covariates.
    *   An additional SARS-CoV-2 model spanning 2021-2022 includes covariates for COVID-19 vaccination coverage and variant circulation.
*   Examine which model covariates were retained in moving model training windows.
*   Compare prediction accuracy across models with different combinations of covariates.
*   Runtime for each model script: < 4 minutes