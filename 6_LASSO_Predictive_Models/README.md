## Short-term forecasts of rhinovirus, adenovirus, and SARS-CoV-2 Rt

### Model inputs
`1_Climate_data/`
*   `noaa_local_climate_data.R`
    *   Import, clean, and interpolate missing values for meteorological data collected in Seattle, Washington. Daily records of precipitation, average wet bulb temperature, and average relative humidity are publicly accessible through the National Centers for Environmental Information’s [U.S. Local Climatological Database](https://www.ncei.noaa.gov/products/land-based-station/local-climatological-data).

`2_Vaccination_coverage_data/`
*   Daily COVID-19 vaccination coverage in King County, Washington. Data are publicly accessible through the [Public Health – Seattle & King County COVID-19 Vaccination dashboard](https://kingcounty.gov/en/dept/dph/health-safety/disease-illness/covid-19/data/vaccination).

`3_VOC_data/`
*   Nextstrain-curated SARS-CoV-2 sequence metadata from GISAID can be downloaded via the [nextstrain-cli tool](https://docs.nextstrain.org/projects/cli/en/stable/). We do not provide this file due to its large size.
    *   Install [nextstrain-cli](https://docs.nextstrain.org/projects/cli/en/stable/).
    *   In your terminal window, change the current directory to `6_LASSO_Predictive_Models/3_VOC_data/`.
    *   Download the Nextstrain-curated sequence metadata and rename the file to `gisaid_metadata.tsv.gz`:
```
nextstrain remote download s3://nextstrain-ncov-private/metadata.tsv.gz
mv metadata.tsv.gz gisaid_metadata.tsv.gz
```        

*   `multinomial_KC_voc_prediction_for_paper.R`
    *   Import `gisaid_metadata.tsv.gz` and filter metadata to sequences collected in King County, Washington with good QC scores. Due to its large size, it may take a few minutes to import`gisaid_metadata.tsv.gz`.
    *   Use a multinomial logistic regression (MLR) model to predict the daily frequencies of SARS-CoV-2 clades circulating in King County.
    *   Figure S22: Predicted SARS-CoV-2 clade frequencies in King County over time.         

### Forecasting models
`4_Forecasting_Models/`
*   Use moving window LASSO regression models to produce 7-day ahead forecasts of rhinovirus, adenovirus, and SARS-CoV-2 Rt.
*   Model covariates include activity of the target virus during the previous 2 weeks (14 autoregressive terms) and lagged cell phone mobility trends, the co-circulation of other viruses, and meteorological data (temperature, precipitation, humidity). The SARS-CoV-2 model also includes masking rates and NPI stringency as covariates.
*   An additional SARS-CoV-2 model spanning 2021-2022 includes covariates for COVID-19 vaccination coverage and variant circulation.
*   Examine which model covariates were retained in moving model training windows.
*   Compare prediction accuracy across models with different combinations of covariates.