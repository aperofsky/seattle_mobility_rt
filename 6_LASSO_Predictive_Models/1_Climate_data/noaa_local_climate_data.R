library(readr)
library(dplyr)
library(zoo)
library(tidyr)
# National Oceanic and Atmospheric Administration Local Climatological Data; downloaded for SeaTac on 2023-01-30

# from: 10.1001/jamanetworkopen.2020.16099
# wet-bulb temperature: captures the complex thermodynamic relationship of temperature and humidity, has been shown to predict human health events with more precision than temperature and humidity separately, and avoids the associated problem of collinearity.

# Special Indicator Appendix:
# s = suspect value (appears with value)
# T = trace precipitation amount or snow depth (an amount too small to measure, usually < 0.005 inches water equivalent) (appears instead of numeric value)
# M = missing value (appears instead of value)
# Blank = value is unreported (appears instead of value)

daily_local_df <- read_csv("6_LASSO_Predictive_Models/1_Climate_data/3216704.csv",
  col_select = c(
    "DATE",
    "DailyAverageWetBulbTemperature",
    "DailyAverageRelativeHumidity",
    "DailyPrecipitation")
) %>%
  filter(!is.na(DailyAverageWetBulbTemperature))
names(daily_local_df)
head(daily_local_df)
range(daily_local_df$DATE)

daily_local_df$DATE <- as.Date(daily_local_df$DATE)

unique(daily_local_df$DailyAverageRelativeHumidity)
unique(daily_local_df$DailyPrecipitation)
unique(daily_local_df$DailyAverageWetBulbTemperature)


daily_local_df <- daily_local_df %>%
  mutate(DailyPrecipitation = ifelse(DailyPrecipitation == "T", 0.00, DailyPrecipitation)) %>%
  mutate(DailyPrecipitation = as.numeric(DailyPrecipitation))

daily_local_df <- daily_local_df %>%
  rename(date=DATE)%>%
  complete(date = seq.Date(as.Date("2018-01-01"), as.Date("2022-09-30"), by = "day")) %>%
  mutate(
    DailyPrecipitation = zoo::na.approx(DailyPrecipitation, maxgap = 4),
    DailyAverageRelativeHumidity = zoo::na.approx(DailyAverageRelativeHumidity, maxgap = 4),
    DailyAverageWetBulbTemperature = zoo::na.approx(DailyAverageWetBulbTemperature, maxgap = 4)
  )
head(daily_local_df)
range(daily_local_df$date)#"2018-01-01" "2022-09-30"
write_csv(daily_local_df, file = "6_LASSO_Predictive_Models/1_Climate_data/seattle_daily_weather_variables_2023_01_30.csv")
