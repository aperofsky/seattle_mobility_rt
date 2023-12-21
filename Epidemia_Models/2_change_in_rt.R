###################################################
## Compare Rt before and after 2 major events:
## 1. Feb 2019 snowstorm
## 2. Feb 28 2020 state of emergency declaration
## Method: t-test for the ratio of two means and Fiellers confidence interval for the ratio of two means (mratios package)
###################################################

library(dplyr)
library(tidyr)
library(readr)
library(mratios)
library(magicfor)

all_rt <- read_csv("Epidemia_Models/rt_all_pathogens_raw.csv")

all_rt %>%
  pull(organism) %>%
  unique()

###################################################
## two weeks before and during major snowstorm in Feb 2019
###################################################
## 2 weeks before snowstorm
rt_pre_snow <- all_rt %>%
  filter(level == 90) %>%
  filter(date < as.Date("2019-02-03") & date > (as.Date("2019-02-03") - 15)) %>%
  dplyr::select(date, organism, median) %>%
  rename(pre_snow_rt = median)

unique(rt_pre_snow$date)
length(unique(rt_pre_snow$date))

## 2 weeks during snowstorm
rt_post_snow <- all_rt %>%
  filter(level == 90) %>%
  filter(date > as.Date("2019-02-02") & date < (as.Date("2019-02-02") + 15)) %>%
  dplyr::select(date, organism, median) %>%
  rename(post_snow_rt = median)
unique(rt_post_snow$date)
length(unique(rt_post_snow$date))

path_list <- unique(rt_post_snow$organism)[!(unique(rt_post_snow$organism) %in% c("IVB", "HPIV-1 + HPIV-2"))]
path_list

## ttestratio() function performs t-test for the ratio of two means and Fiellers confidence interval for the ratio of two means 
magic_for(silent = TRUE, temporary = T)
for (i in path_list) {
  x <- rt_pre_snow %>%
    filter(organism == i) %>%
    pull(pre_snow_rt)
  y <- rt_post_snow %>%
    filter(organism == i) %>%
    pull(post_snow_rt)

  m <- ttestratio(x = x, y = y)
  mean <- 1 - m$estimate[3]
  p_value <- round(m$p.value, 5)
  ci <- 1 - m$conf.int #Fieller's interval for the ratio of two means
  ci <- ci[c(2, 1)]
  put(mean, ci[1], ci[2], p_value)
}
df <- magic_result_as_dataframe()
df <- df %>%
  arrange(mean) %>%
  mutate_at(c("mean", "ci[1]", "ci[2]"), ~ 100 * .)
df$p_value <- format.pval(df$p_value, digits = 2)
df

# i                     mean        ci[1]      ci[2] p_value
# 1                     RSV A -39.838111 -50.35065084 -29.577599  <2e-16
# 2                Adenovirus -38.591973 -47.12772195 -30.553405  <2e-16
# 3                     RSV B -36.485991 -45.20051568 -28.009628  <2e-16
# 4                  IVA/H1N1 -19.381274 -25.27726931 -13.506112   1e-05
# 5  Seasonal CoV HKU1 + NL63 -17.933153 -22.21384720 -13.842456  <2e-16
# 6                  IVA/H3N2 -12.258980 -16.69147134  -7.884669   2e-05
# 7  Seasonal CoV 229E + OC43  -9.526834 -14.51078674  -4.566976   0.001
# 8                Rhinovirus   1.525658  -0.03163833   3.077609   0.054
# 9           HPIV-3 + HPIV-4  13.886475   9.12073125  18.232058   4e-05
# 10                     hMPV  19.492694  15.98544646  22.774193  <2e-16

###################################################
## two weeks before and after WA state of emergency declaration
###################################################

### SARS-CoV-2 Rt leading up to lifting of SAH in early June 2020
all_rt %>%
  filter(organism == "SARS-CoV-2") %>%
  filter(date >= as.Date("2020-04-01") & date< as.Date("2020-06-15")) %>%
  dplyr::select(date,median) %>%
  tail(n=20)

all_rt %>%
  filter(organism == "SARS-CoV-2") %>%
  filter(date >= as.Date("2020-04-01") & date< as.Date("2020-06-01")) %>%
  dplyr::select(date,median) %>%
  pull(median) %>% range()

all_rt %>%
  filter(organism == "SARS-CoV-2") %>%
  filter(date >= as.Date("2020-04-01") & date< as.Date("2020-05-20")) %>%
  dplyr::select(date,median) %>%
  pull(median) %>% range() #0.8104679 0.9197895

all_rt %>%
  filter(organism == "SARS-CoV-2") %>%
  filter(date >= as.Date("2020-04-01") & date< as.Date("2020-06-02")) %>%
  dplyr::select(date,median) %>%
  pull(median) %>% range() #0.8104679 1.0012415; surpasses 1 on June 1

## pathogens circulating during the 2 weeks prior to State of Emergency
all_rt %>%
  filter(level == 90) %>%
  filter(date < as.Date("2020-02-29") & date > (as.Date("2020-02-29") - 15)) %>%
  dplyr::select(date, organism, median) %>%
  pull(organism) %>% unique()

## 2 weeks before State of Emergency
rt_pre_ED <- all_rt %>%
  filter(level == 90) %>%
  filter(date < as.Date("2020-02-29") & date > (as.Date("2020-02-29") - 15)) %>%
  dplyr::select(date, organism, median) %>%
  rename(pre_ED_rt = median)
unique(rt_pre_ED$date)
length(unique(rt_pre_ED$date))

## 2 weeks after State of Emergency
rt_post_ED <- all_rt %>%
  filter(level == 90) %>%
  filter(date > as.Date("2020-02-28") & date < (as.Date("2020-02-28") + 15)) %>%
  dplyr::select(date, organism, median) %>%
  rename(post_ED_rt = median)
unique(rt_post_ED$date)
length(unique(rt_post_ED$date))

path_list <- unique(rt_post_ED$organism)[!(unique(rt_post_ED$organism) %in% c("Influenza A/H3N2"))] %>% as.character()
path_list
magic_for(silent = TRUE, temporary = T)
for (i in path_list) {
  x <- rt_pre_ED %>%
    filter(organism == i) %>%
    pull(pre_ED_rt)
  y <- rt_post_ED %>%
    filter(organism == i) %>%
    pull(post_ED_rt)

  m <- ttestratio(x = x, y = y)
  mean <- 1 - m$estimate[3]
  p_value <- round(m$p.value, 5)
  ci <- 1 - m$conf.int #Fieller's interval for the ratio of two means
  ci <- ci[c(2, 1)]
  put(mean, ci[1], ci[2], p_value)
}
df <- magic_result_as_dataframe()
df <- df %>%
  arrange(mean) %>%
  mutate_at(c("mean", "ci[1]", "ci[2]"), ~ 100 * .)
df$p_value <- format.pval(df$p_value, digits = 2)
df
#                           i       mean      ci[1]      ci[2] p_value
# 1                SARS-CoV-2 -49.512829 -73.253718 -31.427108 < 2e-16
# 2                  IAV/H1N1 -42.444074 -59.205418 -28.771587 < 2e-16
# 3                Adenovirus -32.527591 -45.367890 -21.600040 < 2e-16
# 4                        RV -29.174592 -36.584332 -22.375984 < 2e-16
# 5  Seasonal CoV HKU1 + NL63 -25.170527 -31.165713 -19.591476 < 2e-16
# 6           HPIV-3 + HPIV-4 -21.531402 -27.548681 -15.796916 < 2e-16
# 7                      hMPV -20.487794 -25.616631 -15.472855 < 2e-16
# 8  Seasonal CoV 229E + OC43 -19.947531 -24.504246 -15.574026 < 2e-16
# 9                       IBV -14.486521 -19.199341 -10.103157 < 2e-16
# 10                    RSV A  -8.936362 -13.841364  -4.425130 0.00057
# 11                    RSV B  -3.750857  -5.477746  -2.079393 0.00026
