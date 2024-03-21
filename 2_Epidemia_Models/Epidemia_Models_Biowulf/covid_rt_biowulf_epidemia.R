setwd("//home//perofskyamc//SFS_Rt_Epidemia")
library(epidemia)
library(EpiNow2)
library(ggplot2)
library(dplyr)
library(readr)
library(data.table)
library(tidyr)
library(EpiEstim)
library(lubridate)
library(zoo)
library(lognorm)

seattle_SG_and_covid_cases <- read_csv("wa_doh_king_county_cases.csv") #import king county daily covid-19 cases
load("symptom_onset/covid_reporting_delay.Rdata") # covid_reporting_delay
source("get_quantiles.R")
source("utils.R")

## add 7-day moving average to reduce noise but not oversmooth (covid cases already have 7-day moving average applied by WA DOH)
filter_fun2 <- function(x) {
  df <- x %>%
    mutate(confirm = round(zoo::rollmean(x = confirm, k = 7, fill = F, align = "center"))) %>%
    mutate(confirm = zoo::na.approx(confirm, maxgap = 3, na.rm = F)) %>%
    mutate(confirm = replace_na(confirm, 0)) %>%
    dplyr::select(date, confirm)%>%
    mutate(cum_inc = cumsum(confirm))
}
####################################
## COVID-19
####################################

COVID_2020_SKC <- seattle_SG_and_covid_cases %>% dplyr::select(Specimen.Collection.Date, Total.Cases.7.Day.Count)
names(COVID_2020_SKC) <- c("date", "confirm")
# head(COVID_2020_SKC)
# COVID_2020_SKC %>% filter(confirm>0) #first case 2/21

COVID_2020_SKC_smooth <- filter_fun2(COVID_2020_SKC %>% filter(date >= as.Date("2020-02-19")))
COVID_2020_SKC_smooth <- COVID_2020_SKC_smooth %>% filter(date < as.Date("2022-09-01"))
# head(COVID_2020_SKC_smooth,n=20)

COVID_2020_SKC_smooth <- COVID_2020_SKC_smooth %>% filter(date>as.Date("2020-02-24")) # smoothed data: cum cases reach 10 on 2/24
# head(COVID_2020_SKC_smooth,n=20)
# ggplot(COVID_2020_SKC_smooth)+
#   geom_line(aes(x=date,y=confirm))
# range(COVID_2020_SKC_smooth$date)

date <- as.Date("2020-02-24") + seq(0, along.with = c(NA, COVID_2020_SKC_smooth$confirm))

data <- data.frame(
  city = "Seattle", cases = c(NA, COVID_2020_SKC_smooth$confirm),
  date = date
)

data$day <- lubridate::wday(data$date, label = T) # day of the week effect
# head(data)
####################################
# transmission
####################################

## Rt in a given group at the given date
## simple random walk model
rt <- epirt(formula = R(city, date) ~ rw(time = date, gr = city, prior_scale = 0.01), link = "log",
            prior_intercept = normal(log(3), 0.2))
####################################
# observation model
####################################
# incubation period
mu1 <- EpiNow2::convert_to_logmean(6.3, 3.6)
sigma1 <- EpiNow2::convert_to_logsd(6.3, 3.6)

# symptom onset to report
mu2 <- covid_reporting_delay$mean ## already converted to lognormal
sigma2 <- covid_reporting_delay$sd

## infection to symptom onset + symptom onset to care seeking
combined_dist <- estimateSumLognormal(c(mu1, mu2), c(sigma1, sigma2))

## code borrowed from R0 package
tmax <- 15
t.scale <- c(0, 0.5 + c(0:tmax))
meanlog <- combined_dist[[1]]
sdlog <- combined_dist[[2]]
inf_to_case <- diff(plnorm(t.scale, meanlog = meanlog, sdlog = sdlog))
inf_to_case_adj <- inf_to_case / sum(inf_to_case) # make sure sums to 1

# observations are a function of latent infections in the population
obs <- epiobs(
  formula = cases ~ 0 + factor(day, ordered = FALSE),
  family = "neg_binom",
  link = "logit",
  prior_intercept = normal(0, 1.5),
  prior = normal(0, 0.5),
  prior_aux = normal(10, 5),
  i2o = inf_to_case_adj
)

####################################
# model latent infections
####################################
# Covid_generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")

covid_si <- R0::generation.time("gamma", c(5.2, 1.2), truncate = 15)

# basic model
inf <- epiinf(gen = covid_si$GT)

# extended model
# treat infections as latent parameters
inf_ext <- epiinf(
  gen = covid_si$GT, latent = TRUE,
  prior_aux = normal(10, 2)
) # coefficient of dispersion

####################################
# fitting the model
####################################

args <- list(
  rt = rt, obs = obs, inf = inf, data = data, iter = 40000, warmup = 20000,
  control = list(max_treedepth = 15,adapt_delta = 0.9),
  seed = 12345
)
args_ext <- args
args_ext$inf <- inf_ext

## basic model
# system.time(fm1 <- do.call(epim, args))

# options(mc.cores = parallel::detectCores())
## extended model this takes 9.5 hours to run on locally
fm2 <- do.call(epim, args_ext)

write_rds(fm2, file = "//data//perofskyamc//COVID_2020_2022_biowulf_epidemia_model.rds") #all output

####################################
## extract Rt 
####################################
all_covid_rt = extract_Rt(fm2)
all_covid_rt = all_covid_rt %>% group_by(group,tag,level) %>% tidyr::complete(date=full_seq(date,1))%>% ungroup()

write_rds(all_covid_rt, file = "//data//perofskyamc//COVID_2020_2022_biowulf_epidemia_model_rt_only.rds") #rt only
