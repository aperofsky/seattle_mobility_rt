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
library(R0)

pathogen_df <- readr::read_rds("extended_daily_ILI_per_pos_scaled_for_select_pathogens_2023_10_22.rds")
load("symptom_onset/cold_virus_reporting_delay.Rdata")
source("get_quantiles.R")
source("utils.R")

# 3-week moving average to reduce noise but not oversmooth
filter_fun <- function(x) {
  df <- x %>%
    filter(organism == "Rhinovirus") %>%
    mutate(confirm = all_sites_per_pos_ili_scaled_sum) %>%
    mutate(confirm = round(zoo::rollmean(x = confirm, k = 21, fill = F, align = "center") * 1000)) %>%
    mutate(confirm = zoo::na.approx(confirm, maxgap = 3, na.rm = F)) %>%
    mutate(confirm = tidyr::replace_na(confirm, 0)) %>%
    dplyr::select(date, confirm) %>%
    mutate(cum_inc = cumsum(confirm))
}

####################################
## Rhinovirus 2018 - 2022
####################################
Rhino_2018_2022 <- filter_fun(pathogen_df %>%
  filter(date > as.Date("2018-11-20") & date < as.Date("2022-08-30"))) %>%
  filter(date > as.Date("2018-11-30") & date < as.Date("2022-07-01"))
# range(Rhino_2018_2022$date)
# ggplot(Rhino_2018_2022) +
#   geom_line(aes(x = date, y = confirm))

date <- as.Date("2018-11-30") + seq(0, along.with = c(NA, Rhino_2018_2022$confirm))

data <- data.frame(
  city = "Seattle", cases = c(NA, Rhino_2018_2022$confirm),
  date = date
)

data$day <- lubridate::wday(data$date, label = T) # day of the week effect
# head(data)
####################################
# transmission
####################################

## Rt in a given group at the given date
## simple random walk model
rt <- epirt(
  formula = R(city, date) ~ rw(time = date, gr = city, prior_scale = 0.01), link = "log",
  prior_intercept = normal(log(2), 0.2)
)

####################################
# observation model
####################################
# incubation period
mu1 <- EpiNow2::convert_to_logmean(2.36, 1.1)
sigma1 <- EpiNow2::convert_to_logsd(2.36, 1.1)

# symptom onset to report
mu2 <- cold_virus_reporting_delay$mean ## mean and sd are already on lognormal scale
sigma2 <- cold_virus_reporting_delay$sd

## infection to symptom onset + symptom onset to care seeking
combined_dist <- estimateSumLognormal(c(mu1, mu2), c(sigma1, sigma2))

## code borrowed from R0 package
tmax <- 15
t.scale <- c(0, 0.5 + c(0:tmax))
meanlog <- combined_dist[[1]]
sdlog <- combined_dist[[2]]
inf_to_case <- diff(plnorm(t.scale, meanlog = meanlog, sdlog = sdlog))
inf_to_case_adj <- inf_to_case / sum(inf_to_case) # make sure sums to 1


# observations are a function of latent infections in the population (ex: daily cases)
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

rhino_si <- R0::generation.time("gamma", c(4.4, 2.7), truncate = 21)
## truncate at 21 days because hRV can be shed for a few weeks

# basic model
inf <- epiinf(gen = rhino_si$GT)

# extended model
# treat infections as latent parameters
inf_ext <- epiinf(
  gen = rhino_si$GT, latent = TRUE,
  prior_aux = normal(10, 2)
) # coefficient of dispersion

####################################
# fitting the model
####################################

args <- list(
  rt = rt, obs = obs, inf = inf, data = data, iter = 30000, warmup = 15000,
  control = list(max_treedepth = 15),
  seed = 12345
)
args_ext <- args
args_ext$inf <- inf_ext

## basic model
# system.time(fm1 <- do.call(epim, args))

options(mc.cores = parallel::detectCores())
fm2 <- do.call(epim, args_ext) # this can take days to run on locally

write_rds(fm2, file = "//data//perofskyamc//Rhino_2018_2022_biowulf_epidemia_model.rds")

####################################
## extract Rt
####################################

rhino_rt_2018_2022 <- extract_Rt(fm2)

all_rhino_rt <- rhino_rt_2018_2022 %>%
  group_by(group, tag, level) %>%
  tidyr::complete(date = full_seq(date, 1)) %>%
  ungroup()
all_rhino_rt$organism <- "Rhinovirus"

write_rds(all_rhino_rt, file = "//data//perofskyamc//Rhino_2018_2022_biowulf_epidemia_model_rt_only.rds") # rt only
