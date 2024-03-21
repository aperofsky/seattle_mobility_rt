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
filter_fun <- function(x){df = x %>%
  filter(organism=="IVA/H3N2")%>%
  mutate(confirm = all_sites_per_pos_ili_scaled_sum)%>%
  mutate(confirm = round(zoo::rollmean(x=confirm,k=21,fill=F,align = "center")*1000))%>%
  mutate(confirm = zoo::na.approx(confirm,maxgap = 3, na.rm=F))%>%
  mutate(confirm = tidyr::replace_na(confirm,0))%>%
  dplyr::select(date,confirm)%>%
  mutate(cum_inc = cumsum(confirm))
}

####################################
## H3N2 2019
####################################
H3N2_2019 <- filter_fun(pathogen_df %>% filter(date < as.Date("2019-05-30") & date>as.Date("2018-11-13")))%>% filter(date>as.Date("2018-12-16"))
# ggplot(H3N2_2019)+
#   geom_line(aes(x=date,y=confirm))
# range(H3N2_2019$date)

date <- as.Date("2018-12-15") + seq(0, along.with = c(NA, H3N2_2019$confirm))

data <- data.frame(
  city = "Seattle", cases = c(NA, H3N2_2019$confirm),
  date = date
)

data$day <- lubridate::wday(data$date, label = T) # day of the week effect
# head(data)
####################################
# transmission
####################################

## Rt in a given group at the given date
## simple random walk model
rt <- epirt(formula = R(city, date) ~ 1 + rw(time = date, gr = city, prior_scale = 0.01), link = "log",
            prior_intercept = normal(log(1.2), 0.2))

####################################
# observation model
####################################

## influenza incubation period
mu1 <- EpiNow2::convert_to_logmean(1.9, 1.22)
sigma1 <- EpiNow2::convert_to_logsd(1.9, 1.22)

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

H3N2_si <- R0::generation.time("weibull", c(3.6, 1.6))

# basic model
inf <- epiinf(gen = H3N2_si$GT)

# extended model
# treat infections as latent parameters
inf_ext <- epiinf(
  gen = H3N2_si$GT, latent = TRUE,
  prior_aux = normal(10, 2)
) # coefficient of dispersion

####################################
# fitting the model
####################################

args <- list(
  rt = rt, obs = obs, inf = inf, data = data, iter = 30000, warmup = 15000,
  control = list(max_treedepth = 15,adapt_delta = 0.9),
  seed = 12345
)
args_ext <- args
args_ext$inf <- inf_ext

## basic model
# system.time(fm1 <- do.call(epim, args))

options(mc.cores = parallel::detectCores())
## extended model
system.time(fm2 <- do.call(epim, args_ext))

write_rds(fm2, file = "//data//perofskyamc//H3N2_2019_biowulf_epidemia_model.rds")

####################################
## extract Rt 
####################################

all_H3N2_rt = extract_Rt(fm2)

all_H3N2_rt = all_H3N2_rt %>% group_by(group,tag,level) %>% tidyr::complete(date=full_seq(date,1))%>% ungroup()
all_H3N2_rt$organism = "IVA/H3N2"

write_rds(all_H3N2_rt, file = "//data//perofskyamc//H3N2_biowulf_epidemia_model_rt_only.rds") #rt only

#####################################
## H3N2 2022
####################################
## H3N2 circulated in winter 2021-2022 but not enough cases to estimate Rt = not included in manuscript
# H3N2_2022 <- filter_fun(pathogen_df %>% filter(date > as.Date("2021-12-30") & date < as.Date("2022-08-30"))) %>% filter(date < as.Date("2022-07-01"))
# 
# # ggplot(H3N2_2022)+
# #   geom_line(aes(x=date,y=confirm))
# # not enough cases
# 
# date <- as.Date("2021-12-30") + seq(0, along.with = c(NA, H3N2_2022$confirm))
# 
# data <- data.frame(
#   city = "Seattle", cases = c(NA, H3N2_2022$confirm),
#   date = date
# )
# 
# data$day <- lubridate::wday(data$date, label = T) # day of the week effect
# # ####################################
# # #fitting the model
# # ####################################
# #
# args <- list(
#   rt = rt, obs = obs, inf = inf, data = data, iter = 40000, warmup = 20000,
#   control = list(max_treedepth = 15,adapt_delta = 0.9),
#   seed = 12345
# )
# args_ext <- args
# args_ext$inf <- inf_ext
# 
# ## basic model
# # system.time(fm1 <- do.call(epim, args))
# 
# ## extended model
# system.time(fm2_2022 <- do.call(epim, args_ext))
# 
# write_rds(fm2_2022, file = "//data//perofskyamc//H3N2_2022_biowulf_epidemia_model.rds")
