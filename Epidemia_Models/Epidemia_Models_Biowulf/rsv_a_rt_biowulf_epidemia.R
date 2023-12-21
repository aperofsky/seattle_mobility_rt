setwd('//home//perofskyamc//SFS_Rt_Epidemia')
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

extract_Rt <-function(x){
  post_rt = posterior_rt(x)
  post_rt_quantiles = get_quantiles(post_rt,levels = c(90))
  
  post_rt_median <- data.frame(
    date = post_rt$time,
    median = apply(post_rt$draws, 2, function(x) quantile(x, 0.5)),
    group = post_rt$group)
  
  df = left_join(post_rt_quantiles,post_rt_median,by=c("date","group"))
  return(df)
}

# 3-week moving average to reduce noise but not oversmooth
filter_fun <- function(x){df = x %>%
  filter(organism=="RSV A")%>%
  mutate(confirm = all_sites_per_pos_ili_scaled_sum)%>%
  mutate(confirm = round(zoo::rollmean(x=confirm,k=21,fill=F,align = "center")*1000))%>%
  mutate(confirm = zoo::na.approx(confirm,maxgap = 3, na.rm=F))%>%
  mutate(confirm = tidyr::replace_na(confirm,0))%>%
  dplyr::select(date,confirm)}

####################################
## RSV A 2019
####################################
RSV_A_2019= filter_fun(pathogen_df %>% filter(date<as.Date("2019-06-30") & date > as.Date("2018-11-14")))

date <- as.Date("2018-11-14") + seq(0, along.with = c(NA, RSV_A_2019$confirm))

data <- data.frame(city = "Seattle", cases = c(NA, RSV_A_2019$confirm),
                   date = date)

data$day <- lubridate::wday(data$date, label=T) #day of the week effect
####################################
#transmission
####################################

## Rt in a given group at the given date
## simple random walk model
rt <- epirt(formula = R(city, date) ~ rw(time=date,gr=city,prior_scale = 0.01),link='log',
            prior_intercept = normal(log(2), 0.2)) 

####################################
#observation model
####################################

## RSV incubation period
mu1 = EpiNow2::convert_to_logmean(4.5,0.894)
sigma1 = EpiNow2::convert_to_logsd(4.5,0.894)

# symptom onset to report
mu2 <- cold_virus_reporting_delay$mean ## mean and sd are already on lognormal scale
sigma2 <- cold_virus_reporting_delay$sd

## infection to symptom onset + symptom onset to care seeking
combined_dist = estimateSumLognormal(c(mu1,mu2), c(sigma1,sigma2))

## code borrowed from R0 package
tmax = 15
t.scale = c(0,0.5+c(0:tmax))
meanlog = combined_dist[[1]]
sdlog = combined_dist[[2]]
inf_to_case = diff(plnorm(t.scale,meanlog=meanlog,sdlog=sdlog))
inf_to_case_adj = inf_to_case/sum(inf_to_case) #make sure sums to 1


#observations are a function of latent infections in the population
obs <-  epiobs(formula = cases ~ 0 + factor(day, ordered = FALSE), 
               family = "neg_binom",
               link = "logit",
               prior_intercept = normal(0,1.5), 
               prior = normal(0,0.5),
               prior_aux = normal(10,5),
               i2o = inf_to_case_adj)

####################################
# model latent infections
####################################
rsv_a_si = R0::generation.time(type="gamma",c(7.5,2.5),truncate = 26)

#basic model
inf <- epiinf(gen = rsv_a_si$GT)

#extended model
#treat infections as latent parameters
inf_ext <-  epiinf(gen = rsv_a_si$GT, latent = TRUE, 
                   prior_aux = normal(10,2)) #coefficient of dispersion

####################################
#fitting the model
####################################

args <- list(rt = rt, obs = obs, inf = inf, data = data, iter = 30000, warmup = 15000,
             control = list(max_treedepth = 15),
             seed = 12345)
args_ext <- args; args_ext$inf <- inf_ext

##basic model
# system.time(fm1 <- do.call(epim, args))

options(mc.cores = parallel::detectCores())
##extended model
system.time(fm2 <- do.call(epim, args_ext))

write_rds(fm2,file="//data//perofskyamc//RSV_A_2019_biowulf_epidemia_model.rds")

RSV_A_2019 = extract_Rt(fm2)

#####################################
## RSV A 2020
####################################
RSV_A_2020 = filter_fun(pathogen_df %>% filter(date>as.Date("2019-09-25") & date<as.Date("2020-04-30")))

date <- as.Date("2019-09-25") + seq(0, along.with = c(NA, RSV_A_2020$confirm))

data <- data.frame(city = "Seattle", cases = c(NA, RSV_A_2020$confirm),
                   date = date)

data$day <- lubridate::wday(data$date, label=T) #day of the week effect

####################################
#fitting the model
####################################

args <- list(rt = rt, obs = obs, inf = inf, data = data, iter = 30000, warmup = 15000,
             control = list(max_treedepth = 15),
             seed = 12345)
args_ext <- args; args_ext$inf <- inf_ext

##basic model
# system.time(fm1 <- do.call(epim, args))

options(mc.cores = parallel::detectCores())
##extended model
system.time(fm2_2020 <- do.call(epim, args_ext))

write_rds(fm2_2020,file="//data//perofskyamc//RSV_A_2020_biowulf_epidemia_model.rds")

RSV_A_2020 = extract_Rt(fm2_2020)

####################################
## combine all Rt output
####################################
all_rsv_a_rt = all_rsv_a_rt %>% group_by(group,tag,level) %>% tidyr::complete(date=full_seq(date,1))%>% ungroup()
all_rsv_a_rt$organism = "RSV A"
# range(all_rsv_a_rt$date)
write_rds(all_rsv_a_rt, file = "//data//perofskyamc//RSV_A_epidemia_model_rt_only.rds") #rt only
