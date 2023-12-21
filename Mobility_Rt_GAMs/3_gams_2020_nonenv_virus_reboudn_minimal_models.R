########################################################################################
## GAMS: rhinovirus and adenovirus rebound, June 2020 - Dec 2020
########################################################################################

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(lubridate)
library(mgcv)
library(MuMIn)
library(gratia)
library(patchwork)

########################################################################################
## import data
########################################################################################

## mobility data
combined_mob <- read_rds("Epidemia_Models/mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

## daily rt
rt_df <- read_csv("Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv")

rt_df_wide <- rt_df %>%
  dplyr::select(date, median, organism) %>%
  distinct() %>%
  pivot_wider(names_from = "organism", values_from = "median") %>%
  dplyr::filter(date > as.Date("2018-11-20"))

names(rt_df_wide)[2:14] <- c(
  "rsv_a", "rsv_b",
  "h1n1", "h3n2", "ivb", "rhino",
  "hpiv_3_4", "hpiv_1_2", "hmpv", "adeno", "scov_229E_OC43",
  "scov_HKU1_NL63", "covid"
)
combined = left_join(rt_df_wide,
                     combined_mob_avg,by="date")%>% filter(date<as.Date("2022-06-01"))

combined[is.na(combined)]<- 0 
combined$fb_leaving_home_custom = 100 - combined$fb_stay_put_custom
combined$not_masking = 100 - combined$mask_wearing_final

cols2= c(
  # "within_neighborhood_movement",
         "within_city_movement",
         # "within_state_movement",
         "out_of_state_movement",
         "fb_leaving_home_custom",
         "full_service_restaurants",
         "religious_orgs",
         # "child_day_care",
         # "not_masking",
         # "oxford_stringency_index",
         "elementary_and_secondary_schools"
         # "colleges"
         )

com_df = combined %>% 
  dplyr::select(date,all_of(cols2),rsv_a:covid)%>%
  filter(date>=as.Date("2019-01-08") & date<=as.Date("2022-05-01"))%>%
  complete(date = seq.Date(as.Date("2019-01-08"),as.Date("2022-05-01"),by="day"))

###############################################################
######### June to December 2020 ##################
###############################################################
## first 6 months of rebound
com_2020 = com_df %>%filter(date>=as.Date("2020-06-01")&date<as.Date("2020-12-01"))

com_2020_df_lm = 
  com_2020 %>%
  mutate(time = 1 + (date - min(date)) / ddays()) %>% # Adding numeric times
  mutate(week = isoweek(date))%>%
  mutate(stay_at_home = if_else(date>=as.Date("2020-02-28") & date<=as.Date("2020-06-06"),1,0))

week_df = data.frame(week = unique(com_2020_df_lm$week),
                     week_index = seq(1:length(unique(com_2020_df_lm$week))))
week_df

com_2020_df_lm=left_join(com_2020_df_lm,week_df,by="week")
names(com_2020_df_lm)

com_2020_df_lm_scaled = 
  com_2020_df_lm %>%
  mutate_at(vars(within_city_movement:elementary_and_secondary_schools),~scale(.x)%>% as.vector)

############################
### rhinovirus
############################

rhino_df_lim = com_2020_df_lm_scaled %>% 
  dplyr::select(date,time,week_index,all_of(cols2),rhino) %>% filter(rhino>0)
head(rhino_df_lim)

fm <- paste("s(", names(rhino_df_lim[!(names(rhino_df_lim)%in%c("date","time","rhino"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste('rhino ~', fm))
fm

b1 <- mgcv::gam(formula=fm,data=rhino_df_lim,method = 'ML', gamma = 2, 
                select=TRUE,na.action="na.fail",family=Gamma(link="log"))
gam.check(b1)
dd <- dredge(b1,m.lim = c(1,3),fixed="s(week_index)")
subset(dd, delta < 4)
best_mod = get.models(dd, subset = 1)[[1]]
summary(best_mod)
rhino_b1 <- update(best_mod,method="REML")
summary(rhino_b1)

p = draw(rhino_b1, residuals = TRUE,scales="free",nrow=1)
p = p & theme_bw()
p

p1 <- draw(rhino_b1, residuals=T,select = "s(fb_leaving_home_custom)") & theme_bw()& ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p1
p2 <- draw(rhino_b1, residuals=T,select = "s(out_of_state_movement)") & theme_bw()& ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(rhino_b1, residuals=T,select = "s(week_index)") & theme_bw()& ggtitle("s(week)") & xlab("week")
p3

rhino_p_all = p1 + p2 + p3 + plot_layout(ncol = 3)+ 
  plot_annotation(title = 'Rhinovirus',theme = theme(plot.title = element_text(size = 18,face="bold")))
rhino_p_all

############################
### adenovirus
############################
names(com_2020_df_lm)
adeno_df_lim = com_2020_df_lm_scaled %>% 
  dplyr::select(date,time,week_index,all_of(cols2),adeno) %>% filter(adeno>0)

fm <- paste("s(", names(adeno_df_lim[!(names(adeno_df_lim)%in%c("date","time","adeno"))]), ")", sep = "", collapse = " + ")
fm <- as.formula(paste('adeno ~', fm))
fm
b1 <- mgcv::gam(formula=fm,data=adeno_df_lim,method = 'ML', gamma = 2, select=TRUE,
                family=Gamma(link="log"),na.action="na.fail")
gam.check(b1)
dd <- dredge(b1,m.lim = c(1,3),fixed="s(week_index)")
subset(dd, delta < 4)
best_mod = get.models(dd, subset = 1)[[1]]
summary(best_mod)
adeno_b1 <- update(best_mod,method="REML")
summary(adeno_b1)

p = draw(adeno_b1, residuals = TRUE,scales="free",nrow=1)
p = p & theme_bw()
p

p1 <- draw(adeno_b1, residuals=T,select = "s(fb_leaving_home_custom)") & theme_bw()& ggtitle("s(% devices leaving home)") & xlab("% devices leaving home")
p1
p2 <- draw(adeno_b1, residuals=T,select = "s(out_of_state_movement)") & theme_bw()& ggtitle("s(out-of-state inflow)") & xlab("out-of-state inflow")
p2
p3 <- draw(adeno_b1, residuals=T,select = "s(week_index)") & theme_bw()& ggtitle("s(week)") & xlab("week")

adeno_p_all = p1 + p2 + p3 + plot_layout(ncol = 3)+ 
  plot_annotation(title = 'Adenovirus',theme = theme(plot.title = element_text(size = 18,face="bold")))
adeno_p_all

########################################################################################
## combine all plots into one figure
########################################################################################
non_env_p_2020 = plot_grid(rhino_p_all,adeno_p_all,nrow=2)
save_plot(non_env_p_2020,base_height = 8,base_width = 16,filename="figures/fig_s13_partial_effects_non_env_viruses_june_to_dec_2020.png")

