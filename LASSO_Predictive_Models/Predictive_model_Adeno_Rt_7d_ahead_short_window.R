########################################################################################
## Predictive model of adenovirus Rt
########################################################################################

library(readr)
library(dplyr)
library(tidyverse)
library(glue)

library(forecast)
library(tibbletime)
library(caret)
library(lubridate)
library(zoo)

library(glmnet)
library(Metrics)
library(MLmetrics)

library(ggplot2)
library(cowplot)
library(ggExtra)
library(tidyquant)
library(ggpubr)
########################################################################################
## Import data
########################################################################################

combined_mob <- read_rds("Epidemia_Models/mobility_metrics_for_epidemia.rds") %>% as_tibble()
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)

combined_mob_avg <- combined_mob %>%
  dplyr::mutate(across(where(is.numeric), ~ zoo::rollmean(.x, k = 15, align = "center", fill = "extend")))

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
combined <- left_join(rt_df_wide,
                      combined_mob_avg,
                      by = "date"
) %>% filter(date < as.Date("2022-06-01"))


cols2= c("within_neighborhood_movement",
         "within_city_movement",
         "within_state_movement",
         "out_of_state_movement",
         "fb_leaving_home_custom",
         "full_service_restaurants",
         "groceries_and_pharmacies",
         "full_service_restaurants",
         "transit",
         "religious_orgs",
         "child_day_care",
         "elementary_and_secondary_schools",
         "colleges")

cols2
max_date = combined %>%
  filter(rhino>0)%>%
  slice_max(date) %>%
  pull(date)
max_date

min_date = combined %>%
  filter(rhino>0)%>%
  slice_min(date) %>%
  pull(date)
min_date

combined[is.na(combined)]<- 0 
combined$fb_leaving_home_custom = 100 - combined$fb_stay_put_custom
combined[!complete.cases(combined),] %>% data.table::data.table()

adeno_df = combined %>% 
  dplyr::select(date,all_of(cols2),adeno,covid)%>%
  filter(date>=as.Date("2019-01-08") & date<=as.Date("2022-04-30"))%>%
  complete(date = seq.Date(as.Date("2019-01-08"),as.Date("2022-04-30"),by="day"))

########################################################################################
## Moving window LASSO predictive model
########################################################################################

adeno_df_scaled=adeno_df

window_train <- 30 # number of days of the training window
w = c(rep(1, window_train-28), rep(2, 28)) #observation weights (last month has double the weight)
lag_AR <- 1 # minimum lag (in days) of the autoregressive term (sometimes surveillance data may have more than 1 day lag)
window_AR <- 14 # time window of the autoregressive terms (in days)
t0_train <- 15 # starting day of the training
dt_predict <- 7 # prediction ahead relative to the last data point of the training set
alpha <- 1 # weight between L1 and L2 penalty alpha = 1 is Lasso regression
nfolds <- 10 # nfolds crossvalidation to evaluate hyperparameter lambda for elastic net regression
s <- "lambda.min" # choose lambda either min cross validation error (lambda.min) or most regularized (lambda.1se)

var_target <- "adeno" # dependent variable
var_regressors_GT <- select(adeno_df_scaled,-date,-adeno,-covid) %>% names()
var_regressors_GT_COVID <- select(adeno_df_scaled,-date,-adeno) %>% names()
var_regressors_AR <- sprintf("AR%02d", lag_AR:(lag_AR + window_AR - 1)) # regressors, autoregressive
var_regressors_ARGO <- c(var_regressors_AR, var_regressors_GT_COVID) # regressors, ARGO
var_predicts <- sprintf("predict_%dd_ahead", 1:dt_predict)
var_data <- sprintf("data_%dd_ahead", 1:dt_predict)

min_date = min(adeno_df_scaled$date)

# make reference dates
reference_dates = data.frame(num_days = seq(from = -1500, to = 1500, by = 1))
reference_dates[["date"]] = reference_dates[["num_days"]] + ymd(min_date)
reference_dates[["epi_year"]] = epiyear(reference_dates[["date"]])
reference_dates[["epi_week"]] = epiweek(reference_dates[["date"]])

rdf = left_join(adeno_df_scaled,reference_dates)
tail(rdf[,c("num_days")])
day_max =max(rdf$num_days)

t0_predict_start <- t0_train + window_train
tmax_predict_start <- day_max - dt_predict + 1

# Set up lags
lags <- c(1:14)

# Name the columns that will contain the lags, with appropriate index
lag_names <- glue::glue('AR{str_pad(lags, nchar(max(lags)), pad = "0")}')

# Create list of lag functions, as eval-ed/parse-d labmda functions
lag_functions <-
  map(lags, ~ eval(parse(text=glue::glue("~ dplyr::lag(.x, {.x})")))) %>%
  set_names(lag_names)

adeno_agg = rdf %>% 
  dplyr::select(date,num_days,epi_year,epi_week,adeno)%>%
  mutate_at(vars(adeno), .funs=lag_functions)

names(adeno_agg)
head(adeno_agg)
tail(adeno_agg)

adeno_agg <- gather(adeno_agg, key = "metric", value = "adeno",adeno:AR14)
adeno_agg[complete.cases(adeno_agg),]

temp_adeno <- adeno_agg %>% 
  rename(val = adeno) %>%
  filter(num_days >=15 & num_days <= day_max)

sgtrends_adeno = adeno_df_scaled%>% 
  left_join(reference_dates)%>%
  dplyr::select(date,all_of(unique(c(var_regressors_GT,var_regressors_GT_COVID))),
                num_days,epi_year,epi_week)%>%
  pivot_longer(cols = all_of(unique(c(var_regressors_GT,var_regressors_GT_COVID))),
               names_to = "metric",values_to = "val")%>%
  filter(num_days >=15 & num_days <= day_max)
unique(sgtrends_adeno$metric)

sg = unique(sgtrends_adeno$metric)
sg

col_order = c("date","num_days","epi_year","epi_week", "adeno", sprintf("AR%02d", 1:14), sg)
glm_adeno <- rbind(temp_adeno, sgtrends_adeno)
glm_adeno <- spread(glm_adeno, key = "metric", value = "val")
glm_adeno <- glm_adeno[col_order]
head(glm_adeno)
tail(glm_adeno)

res_GT_COVID = data.frame()
res_GT = data.frame()
res_AR = data.frame()
res_ARGO = data.frame()
for (t_predict_start in t0_predict_start:tmax_predict_start) {
  print(t_predict_start)
  t_predict_end <- t_predict_start + dt_predict - 1
  t_train_start <- t_predict_start - window_train 
  t_train_end <- t_predict_start - 1
  date_train_end <- reference_dates[reference_dates$num_days == t_train_end, c("date")]
  
  glm_train <- filter(glm_adeno, num_days >= t_train_start & num_days <= t_train_end)
  glm_predict <- filter(glm_adeno, num_days >= t_predict_start & num_days <= t_predict_end)
  
  x_train_GT <- data.matrix(glm_train[, var_regressors_GT])
  x_train_GT_COVID <- data.matrix(glm_train[, var_regressors_GT_COVID])
  x_train_AR <- data.matrix(glm_train[, var_regressors_AR])
  x_train_ARGO <- data.matrix(glm_train[, var_regressors_ARGO])
  
  x_predict_GT <- data.matrix(glm_predict[, var_regressors_GT])
  x_predict_GT_COVID <- data.matrix(glm_predict[, var_regressors_GT_COVID])
  x_predict_AR <- data.matrix(glm_predict[, var_regressors_AR])
  x_predict_ARGO <- data.matrix(glm_predict[, var_regressors_ARGO])
  
  y_train <- glm_train[, var_target] %>% pull(adeno)
  y_train = as.numeric(y_train)
  
  cvfit_GT <- cv.glmnet(x_train_GT, y_train,weights=w,type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_GT_COVID <- cv.glmnet(x_train_GT_COVID, y_train,weights=w,type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_AR <- cv.glmnet(x_train_AR, y_train,weights=w, type.measure = "mse", nfolds = nfolds, alpha = alpha)
  cvfit_ARGO <- cv.glmnet(x_train_ARGO, y_train, weights=w,type.measure = "mse", nfolds = nfolds, alpha = alpha)
  
  coef_GT <- cvfit_GT %>% coef(s = s) %>% as.matrix %>% t %>% as.data.frame %>% rename(intercept = "(Intercept)")
  coef_GT_COVID <- cvfit_GT_COVID %>% coef(s = s) %>% as.matrix %>% t %>% as.data.frame %>% rename(intercept = "(Intercept)")
  coef_AR <- coef(cvfit_AR, s = s) %>% as.matrix %>% t %>% as.data.frame %>% rename(intercept = "(Intercept)")
  coef_ARGO <- coef(cvfit_ARGO, s = s) %>% as.matrix %>% t %>% as.data.frame %>% rename(intercept = "(Intercept)")
  
  predict_GT <- cvfit_GT %>% predict(newx = x_predict_GT, s = s) %>% t %>% data.frame
  predict_GT_COVID <- cvfit_GT_COVID %>% predict(newx = x_predict_GT_COVID, s = s) %>% t %>% data.frame
  predict_AR <- cvfit_AR %>% predict(newx = x_predict_AR, s = s) %>% t %>% data.frame
  predict_ARGO <- cvfit_ARGO %>% predict(newx = x_predict_ARGO, s = s) %>% t %>% data.frame
  colnames(predict_GT) <- var_predicts 
  colnames(predict_GT_COVID) <- var_predicts 
  colnames(predict_AR) <- var_predicts 
  colnames(predict_ARGO) <- var_predicts 
  
  data_epi <- glm_predict[, "adeno"] %>% t %>% data.frame 
  colnames(data_epi) <- var_data
  
  res_GT <- bind_rows(res_GT, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_GT, coef_GT))
  head(res_GT)
  
  res_GT_COVID <- bind_rows(res_GT_COVID, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_GT_COVID, coef_GT_COVID))
  head(res_GT_COVID)
  
  res_AR <- bind_rows(res_AR, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_AR, coef_AR))
  res_ARGO <- bind_rows(res_ARGO, cbind(data.frame("num_days" = t_train_end, "date" = date_train_end), data_epi, predict_ARGO, coef_ARGO))
}

########################################################################################
## Summarize model results
########################################################################################

head(res_ARGO)
tail(res_ARGO)

coef_ARGO <- select(res_ARGO, -c("num_days"), -ends_with("_ahead"))
coef_ARGO <- coef_ARGO %>% gather(key = "regressor", val = "coef", -date)
range(coef_ARGO$coef)
coef_ARGO$regressor <- gsub("_",".", coef_ARGO$regressor)
orders = rev(c("intercept", gsub("_", ".", var_regressors_ARGO)))

# look at coefficients by correlation
summary_coef <- coef_ARGO %>% 
  group_by(regressor) %>% 
  summarise(avg = mean(coef, na.rm = TRUE)) %>% 
  ungroup()
summary_coef %>% arrange(-abs(avg)) %>% print(n=30)
summary_coef[summary_coef$avg==0,]

summary_coef <- coef_ARGO %>% 
  group_by(regressor) %>% 
  summarise(avg = mean(abs(coef), na.rm = TRUE)) %>% 
  ungroup()
summary_coef %>% arrange(-abs(avg)) %>% filter(!grepl("AR",regressor)) %>% print(n=30)
summary_coef[summary_coef$avg==0,]

summary_coef <- coef_ARGO %>% 
  filter(coef!=0)%>%
  group_by(regressor) %>% 
  tally() %>% arrange(-n)
summary_coef %>% head(n=10)

########################################################################################
## Plots
########################################################################################

########################################################################################
## Plot coefficient values at each time point
########################################################################################

patterns <- c("AR","intercept")
keep <- unique(coef_ARGO$regressor)[!grepl(paste(patterns, collapse="|"), unique(coef_ARGO$regressor))]

ar_terms = unique(coef_ARGO$regressor)[c(1:15)]
coef_ARGO %>% filter(regressor %in% ar_terms) %>% pull(coef)%>% range

label_df = data.frame(var = unique(coef_ARGO$regressor),
                      label = c(unique(coef_ARGO$regressor)[1:15],
                                "within-neighborhood movement",
                                "between-neighborhood movement",
                                "influx visitors other WA counties",
                                "influx out-of-state visitors",
                                "% devices leaving home",
                                "restaurants",
                                "groceries and pharmacies",
                                "transit",
                                # "hospitals",
                                "religious organizations",
                                "child daycare",
                                "elementary and high schools",
                                "colleges",
                                "COVID Rt"))

coef_ARGO <- select(res_ARGO, -c("num_days"), -ends_with("_ahead"))
coef_ARGO <- coef_ARGO %>% gather(key = "regressor", val = "coef", -date)
range(coef_ARGO$coef)

coef_ARGO[coef_ARGO$coef < -0.006, "coef"] <- -0.006
coef_ARGO[coef_ARGO$coef > 0.006, "coef"] <- 0.006
coef_ARGO$regressor <- gsub("_",".", coef_ARGO$regressor)
plot_df = left_join(coef_ARGO,label_df,by=c("regressor"="var"))

orders = c(unique(plot_df$label)[1:15],
           "within-neighborhood movement",
           "between-neighborhood movement",
           "influx visitors other WA counties",
           "influx out-of-state visitors",
           "% devices leaving home",
           "restaurants",
           "groceries and pharmacies",
           "transit",
           # "hospitals",
           "religious organizations",
           "child daycare",
           "elementary and high schools",
           "colleges",
           "COVID Rt")   

p <-ggplot(plot_df %>% filter(regressor %in% unique(coef_ARGO$regressor)[16:28]), 
           aes(x = date, y = factor(label, level = orders), fill = coef)) +
  geom_tile(color= "white",linewidth=0.1) + 
  scale_fill_gradient2(name="Coefficient",high = scales::muted("red"), mid = "white",
                       low= scales::muted("blue"))+
  ylab("")+xlab("Date")+
  scale_x_date(expand=c(0,0),date_labels = "%b %y",date_breaks = "4 months")+
  theme_minimal(base_size = 18)
p

p <-p + theme(legend.position = "right")+
  theme(strip.background = element_rect(colour="white"))+
  theme(legend.title=element_text(size=14))+
  theme(legend.text=element_text(size=10))+
  # ggtitle("adenovirus AR-SG: time-varying coefficients") +
  labs(x = "", y = "") +
  removeGrid()#ggExtra
p

coef_ARGO <- select(res_ARGO, -c("num_days"), -ends_with("_ahead"))
coef_ARGO <- coef_ARGO %>% gather(key = "regressor", val = "coef", -date)
range(coef_ARGO$coef)
coef_ARGO[coef_ARGO$coef < -0.08, "coef"] <- -0.08
coef_ARGO[coef_ARGO$coef > 0.08, "coef"] <- 0.08
coef_ARGO$regressor <- gsub("_",".", coef_ARGO$regressor)
plot_df = left_join(coef_ARGO,label_df,by=c("regressor"="var"))

q <-ggplot(plot_df %>% filter(!regressor %in% keep), 
           aes(x = date, y = factor(label, level = rev(orders)), fill = coef)) +
  geom_tile(color= "white",linewidth=0.1) + 
  scale_fill_gradient2(name="Coefficient",high = scales::muted("red"), mid = "white",
                       low= scales::muted("blue"))+
  ylab("")+xlab("Date")+
  scale_x_date(expand=c(0,0),date_labels = "%b %y",date_breaks = "4 months")+
  theme_minimal(base_size = 18)
q

q <-q + theme(legend.position = "right")+
  theme(strip.background = element_rect(colour="white"))+
  theme(legend.title=element_text(size=14))+
  theme(legend.text=element_text(size=10))+
  # ggtitle("adenovirus AR-SG: time-varying coefficients") +
  labs(x = "", y = "") +
  removeGrid()#ggExtra
q
plot_grid(q,p,nrow=2,align = "v")

########################################################################################
## Figure S20: Overlay model predictions on observed AdV Rt and incidence
########################################################################################

res_ARGO2 <- res_ARGO %>% 
  mutate(data_7d_ahead = as.numeric(as.character(data_7d_ahead)))

ci_df = rt_df %>% filter(organism=="Adenovirus") %>%
  dplyr::select(date,median,lower,upper,level)%>%
  filter(level==90)
ci_df$date = as.Date(ci_df$date)

res_GT[!complete.cases(res_GT$predict_7d_ahead),]

full_df = readRDS("Epidemia_Models/Epidemia_Models_Biowulf/extended_daily_ILI_per_pos_scaled_for_select_pathogens_2023_10_22.rds")

adeno_results = full_df %>% filter(organism =="Adenovirus")%>%droplevels()

adeno_results = adeno_results %>%
  arrange(organism,date)%>%
  group_by(organism)%>%
  mutate(all_sites_per_pos_ili_scaled_sum_mv_avg = 
           zoo::rollmean(all_sites_per_pos_ili_scaled_sum,k=15,align="center",
                         fill=NA))%>%
  ungroup()

r = ggplot() +
  geom_hline(aes(yintercept=1),lty="dashed")+
  geom_rect(aes(xmin = as.Date("2019-01-01"), 
                xmax = as.Date("2019-02-21"), ymin = -Inf, ymax = Inf),
            alpha = 0.1,
            fill = "yellow")+
  geom_rect(aes(xmin = as.Date("2020-03-23"), 
                xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf),
            alpha = 0.1,
            fill = "orange")+
  geom_line(data=adeno_results %>%
              filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")),
            aes(x=date,y=(100*all_sites_per_pos_ili_scaled_sum_mv_avg)/40),
            alpha=0.8,lwd=1,lty="solid",color="#117733")+
  geom_ribbon(data=adeno_results %>% 
                filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")),
              aes(x=date,ymin=0,ymax=(100*all_sites_per_pos_ili_scaled_sum_mv_avg)/40),
              fill="#117733",alpha=0.4)+
  geom_line(data = res_GT %>% 
              filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")) , 
            aes(x = date, y = predict_7d_ahead, colour = "SG"), linewidth = 1,alpha=0.9) +
  geom_line(data = res_GT_COVID %>% 
              filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")) , 
            aes(x = date, y = predict_7d_ahead, colour = "SG-COVID"), linewidth = 1,alpha=0.9) +
  geom_line(data = res_ARGO %>%
              filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")), 
            aes(x = date,  y = predict_7d_ahead, colour = "AR-SG-COVID"), 
            linewidth = 1,alpha=0.9)+
  geom_line(data = res_AR %>% 
              filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")), 
            aes(x = date , y = predict_7d_ahead, colour = "AR"), linewidth = 1,alpha=0.9) +
  geom_line(data = ci_df %>% 
              filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")),
            aes(x = date, y=median,colour = "Adeno Rt"),linewidth=1.2) +
  geom_ribbon(data = ci_df %>% 
                filter(date>=as.Date("2019-01-01")&date<=as.Date("2022-05-01")),
              aes(x = date, ymin = lower,ymax=upper),fill="#ABDDA4",alpha=0.5) +
  labs(x = "Date", y = "Rt", color = "Model") +
  theme_classic(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),legend.position = "bottom"
  ) + 
  geom_vline(aes(xintercept=as.Date("2019-02-21")),lty="dashed",color="yellow")+
  scale_color_manual(name=NULL,
                     values=c("#ABDDA4","#D7191C","#2B83BA","purple", "#FDAE61"),
                     breaks=c("Adeno Rt","AR-SG-COVID","SG-COVID","SG","AR"),
                     labels=c("AdV Rt","AR-Mobility-COVID","Mobility-COVID","Mobility","AR"))+
  background_grid()+
  scale_x_date(expand=c(0,0),date_breaks = "3 months",date_labels = "%b %y")+
  scale_y_continuous(
    # Features of the first axis
    name = "Rt",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*40, name="Incidence"),
    expand=c(0,0)
  )
r
# save_plot(r,filename = "figures/fig_s20_adenovirus_rt_predictive_model_fit_7d_ahead.png",base_width = 16,base_height = 6)

########################################################################################
## Predictive accuracy
########################################################################################

rmse(predicted=res_ARGO$predict_7d_ahead, actual=res_ARGO2$data_7d_ahead)# 0.035
rmse(predicted=res_GT_COVID$predict_7d_ahead,actual=res_ARGO2$data_7d_ahead)#0.174
rmse(predicted=res_GT$predict_7d_ahead, actual=res_ARGO2$data_7d_ahead)#0.174
rmse(predicted=res_AR$predict_7d_ahead, actual=res_ARGO2$data_7d_ahead)# 0.029

mae(predicted=res_ARGO$predict_7d_ahead, actual=res_ARGO2$data_7d_ahead)# 0.026
mae(predicted=res_GT_COVID$predict_7d_ahead,actual=res_ARGO2$data_7d_ahead)#0.126
mae(predicted=res_GT$predict_7d_ahead, actual=res_ARGO2$data_7d_ahead)#0.125
mae(predicted=res_AR$predict_7d_ahead, actual=res_ARGO2$data_7d_ahead)# 0.021

MAPE(res_ARGO$predict_7d_ahead, res_ARGO2$data_7d_ahead)#0.028
MAPE(res_GT_COVID$predict_7d_ahead, res_ARGO2$data_7d_ahead)#0.132
MAPE(res_GT$predict_7d_ahead, res_ARGO2$data_7d_ahead)#0.132
MAPE(res_AR$predict_7d_ahead, res_ARGO2$data_7d_ahead)#0.022

cor(res_ARGO$predict_7d_ahead, res_ARGO2$data_7d_ahead,  method = "pearson", use = "complete.obs")#0.99
cor(res_GT_COVID$predict_7d_ahead, res_ARGO2$data_7d_ahead,  method = "pearson", use = "complete.obs")#0.83
cor(res_GT$predict_7d_ahead, res_ARGO2$data_7d_ahead,  method = "pearson", use = "complete.obs")#0.82
cor(res_AR$predict_7d_ahead, res_ARGO2$data_7d_ahead,  method = "pearson", use = "complete.obs")#0.99

## mobility only
predicted <- res_GT$predict_7d_ahead
actual <- res_ARGO2$data_7d_ahead
error <- actual - predicted

GT_res_df = data.frame(date = res_GT$date,
                       predicted = predicted,
                       model="Mobility",
                       actual = actual , 
                       error = error,
                       year_mon = as.yearmon(res_GT$date))
head(GT_res_df)

## mobility + COVID Rt
predicted <- res_GT_COVID$predict_7d_ahead
actual <- res_ARGO2$data_7d_ahead
error <- actual - predicted

GT_COVID_res_df = data.frame(date = res_GT_COVID$date,
                             predicted = predicted,
                             model="Mobility-COVID",
                             actual = actual , 
                             error = error,
                             year_mon = as.yearmon(res_GT_COVID$date))

## AR only
predicted <- res_AR$predict_7d_ahead
actual <- res_ARGO2$data_7d_ahead
error <- actual - predicted
AR_res_df = data.frame(date = res_AR$date,
                       model = "AR",
                       predicted = predicted,
                       actual = actual , 
                       error = error,
                       year_mon = as.yearmon(res_AR$date))

## full model
predicted <- res_ARGO$predict_7d_ahead
actual <- res_ARGO2$data_7d_ahead
error <- actual - predicted
ARGO_res_df = data.frame(date = res_ARGO$date,
                         model="AR-Mobility-COVID",
                         predicted = predicted,
                         actual = actual , 
                         error = error,
                         year_mon = as.yearmon(res_ARGO$date))

all_resid = bind_rows(GT_res_df,AR_res_df,ARGO_res_df, GT_COVID_res_df)
all_resid$model = factor(all_resid$model,levels=c("AR-Mobility-COVID","Mobility-COVID","Mobility","AR"))


p = ggplot()+
  geom_rect(aes(xmin = as.Date("2019-01-01"), 
                xmax = as.Date("2019-02-21"), ymin = -Inf, ymax = Inf),
            alpha = 0.1,
            fill = "yellow")+
  geom_rect(aes(xmin = as.Date("2020-03-23"), 
                xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf),
            alpha = 0.1,
            fill = "orange")+
  geom_line(aes(x=date,y=error,color=model,group=model),
            data=all_resid %>% filter(date<as.Date("2022-05-01")),lwd=0.8)+
  scale_x_date(expand=c(0,0),date_breaks = "3 months",date_labels = "%b %y",limits = c(as.Date("2019-01-01"),as.Date("2022-04-30")))+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom")+
  scale_color_manual(name="Model",
                     values=c("#D7191C","#2B83BA","purple", "#FDAE61"),
                     breaks=c("AR-Mobility-COVID","Mobility-COVID","Mobility","AR"))+
  xlab("Date")+
  ylab("Estimation error")+
  ylim(-1,1)
p
# save_plot(p,filename="figures/estimation_error_adeno_models.png",base_width = 12,base_height = 4)

com = plot_grid(r,p,nrow=2,labels="AUTO",align = "v")
com
save_plot(com, filename = "figures/figure_s20_adeno_predictive_model_fit_and_error.png",base_width = 12,base_height = 10)


##### stay at home and initial rebound
res_ARGO_sah = res_ARGO %>% filter(date>=as.Date("2020-02-28")&date<as.Date("2020-06-30"))
res_AR_sah = res_AR %>% filter(date>=as.Date("2020-02-28")&date<as.Date("2020-06-30"))
res_GT_sah = res_GT %>% filter(date>=as.Date("2020-02-28")&date<as.Date("2020-06-30"))
res_GT_COVID_sah = res_GT_COVID %>% filter(date>=as.Date("2020-02-28")&date<as.Date("2020-06-30"))
res_ARGO2_sah = res_ARGO2 %>% filter(date>=as.Date("2020-02-28")&date<as.Date("2020-06-30"))
range(res_ARGO2_sah$date)

rmse(res_ARGO_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)# 0.021
rmse(res_GT_COVID_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.096
rmse(res_GT_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.09
rmse(res_AR_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.019

mae(res_ARGO_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)# 0.017
mae(res_GT_COVID_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.077
mae(res_GT_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.072
mae(res_AR_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.014

MAPE(res_ARGO_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.025
MAPE(res_GT_COVID_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.106
MAPE(res_GT_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.099
MAPE(res_AR_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead)#0.021

cor(res_ARGO_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead,  method = "pearson", use = "complete.obs")#1.0
cor(res_GT_COVID_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead,  method = "pearson", use = "complete.obs")#0.95
cor(res_GT_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead,  method = "pearson", use = "complete.obs")#0.95
cor(res_AR_sah$predict_7d_ahead, res_ARGO2_sah$data_7d_ahead,  method = "pearson", use = "complete.obs")#1.0