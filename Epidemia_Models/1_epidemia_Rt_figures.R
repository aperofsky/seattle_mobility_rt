## Input data are provided in the Epidemia_Models/epidemia_rt_output folder, so it's not necessary to run the Epidemia models to run this code
###################################################
## Combine Rt estimates from all pathogens
## Compare SCAN COVID incidence vs King County cases
## Plot Rt and incidence (Figure 1)
## Save all Rt combined for downstream analyses
###################################################

library(dplyr)
library(tidyr)
library(readr)

library(ggplot2)
library(cowplot)
library(jcolors)
library(tidyquant)

###################################################
##### RHINO
###################################################
all_rhino_rt <- read_rds("Epidemia_Models/epidemia_rt_output/Rhino_2018_2022_biowulf_epidemia_model_rt_only.rds")
head(all_rhino_rt)
###################################################
##### Adenovirus
###################################################
all_adeno_rt <- read_rds("Epidemia_Models/epidemia_rt_output/Adeno_2018_2022_biowulf_epidemia_model_rt_only.rds")
head(all_adeno_rt)
###################################################
##### RSV A
###################################################
all_rsv_a_rt <- read_rds("Epidemia_Models/epidemia_rt_output/RSV_A_biowulf_epidemia_model_rt_only.rds")
head(all_rsv_a_rt)
###################################################
##### RSV B
###################################################
all_rsv_b_rt <- read_rds("Epidemia_Models/epidemia_rt_output/RSV_B_biowulf_epidemia_model_rt_only.rds") %>% filter(date <= as.Date("2022-02-20"))
head(all_rsv_b_rt)
###################################################
##### Flu A/H3N2
###################################################
all_H3N2_rt <- read_rds("Epidemia_Models/epidemia_rt_output/H3N2_biowulf_epidemia_model_rt_only.rds")
head(all_H3N2_rt)
###################################################
##### Flu A/H1N1
###################################################
all_H1N1_rt <- read_rds("Epidemia_Models/epidemia_rt_output/H1N1_biowulf_epidemia_model_rt_only.rds")
head(all_H1N1_rt)
###################################################
##### Flu B
###################################################
all_IVB_rt <- read_rds("Epidemia_Models/epidemia_rt_output/IBV_2020_biowulf_epidemia_model_rt_only.rds")
head(all_IVB_rt)
###################################################
##### HPIV 3 + 4
###################################################
all_HPIV_3_4_rt <- read_rds("Epidemia_Models/epidemia_rt_output/HPIV_3_4_biowulf_epidemia_model_rt_only.rds")
head(all_HPIV_3_4_rt)
###################################################
##### HPIV 1 + 2
###################################################
all_HPIV_1_2_rt <- read_rds("Epidemia_Models/epidemia_rt_output/HPIV_1_2_biowulf_epidemia_model_rt_only.rds")
head(all_HPIV_1_2_rt)
###################################################
##### seasonal CoV 229E + OC43
###################################################
all_seasonal_CoV_229E_OC43_rt <- read_rds("Epidemia_Models/epidemia_rt_output/seasonal_CoV_229E_OC43_biowulf_epidemia_model_rt_only.rds")
head(all_seasonal_CoV_229E_OC43_rt)
###################################################
##### seasonal CoV 229E + OC43
###################################################
all_seasonal_CoV_HKU1_NL63_rt <- read_rds("Epidemia_Models/epidemia_rt_output/seasonal_CoV_HKU1_NL63_biowulf_epidemia_model_rt_only.rds")
head(all_seasonal_CoV_HKU1_NL63_rt)
###################################################
##### hMPV
###################################################
all_hMPV_rt <- read_rds("Epidemia_Models/epidemia_rt_output/hMPV_biowulf_epidemia_model_rt_only.rds")
head(all_hMPV_rt)
###################################################
## combine all pathogen Rt estimates
###################################################

all_endemic_results <- bind_rows(
  all_rsv_a_rt, all_rsv_b_rt, all_H1N1_rt, all_H3N2_rt,
  all_IVB_rt, all_rhino_rt, all_HPIV_3_4_rt, all_HPIV_1_2_rt,
  all_hMPV_rt, all_adeno_rt, all_seasonal_CoV_229E_OC43_rt,
  all_seasonal_CoV_HKU1_NL63_rt
)
head(all_endemic_results)
unique(all_endemic_results$organism)

all_endemic_results %>%
  group_by(organism, date) %>%
  tally() %>%
  filter(n != 1)

all_endemic_results$organism <- as.factor(all_endemic_results$organism)
levels(all_endemic_results$organism)[c(5:8)]
levels(all_endemic_results$organism)[c(5:8)] <- c("IAV/H1N1", "IAV/H3N2", "IBV", "RV")

levels(all_endemic_results$organism)
all_endemic_results$organism <- factor(all_endemic_results$organism)
levels(all_endemic_results$organism)
levels(all_endemic_results$organism) <- c(
  "Adenovirus", "Human metapneumovirus",
  "Human parainfluenza 1 + 2", "Human parainfluenza 3 + 4",
  "Influenza A/H1N1", "Influenza A/H3N2", "Influenza B",
  "Rhinovirus",
  "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B",
  "Seasonal CoV 229E + OC43", "Seasonal CoV HKU1 + NL63"
)

all_endemic_results_avg <- all_endemic_results %>%
  group_by(group, tag, level, organism) %>%
  mutate_at(
    c("lower", "upper", "median"),
    ~ zoo::rollmean(., k = 15, align = "center", fill = NA)
  ) %>%
  ungroup()

all_endemic_results_avg <- all_endemic_results_avg %>% filter(date < as.Date("2022-07-01"))

###################################################
## endemic pathogen df
###################################################
full_df <- read_rds("Epidemia_Models/Epidemia_Models_Biowulf/extended_daily_ILI_per_pos_scaled_for_select_pathogens_2023_10_22.rds")
head(full_df)
unique(full_df$organism)

full_df$organism <- as.factor(full_df$organism)
range(full_df$date)
unique(full_df$organism)
levels(full_df$organism)
levels(full_df$organism) <- c(
  "Adenovirus", "Seasonal CoV 229E + OC43", "Seasonal CoV HKU1 + NL63",
  "Enterovirus", "Human bocavirus", "Human metapneumovirus",
  "Human parainfluenza 1 + 2",
  "Human parainfluenza 3 + 4", "Influenza A/H1N1", "Influenza A/H3N2", "Influenza B",
  "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B", "Rhinovirus", "SARS-CoV-2",
  "Streptococcus pneumoniae"
)
full_df$organism <- factor(full_df$organism, levels = c(
  "Influenza A/H1N1", "Influenza A/H3N2", "Influenza B",
  "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B",
  "Human metapneumovirus",
  "Seasonal CoV 229E + OC43", "Seasonal CoV HKU1 + NL63",
  "Human parainfluenza 1 + 2", "Human parainfluenza 3 + 4",
  "Rhinovirus",
  "Enterovirus",
  "Adenovirus",
  "Human bocavirus",
  "SARS-CoV-2",
  "Streptococcus pneumoniae"
))

unique(full_df$organism)

full_df <- full_df %>%
  arrange(organism, date) %>%
  group_by(organism) %>%
  mutate(
    all_sites_per_pos_ili_scaled_mean_mv_avg = zoo::rollmean(all_sites_per_pos_ili_scaled_mean,
      k = 15, align = "center",
      fill = NA
    ),
    all_sites_per_pos_ili_scaled_sum_mv_avg = zoo::rollmean(all_sites_per_pos_ili_scaled_sum,
      k = 15, align = "center",
      fill = NA
    )
  ) %>%
  ungroup() %>%
  filter(date < as.Date("2022-07-01"))

# covid-19 incidence based on SFS/SCAN data
sfs_covid_incidence <- full_df %>%
  filter(organism == "SARS-CoV-2") %>%
  dplyr::select(date, comm_all_ages_adj_per_pos_ili_scaled, hosp_all_ages_adj_per_pos_ili_scaled, all_sites_per_pos_ili_scaled_sum) %>%
  filter(date > as.Date("2020-01-15"))

sfs_covid_incidence_long <-
  sfs_covid_incidence %>%
  pivot_longer(cols = c(comm_all_ages_adj_per_pos_ili_scaled:all_sites_per_pos_ili_scaled_sum), names_to = "incidence_measure", values_to = "incidence")

###################################################
## King County COVID-19 cases
###################################################

covid_cases <- read_csv("Epidemia_Models/Epidemia_Models_Biowulf/wa_doh_king_county_cases.csv")
head(covid_cases)

covid_cases <- covid_cases %>%
  rename(
    total_cases = Total.Cases.7.Day.Count,
    date = Specimen.Collection.Date
  )

covid_cases_all_voc <- covid_cases %>%
  dplyr::select(date, total_cases) %>%
  distinct() %>%
  filter(date < as.Date("2022-07-01"))

###################################################
## Fig. S22 Comparison between SCAN COVID incidence vs county cases
###################################################

pre_omicron_comparison <- ggplot() +
  geom_ribbon(aes(x = date, ymin = 0, ymax = incidence, fill = incidence_measure),
    alpha = 0.1,
    data = sfs_covid_incidence_long %>%
      filter(incidence_measure == "comm_all_ages_adj_per_pos_ili_scaled") %>% filter(date < as.Date("2021-12-01"))
  ) +
  geom_ma(aes(x = date, y = incidence, color = incidence_measure), data = sfs_covid_incidence_long %>%
    filter(incidence_measure == "comm_all_ages_adj_per_pos_ili_scaled") %>% filter(date < as.Date("2021-12-01")), lty = "solid", n = 7, size = 2) +
  geom_ribbon(aes(x = date, ymin = 0, ymax = total_cases / 65000, fill = "cases"),
    alpha = 0.1,
    data = covid_cases_all_voc %>% filter(date < as.Date("2021-12-01"))
  ) +
  geom_line(aes(x = date, y = total_cases / 65000, color = "cases"), data = covid_cases_all_voc %>% filter(date < as.Date("2021-12-01")), lty = "solid") +
  scale_y_continuous(
    # Features of the first axis
    name = "SCAN COVID-19 Incidence",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * 65000, name = "King County COVID-19 Cases"),
    expand = c(0, 0)
  ) +
  scale_color_jcolors(
    palette = "pal8", breaks = c("cases", "comm_all_ages_adj_per_pos_ili_scaled"),
    labels = c("King County COVID-19 Cases", "SCAN COVID-19 Community Incidence"), name = "Incidence Measure"
  ) +
  scale_fill_jcolors(
    palette = "pal8", breaks = c("cases", "comm_all_ages_adj_per_pos_ili_scaled"),
    labels = c("King County COVID-19 Cases", "SCAN COVID-19 Community Incidence"), name = "Incidence Measure"
  ) +
  theme_bw(base_size = 16) +
  scale_x_date(expand = c(0, 0), date_labels = "%b-%y", date_breaks = "3 months") +
  theme(legend.position = "bottom")
pre_omicron_comparison <- pre_omicron_comparison + ggtitle("Pre-Omicron Period")
pre_omicron_comparison

post_omicron_comparison <- ggplot() +
  geom_ribbon(aes(x = date, ymin = 0, ymax = incidence, fill = incidence_measure),
    alpha = 0.1,
    data = sfs_covid_incidence_long %>%
      filter(incidence_measure == "comm_all_ages_adj_per_pos_ili_scaled") %>%
      filter(date >= as.Date("2021-12-01") & date < as.Date("2022-07-01"))
  ) +
  geom_ma(aes(x = date, y = incidence, color = incidence_measure),
    data = sfs_covid_incidence_long %>%
      filter(incidence_measure == "comm_all_ages_adj_per_pos_ili_scaled") %>% filter(date >= as.Date("2021-12-01") & date < as.Date("2022-07-01")),
    lty = "solid", n = 7, size = 2
  ) +
  geom_ribbon(aes(x = date, ymin = 0, ymax = total_cases / 65000, fill = "cases"),
    alpha = 0.1,
    data = covid_cases_all_voc %>% filter(date >= as.Date("2021-12-01") & date < as.Date("2022-07-01"))
  ) +
  geom_line(aes(x = date, y = total_cases / 65000, color = "cases"), data = covid_cases_all_voc %>%
    filter(date >= as.Date("2021-12-01") & date < as.Date("2022-07-01")), lty = "solid") +
  scale_y_continuous(
    # Features of the first axis
    name = "SCAN COVID-19 Incidence",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * 65000, name = "King County COVID-19 Cases"),
    expand = c(0, 0)
  ) +
  scale_color_jcolors(
    palette = "pal8", breaks = c("cases", "comm_all_ages_adj_per_pos_ili_scaled"),
    labels = c("King County COVID-19 Cases", "SCAN COVID-19 Community Incidence"), name = "Incidence Measure"
  ) +
  scale_fill_jcolors(
    palette = "pal8", breaks = c("cases", "comm_all_ages_adj_per_pos_ili_scaled"),
    labels = c("King County COVID-19 Cases", "SCAN COVID-19 Community Incidence"), name = "Incidence Measure"
  ) +
  theme_bw(base_size = 16) +
  scale_x_date(expand = c(0, 0), date_labels = "%b-%y", date_breaks = "3 months") +
  theme(legend.position = "bottom")
post_omicron_comparison <- post_omicron_comparison +
  ggtitle("Omicron Period")
post_omicron_comparison

covid_comparison_plot <- plot_grid(pre_omicron_comparison + theme(legend.position = "none"), post_omicron_comparison + theme(legend.position = "none"), nrow = 1, rel_widths = c(0.75, 0.5))
leg <- get_legend(post_omicron_comparison)
covid_comparison2 <- plot_grid(covid_comparison_plot, leg, nrow = 2, rel_heights = c(1, 0.1))
covid_comparison2
save_plot(covid_comparison2, filename = "figures/fig_s22_scan_vs_county_covid_comparison.png", base_width = 20, base_height = 8)

###################################################
## Plot Rt and incidence
###################################################

## limit incidence dataset to pathogens with Rt estimates
full_df2 <- full_df %>%
  filter(organism %in% levels(all_endemic_results_avg$organism)) %>%
  droplevels()

unique(full_df2$organism)
setdiff(unique(full_df2$organism), unique(all_endemic_results_avg$organism))
setdiff(unique(all_endemic_results_avg$organism), unique(full_df2$organism))

range(full_df2$date)

###################################################
# endemic pathogen incidence and Rt
###################################################

color_vec <- c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
length(color_vec)

plot_coef <- 20
all_endemic_plot <- ggplot(full_df2 %>% filter(date > as.Date("2018-12-01"))) +
  geom_rect(aes(xmin = as.Date("2019-02-03"), xmax = as.Date("2019-02-15"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "blue", data = data.frame(x = 0, y = 0)
  ) +
  geom_rect(aes(xmin = as.Date("2020-03-23"), xmax = as.Date("2020-06-05"), ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "orange", data = data.frame(x = 0, y = 0)
  ) +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen") +
  geom_line(
    data = full_df2,
    aes(x = date, y = 100 * all_sites_per_pos_ili_scaled_sum_mv_avg / plot_coef, color = organism),
    alpha = 0.5, lwd = 1.1, lty = "solid"
  ) +
  geom_ribbon(
    data = full_df2,
    aes(x = date, ymin = 0, ymax = 100 * all_sites_per_pos_ili_scaled_sum_mv_avg / plot_coef, fill = organism), alpha = 0.1
  ) +
  geom_hline(yintercept = 1, lty = "dashed", color = "#004488", lwd = 1) +
  geom_line(
    data = all_endemic_results_avg %>%
      dplyr::select(date, median, organism) %>%
      distinct() %>%
      filter(date > as.Date("2018-12-01")) %>%
      filter(!(organism == "Influenza A/H3N2" & date > as.Date("2021-01-01"))),
    aes(x = date, y = median, color = organism), lwd = 1.5, lty = "solid"
  ) +
  geom_ribbon(
    data = all_endemic_results_avg %>%
      filter(date > as.Date("2018-12-01")) %>%
      dplyr::select(date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>%
      filter(!(organism == "Influenza A/H3N2" & date > as.Date("2021-01-01"))),
    aes(x = date, ymin = lower, ymax = upper, fill = organism), alpha = 0.75
  ) +
  theme_bw(base_size = 16) +
  scale_x_date(expand = c(0, 0)) +
  facet_wrap(~organism, scales = "free_y", ncol = 2, dir = "v") +
  scale_y_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6),
    # Features of the first axis
    name = "Rt",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * plot_coef, name = "Scaled Incidence"),
    expand = c(0, 0)
  ) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 16)) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = color_vec)
all_endemic_plot

###################################################
# Endemic pathogen incidence and Rt for 2018-2019 season (snowstorm figure)
###################################################

full_df2 <- full_df %>%
  filter(date < as.Date("2019-05-01") & date > as.Date("2018-12-20")) %>%
  droplevels()

full_df2$organism <- factor(full_df2$organism, levels = c(
  "Influenza A/H1N1", "Influenza A/H3N2", "Influenza B",
  "Seasonal CoV 229E + OC43", "Seasonal CoV HKU1 + NL63",
  "Human metapneumovirus",
  "Human parainfluenza 1 + 2", "Human parainfluenza 3 + 4",
  "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B",
  "Rhinovirus",
  "Enterovirus",
  "Adenovirus",
  "Human bocavirus",
  "SARS-CoV-2",
  "Streptococcus pneumoniae"
))

all_endemic_results$organism <- factor(all_endemic_results$organism, levels = c(
  "Influenza A/H1N1", "Influenza A/H3N2",
  "Influenza B", "Respiratory syncytial virus (RSV) A",
  "Respiratory syncytial virus (RSV) B", "Human metapneumovirus",
  "Seasonal CoV 229E + OC43", "Seasonal CoV HKU1 + NL63",
  "Human parainfluenza 1 + 2", "Human parainfluenza 3 + 4",
  "Adenovirus", "Rhinovirus", "SARS-CoV-2"
))

keep <- all_endemic_results %>%
  filter(date < as.Date("2019-05-01") & date > as.Date("2018-12-20")) %>%
  pull(organism) %>%
  unique() %>%
  as.character()
keep

color_vec <- c(
  "#332288",
  "#6699CC",
  "#88CCEE",
  "#44AA99", "#117733", "#999933", "#DDCC77",
  "#CC6677", "#AA4466", "#882255", "#AA4499"
)

length(color_vec)
coeff <- 10
p <- ggplot() +
  geom_rect(aes(xmin = as.Date("2019-02-03"), xmax = as.Date("2019-02-15"), ymin = -Inf, ymax = Inf),
            alpha = 0.2,
            fill = "blue", data = data.frame(x = 0, y = 0)
  ) +
  geom_line(
    data = full_df2 %>%
      filter(date < as.Date("2019-05-01") & date > as.Date("2018-12-20")) %>%
      filter(organism %in% keep) %>% droplevels(),
    aes(x = date, y = 100 * all_sites_per_pos_ili_scaled_sum_mv_avg / coeff, color = organism),
    alpha = 0.5, lwd = 1, lty = "solid"
  ) +
  geom_ribbon(
    data = full_df2 %>%
      filter(date < as.Date("2019-05-01") & date > as.Date("2018-12-20")) %>%
      filter(organism %in% keep) %>%
      droplevels(),
    aes(x = date, ymin = 0, ymax = 100 * all_sites_per_pos_ili_scaled_sum_mv_avg / coeff, fill = organism), alpha = 0.1
  ) +
  geom_line(
    data = all_endemic_results %>%
      dplyr::select(date, median, organism) %>%
      distinct() %>%
      filter(date < as.Date("2019-05-01") & date > as.Date("2018-12-20")),
    aes(x = date, y = median, color = organism), lwd = 1.5, lty = "solid"
  ) +
  geom_ribbon(
    data = all_endemic_results %>%
      filter(date < as.Date("2019-05-01") & date > as.Date("2018-12-20")) %>%
      dplyr::select(date, level, organism, lower, upper) %>%
      filter(level == 90) %>%
      distinct(),
    aes(x = date, ymin = lower, ymax = upper, fill = organism), alpha = 0.5
  ) +
  scale_y_continuous(
    # Features of the first axis
    name = "Rt",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * coeff, name = "Incidence"),
    expand = c(0, 0)
  ) +
  geom_hline(yintercept = 1, lty = "dashed") +
  # ylab("Rt")+
  xlab("Date") +
  scale_color_manual(values = color_vec, name = "Pathogen") +
  scale_fill_manual(values = color_vec, name = "Pathogen") +
  theme_bw(base_size = 18) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 months", date_labels = "%b %Y") +
  theme(legend.position = "none", strip.background = element_blank()) +
  facet_wrap(~organism, scales = "free_y")
p
save_plot(p, filename = "figures/fig_s5_snowstorm_rt_closeup.png", base_width = 18, base_height = 12)

###################################################
##### COVID-19
###################################################

all_covid_rt <- read_rds("Epidemia_Models/epidemia_rt_output/COVID_2020_2022_biowulf_epidemia_model_rt_only.rds")
all_covid_rt$organism="SARS-CoV-2"
head(all_covid_rt)
all_covid_rt_avg = all_covid_rt %>%
  group_by(group,tag,level,organism)%>%
  mutate_at(c("lower","upper","median"),~zoo::rollmean(.,k=15,align="center", fill=NA))%>%ungroup()%>%
  filter(date < as.Date("2022-07-01"))

combined_mob <- read_rds("Seattle_Mobility_Data/mobility_data/mobility_metrics_for_epidemia.rds") %>% dplyr::select(date,epi_date,oxford_stringency_index)
head(combined_mob)
combined_mob$date <- as.Date(combined_mob$date)
combined_mob$epi_date <- as.Date(combined_mob$epi_date)
range(combined_mob$date)

## Rt from all cases combined
all_covid_plot <- ggplot() +
  geom_vline(xintercept = as.Date("2020-02-29"), lty = "dashed", color = "darkgreen") +
  geom_rect(
    aes(
      xmin = as.Date("2020-03-23"),
      xmax = as.Date("2020-06-05"),
      ymin = -Inf, ymax = Inf
    ),
    alpha = 0.2,
    fill = "orange", data = data.frame(x = 0, y = 0)
  ) +
  geom_line(
    data = combined_mob %>% filter(date >= as.Date("2020-02-20") & date < as.Date("2022-07-01")),
    aes(x = date, y = oxford_stringency_index / 10, color = "OSI"), alpha = 0.8, lwd = 1
  ) +
  geom_line(
    data = covid_cases %>% dplyr::select(date, total_cases) %>% distinct() %>% filter(date>=as.Date("2020-02-20") & date < as.Date("2022-07-01")),
    aes(x = date, y = total_cases / 5000, col = "All"),
    alpha = 0.8, lwd = 1, lty = "solid"
  ) +
  geom_ribbon(
    data = covid_cases %>% dplyr::select(date, total_cases) %>% distinct() %>% filter(date>=as.Date("2020-02-20") & date < as.Date("2022-07-01")),
    aes(x = date, ymin = 0, ymax = (total_cases) / 5000, fill = "All"), alpha = 0.2
  ) +
  geom_line(
    data = all_covid_rt %>%
      dplyr::select(date, median) %>% distinct() %>% filter(date>=as.Date("2020-02-20") & date < as.Date("2022-07-01")),
    aes(x = date, y = median, col = "All"), lwd = 1.2, lty = "solid"
  ) +
  geom_ribbon(
    data = all_covid_rt %>%
      dplyr::select(date, level, lower, upper) %>%
      filter(level == 90) %>%
      distinct() %>% filter(date>=as.Date("2020-02-20") & date < as.Date("2022-07-01")),
    aes(x = date, ymin = lower, ymax = upper, fill = "All"), alpha = 0.7
  ) +
  theme_bw(base_size = 16) +
  geom_hline(aes(yintercept = 1), lty = "dashed", color = "#004488", lwd = 1) +
  scale_y_continuous(
    # Features of the first axis
    name = "Rt or OSI",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ . * 5000, name = "COVID-19 cases"),
    expand = c(0, 0.1)
  ) +
  scale_fill_manual(values = c("#ABB065", "#DB9D85"), name = NULL, breaks = c("All", "OSI"), labels = c("SARS-CoV-2 Rt", "Oxford Stringency Index")) +
  scale_color_manual(values = c("#ABB065", "#DB9D85"), name = NULL, breaks = c("All", "OSI"), labels = c("SARS-CoV-2 Rt", "Oxford Stringency Index")) +
  guides(fill = "none") +
  theme(legend.position = c(0.65, 0.9), legend.background = element_blank()) +
  xlab("date") +
  scale_x_date(expand = c(0, 0), date_breaks = "4 months", date_labels = "%b %Y")
all_covid_plot

end_and_cov <- plot_grid(all_endemic_plot,
  all_covid_plot + ggtitle("SARS-CoV-2") + theme(plot.title = element_text(hjust = 0.5)),
  nrow = 2,
  labels = "AUTO", rel_heights = c(6, 3)
)
end_and_cov
save_plot(end_and_cov, filename = "figures/fig_1_rt_endemic_and_covid_combined.png", base_width = 14, base_height = 12)

#####################################################################
## Save combined Rt output
#####################################################################

## remake df to save with pathogen abbreviations
all_endemic_results <- bind_rows(
  all_rsv_a_rt, all_rsv_b_rt, all_H1N1_rt, all_H3N2_rt,
  all_IVB_rt, all_rhino_rt, all_HPIV_3_4_rt, all_HPIV_1_2_rt,
  all_hMPV_rt, all_adeno_rt, all_seasonal_CoV_229E_OC43_rt,
  all_seasonal_CoV_HKU1_NL63_rt
)
head(all_endemic_results)
unique(all_endemic_results$organism)
# keep <- unique(all_endemic_results$organism)

all_endemic_results %>%
  group_by(organism, date) %>%
  tally() %>%
  filter(n != 1)

all_endemic_results$organism <- as.factor(all_endemic_results$organism)
levels(all_endemic_results$organism)[c(5:8)]
levels(all_endemic_results$organism)[c(5:8)] <- c("IAV/H1N1", "IAV/H3N2", "IBV", "RV")

all_endemic_results_avg <- all_endemic_results %>%
  group_by(group, tag, level, organism) %>%
  mutate_at(
    c("lower", "upper", "median"),
    ~ zoo::rollmean(., k = 15, align = "center", fill = NA)
  ) %>%
  ungroup()

all_rt <- bind_rows(all_endemic_results, all_covid_rt)
all_rt
write.csv(all_rt, file = "Epidemia_Models/rt_all_pathogens_raw.csv", row.names = F)
write.csv(all_rt, file = "Block_Bootstrap_Mobility_vs_Rt/rt_all_pathogens_raw.csv", row.names = F)


all_rt_mv_avg <- bind_rows(all_endemic_results_avg, all_covid_rt_avg)
head(all_rt_mv_avg)
write.csv(all_rt_mv_avg, file = "Epidemia_Models/rt_all_pathogens_15day_mv_avg.csv", row.names = F)
write.csv(all_rt_mv_avg, file = "Block_Bootstrap_Mobility_vs_Rt/rt_all_pathogens_15day_mv_avg.csv", row.names = F)
