rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(readr) #readr_2.1.5
library(ggplot2) #ggplot2_3.5.1
library(lubridate) #lubridate_1.9.4

load("./section_4/wastewater_model_newzealand_censored_comp_fullAR.RData")
source("./functions_general/process_results_original.R")

## Get site weights 
sites <- read_csv("./section_4/data/sites.csv")
sites <- sites %>% 
  dplyr::select(SampleLocation, Region, Population,DisplayName) %>% 
  rename("site_id"= SampleLocation,
         "weights"=Population)

sites2 <- sites %>% 
  group_by(Region) %>% 
  mutate(region_weight = sum(weights)) %>% 
  mutate(region_weight = ifelse(Region == "Auckland", region_weight-weights[which(site_id == "AU_Mangere")], region_weight)) %>% 
  mutate(weights = ifelse(site_id == "AU_Mangere", 1.716*1000000- region_weight, weights))

# If your Rstudio doesn't have enough memory
# Step 1: Open terminal,
# 
# Step 2:
#   
#   cd ~
#   touch .Renviron
# open .Renviron
# Step 3: Save the following as the first line of .Renviron:
#   
#   R_MAX_VSIZE=100Gb 
# Step 4: Close RStudio and reopen

set.seed(2917)
ss <- sample(1:3000,2000)
samps1$samps <- samps1$samps[, ss]

results <- process_results(df_full =df_full,
                           tmbdat = tmbdat,
                           samps1 = samps1,
                           polyOrder=3,
                           id_group = 1,
                           weights_df = sites2 %>% dplyr::select(site_id, weights) )

save(file="./section_4/results_forplots_model_NZ.RData", list = "results")

##############################################
load(file="./section_4/results_forplots_model_NZ.RData")
results$station_ave_df  %>% arrange(-desc(sample_date)) %>% filter(ave_exp_v_u_fixed_deriv_lwr  >0) %>% dplyr::select(sample_date,ave_exp_v_u_fixed_deriv_lwr ) %>% slice(1)
results$df_full %>% group_by(sample_date) %>% slice(1) %>% ungroup()%>% filter(exp_v_deriv_lwr >0) %>% dplyr::select(sample_date, exp_v, exp_v_deriv_lwr) %>% slice(1)

gg1 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date, ave_exp_v_fixed_med)) +
  geom_line(size=0.2) + 
  # geom_line(data=results$df_full, aes(x= sample_date, exp_v), col = "red")+
  geom_ribbon(aes(ymax = ave_exp_v_fixed_upr, ymin = ave_exp_v_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),
                     limits = c(0,150000))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())


gg2 <- results$df_full %>% 
  group_by(sample_date) %>%
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(x = sample_date,y)) + 
  geom_line(aes(y = inst_repro),size=0.2) + 
  geom_ribbon(aes(ymax = inst_repro_upr, ymin = inst_repro_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste("exp[", bar(mu),"'(t)/",bar(mu),"(t)]")), breaks = scales::pretty_breaks(n=5), limits = c(0.8,1.5))+ 
  # geom_hline(yintercept=0, lty="dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())


gg3 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date, ave_exp_v_u_fixed_med)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_v_u_fixed_upr, ymin = ave_exp_v_u_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5),
                     limits = c(0,150000))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())


gg4 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date,inst_repro_med)) + 
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = inst_repro_upr, ymin = inst_repro_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste("exp[", bar(mu),"'(t)/",bar(mu),"(t)]")), breaks = scales::pretty_breaks(n=5),
                     limits = c(0.8,1.5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())


fest1 = cowplot::plot_grid(add_sub(gg1,expression("a) Signal, " * U[i](t) * " omitted"),size = 11),
                           add_sub(gg3,expression("b) Signal, " * U[i](t) * " included"),size = 11),
                           add_sub(gg2,expression("c) Geometric derivative of signal, " * U[i](t) * " omitted"),size = 11),
                           add_sub(gg4,expression("d) Geometric derivative of signal, " * U[i](t) * " included"),size = 11),
                           ncol=2, align="v", byrow = TRUE) 



ggsave(filename = paste0("./section_4/plots/allsignals_NZ_wastewatermodel.pdf"),
       plot = grid.arrange(fest1), 
       device = "pdf",
       width = 8, 
       height = 16/3,
       dpi = 300)

rstudioapi::viewer(paste0("./section_4/plots/allsignals_NZ_wastewatermodel.pdf"))

a = gg1+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_4/plots/allsignals_NZ_wastewatermodel_a.pdf"),
       plot = a, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)
rstudioapi::viewer(paste0("./section_4/plots/allsignals_NZ_wastewatermodel_a.pdf"))

b = gg3+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_4/plots/allsignals_NZ_wastewatermodel_b.pdf"),
       plot = b, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

c = gg2+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_4/plots/allsignals_NZ_wastewatermodel_c.pdf"),
       plot = c, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

d = gg4+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_4/plots/allsignals_NZ_wastewatermodel_d.pdf"),
       plot = d, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)


#### Do the plot with common + station-specific
to_plot <- sites2 %>% arrange(desc(weights)) %>% 
  ungroup() %>% 
  slice(1:15) 

to_plot$DisplayName[4] <- "Mangere"
  
dd = scales::pretty_breaks(n=10)
round_custom <- function(x) {
  round(x, -floor(log10(x)) * (x >= 100)) + (x < 100) * (round(x, -1) - x)
}

gg1 <- results$df_full %>% 
  filter(site_id %in% to_plot$site_id) %>% 
  left_join(to_plot, by = "site_id") %>% 
  mutate(DisplayName = factor(DisplayName, levels = to_plot$DisplayName)) %>% 
  ggplot(aes(x = sample_date, y)) +
  facet_wrap(~DisplayName, nrow = 3)+
  geom_line(aes(y = exp_v_u_fixed),size=0.2) + 
  geom_ribbon(aes(ymax = exp_v_u_fixed_upr, ymin = exp_v_u_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  geom_ribbon(aes(ymax =exp_v_fixed_upr, ymin = exp_v_fixed_lwr), alpha = 0.4, size = 0.2, fill = "red") +
  theme_bw()+
  geom_point(alpha =0.5, shape = 16, size = 0.2)+
  geom_point(data=results$df_full %>%  
               filter(site_id %in% to_plot$site_id) %>%  
               left_join(to_plot, by = "site_id") %>% 
               mutate(DisplayName = factor(DisplayName, levels = to_plot$DisplayName)) %>% 
               filter(censored_y == TRUE),col = "blue",alpha =0.5, shape = 16, size = 0.2)+
  # scale_y_continuous(name = expression(paste(mu[i],"(t)")), 
  #                    trans="log",
  #                    labels = function(x)format(x,digits=2),
  #                    breaks = c(0.14,1,8,50,400,3000))+ 
  scale_y_continuous(name = '', 
                     trans="log",
                     labels = function(x)format(x,digits=2,big.mark = ","),
                     breaks = exp(dd(log(c(results$df_full$exp_v_u_fixed_lwr,results$df_full$exp_v_u_fixed_upr)),8)))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               date_labels ="%b",
               name = "",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                     years <- year(x)
                                     years[duplicated(years)] <- ""  # Remove duplicate year labels
                                     years}))+
  theme(axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_text(vjust = -163),
        axis.title.y = element_blank())



ggsave(filename = "./section_4/plots/time_trend_fixed_AR_data.pdf",
       plot = gg1, 
       device = "pdf",
       width = 8.5, 
       height = 5,
       dpi = 300)

rstudioapi::viewer(paste0("./section_4/plots/time_trend_fixed_AR_data.pdf"))


