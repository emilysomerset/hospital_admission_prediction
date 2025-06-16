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
library(ggplot2)
library(lubridate)

load("./section_3/wastewater_model_ON_phachist.RData")
source("./functions_general/process_results_original.R")

weights_datadic = data.frame(site_id = c("TAB","THC","THU","TNT"), 
                             weights = c(0.499/(0.177+0.232+0.064+0.499),
                                        0.177/(0.177+0.232+0.064+0.499),
                                        0.232/(0.177+0.232+0.064+0.499),
                                        0.064/(0.177+0.232+0.064+0.499)))

results <- process_results(df_full =df_full,
                           tmbdat = tmbdat,
                           samps1 = samps1,
                           polyOrder=3,
                           id_group = 1,
                           weights_df = weights_datadic)


##############################################
gg1 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date, ave_exp_v_fixed_med)) +
  geom_line(size=0.2) + 
  # geom_line(data=results$df_full, aes(x= sample_date, exp_v), col = "red")+
  geom_ribbon(aes(ymax = ave_exp_v_fixed_upr, ymin = ave_exp_v_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5), limits = c(0,565))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())


gg2 <- results$df_full %>% 
  group_by(sample_date) %>%
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(x = sample_date,y)) + 
  geom_line(aes(y = inst_repro),size=0.2) + 
  geom_ribbon(aes(ymax = inst_repro_upr, ymin = inst_repro_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste("exp[", bar(mu),"'(t)/",bar(mu),"(t)]")), breaks = scales::pretty_breaks(n=5), limits = c(0.7,1.3))+ 
  # geom_hline(yintercept=0, lty="dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())


gg3 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date, ave_exp_v_u_fixed_med)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_v_u_fixed_upr, ymin = ave_exp_v_u_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5), limits = c(0,565))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())

gg4 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date,inst_repro_med)) + 
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = inst_repro_upr, ymin = inst_repro_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste("exp[", bar(mu),"'(t)/",bar(mu),"(t)]")), breaks = scales::pretty_breaks(n=5),
                     limits = c(0.7,1.3))+ 
  # geom_hline(yintercept=0, lty="dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x.top = element_text(vjust = -65),
        axis.ticks.x.top = element_blank())


fest1 = cowplot::plot_grid(add_sub(gg1,expression("a) Signal, " * U[i](t) * " omitted"),size = 11),
                           add_sub(gg3,expression("b) Signal, " * U[i](t) * " included"),size = 11),
                           add_sub(gg2,expression("c) Geometric derivative of signal, " * U[i](t) * " omitted"),size = 11),
                           add_sub(gg4,expression("d) Geometric derivative of signal, " * U[i](t) * " included"),size = 11),
                           ncol=2, align="v", byrow = TRUE) 



ggsave(filename = paste0("./section_3/plots/allsignals_ON_wastewatermodel.pdf"),
       plot = grid.arrange(fest1), 
       device = "pdf",
       width = 8, 
       height = 16/3,
       dpi = 300)

rstudioapi::viewer(paste0("./section_3/plots/allsignals_ON_wastewatermodel.pdf"))



a = gg1+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_ON_wastewatermodel_a.pdf"),
       plot = a, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)
rstudioapi::viewer(paste0("./section_3/plots/allsignals_ON_wastewatermodel_a.pdf"))

b = gg3+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_ON_wastewatermodel_b.pdf"),
       plot = b, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

c = gg2+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_ON_wastewatermodel_c.pdf"),
       plot = c, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)

d = gg4+ theme(axis.title.y = element_blank(),
               axis.text.x.top = element_text(vjust = -78))
ggsave(filename = paste0("./section_3/plots/allsignals_ON_wastewatermodel_d.pdf"),
       plot = d, 
       device = "pdf",
       width = 8/2, 
       height = 8/3,
       dpi = 300)




proper_df = data.frame(site_id = c("TAB","THC","THU","TNT"),
                       site_name = c("Ashbridges Bay","Highland Creek",
                                     "Humber","North Toronto"))

gg5 <- results$df_full %>% 
  left_join(proper_df, by = "site_id") %>% 
  group_by(sample_date) %>% 
  ggplot(aes(x = sample_date, y)) +
  geom_line(aes(y = exp_v_u_fixed),size=0.2) + 
  geom_ribbon(aes(ymax = exp_v_u_fixed_upr, ymin = exp_v_u_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  geom_point(size = 0.05)+
  facet_wrap(~site_name, ncol = 1)+
  scale_y_continuous(name = expression(paste(mu[i],"(t)")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 16),
        axis.text.x.top = element_text(vjust = -235),
        axis.ticks.x.top = element_blank())

gg6 <- results$df_full %>% 
  left_join(proper_df, by = "site_id") %>% 
  group_by(sample_date) %>% 
  ggplot(aes(x = sample_date, y)) +
  # geom_point(size = 0.1)+
  geom_line(aes(y = exp_v_u_fixed_deriv),size=0.2) + 
  geom_ribbon(aes(ymax = exp_v_u_fixed_deriv_upr, ymin = exp_v_u_fixed_deriv_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  facet_wrap(~site_name, ncol = 1)+
  scale_y_continuous(name = expression(paste(mu[i],"(t)'")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 16),
        axis.text.x.top = element_text(vjust = -235),
        axis.ticks.x.top = element_blank())

ggsave(filename = paste0("./section_3/plots/indvsignals_ON_wastewatermodel.pdf"),
       plot = grid.arrange(gg5, gg6, ncol = 2), 
       device = "pdf",
       width = 8, 
       height = 6,
       dpi = 300)


#### Do the plot with common + station-specific
dd = scales::pretty_breaks(n=10)
round_custom <- function(x) {
  round(x, -floor(log10(x)) * (x >= 10)) + (x < 10) * (round(x, -1) - x)
}

gg1 <- results$df_full %>% 
  left_join(proper_df, by = "site_id") %>% 
  ggplot(aes(x = sample_date, y)) +
  facet_wrap(~site_name, nrow = 3)+
  geom_line(aes(y = exp_v_u_fixed),size=0.2) + 
  geom_ribbon(aes(ymax = exp_v_u_fixed_upr, ymin = exp_v_u_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  geom_ribbon(aes(ymax =exp_v_fixed_upr, ymin = exp_v_fixed_lwr), alpha = 0.4, size = 0.2, fill = "red") +
  theme_bw()+
  geom_point(alpha =0.5, shape = 16, size = 0.2)+
  # scale_y_continuous(name = expression(paste(mu[i],"(t)")),
  #                    # trans="log",
  #                    # labels = function(x)format(x,digits=2),
  #                    breaks = scales::pretty_breaks(n=8))+
  scale_y_continuous(name =expression(mu[i]^"*"*(t)),
                     trans="log",
                     labels = function(x)format(x,digits=2,big.mark = ","),
                     breaks = c(1,3,8,20,55,150,400,1100))+
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
        axis.text.x.top = element_text(vjust = -163))

ggsave(filename = "./section_3/plots/time_trend_fixed_AR_data.pdf",
       plot = gg1, 
       device = "pdf",
       width = 8.5, 
       height = 5,
       dpi = 300)

rstudioapi::viewer(paste0("./section_3/plots/time_trend_fixed_AR_data.pdf"))
