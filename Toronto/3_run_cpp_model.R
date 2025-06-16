library(rstan) # rstan_2.32.7 
library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(forcats)
library(bayesplot)
library(ggplot2)
library(aghq)
library(TMB)

load('./Toronto/work_d.RData')


### Get wastewater reproduction 
#############################################################
load(file="./Toronto/results_model_ON_phachist.RData")
source("./functions_general/prep_data_covid_with_fittedwastewater.R")

weights_datadic = data.frame(site_id = c("TAB","THC","THU","TNT"), 
                             weight = c(0.499/(0.177+0.232+0.064+0.499),
                                        0.177/(0.177+0.232+0.064+0.499),
                                        0.232/(0.177+0.232+0.064+0.499),
                                        0.064/(0.177+0.232+0.064+0.499)))

data_foranalysis_full <- prep_data(case_data = case_age %>% group_by(earliest_week_end_date) %>% summarise(totcase = mean(totcase)),
                                   y_var = "totcase",
                                   results = results,
                                   AR=TRUE,
                                   weight_ratio = TRUE,
                                   weights_datadic = weights_datadic)

case_age %>% group_by(earliest_week_end_date) %>% summarise(totcase = mean(totcase))

data_foranalysis_full$analysis_d[[1]] %>% 
  arrange(-desc(earliest_week_end_date)) %>% 
  filter(!is.na(ratio_cases)) %$% range(earliest_week_end_date)

ratio = data_foranalysis_full$analysis_d[[1]] %>% 
  arrange(-desc(earliest_week_end_date)) %>% 
  filter(!is.na(ratio_cases)) %$% ratio_v_u_fixed


###############

case_merged_age = case_age %>% 
  mutate(age = fct_collapse(age, `<39` = c("0 to 4","5 to 11", "12 to 19", "20 to 39"))) %>% 
  mutate(age = fct_collapse(age, `<59` = c("<39","40 to 59"))) %>% 
  mutate(age = fct_collapse(age, `60+` = c("60 to 79","80+"))) %>% 
  group_by(age, earliest_week_end_date) %>% 
  mutate(case = sum(case),
         integer_number = sum(integer_number)) %>% 
  select(-'average_bed_occupancy',-'age_orig_hosp') %>% 
  slice(1) %>% 
  arrange(-desc(earliest_week_end_date))

init = case_merged_age$case[1:2]

case_merged_age<-case_merged_age %>% 
filter(earliest_week_end_date >= "2021-10-23" & earliest_week_end_date <= "2024-04-27")

case_merged_age <- case_merged_age %>%
  filter(!is.na(hospitalized_cases))

case_age_matrix2 <- case_merged_age  %>% dcast(earliest_week_end_date~age, value.var = 'case') %>% 
  dplyr::select(-'Unknown')

hosp_age_matrix2 <- case_merged_age  %>% dcast(earliest_week_end_date~age, value.var = 'integer_number') %>% 
  dplyr::select(-'Unknown')

adm_age_matrix2 <- case_merged_age %>% 
  filter(age != "Unknown") %>% 
  mutate(adm_tilde = integer_number) %>% 
  group_by(earliest_week_end_date) %>% 
  mutate(adm_star = ceiling(adm_tilde*hospitalized_cases[1]/sum(adm_tilde))) %>% 
  dcast(earliest_week_end_date~age, value.var = 'adm_star') 

case_age_matrix2
hosp_age_matrix2
adm_age_matrix2

case_age_matrix2 %>% dim()
hosp_age_matrix2 %>% dim()

case_merged_age %>% 
  filter(age == "Unknown")%$% case %>% sum() 


tmbdat <- list(
  J = nrow(case_age_matrix2),
  I = ncol(case_age_matrix2[,-1]),

  # Infection model  
  Y = as.matrix(case_age_matrix2[,-1]) %>% t(),
  ratio = ratio[1:56],
  mean_X0 = init*2,
  # sd_X0 = c(200,50,20,3), 
  sd_X0 = c(100, 10),
  phi_p = -log(0.5) / 1,
  
  # Hospital 
  A = as.matrix(adm_age_matrix2[,-1]) %>% t(),
  phi_tau = -log(0.5)/1
) 


compile(file="./Toronto/cpp/model4_cases_admission_reparam_work.cpp")
try(dyn.unload(dynlib("./Toronto/cpp/model4_cases_admission_reparam_work")),silent = TRUE)
dyn.load(dynlib("./Toronto/cpp/model4_cases_admission_reparam_work"))

tmbparams <- list(
  X = matrix(10, nrow = tmbdat$I, ncol = tmbdat$J),                   # latent infections
  logit_pi_init = 0,                                    # logit(pi_1)
  z_logit_pi = rep(0, tmbdat$J - 1),                           # RW std normal increments
  psi_p = log(0.5),                                     # log sd of pi RW
  
  C_raw = rep(1.5, tmbdat$I * (tmbdat$I + 1) / 2),                    # symmetric contact matrix
  
  tau_logit = matrix(-2, nrow = tmbdat$I, ncol = tmbdat$J),           # logit(tau_ij)
  psi_tau = log(0.5),                                   # log sd of tau AR(1)
  rho = 0.5                                             # AR(1) correlation
)


ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = c("X", "z_logit_pi", "tau_logit"),
  DLL = "model4_cases_admission_reparam_work",
  silent = TRUE
)

aghq_k = 10

mdl1 <- aghq::marginal_laplace_tmb(ff,k=aghq_k,startingvalue = startingvalue <- c(
  0,             # logit_pi_init
  log(0.5),      # psi_p
  rep(1.5, 3),  # C_raw
  log(0.5),      # psi_tau
  0.5            # rho
))
samps1 <- aghq::sample_marginal(mdl1, M = 500) 






fit <- sampling(sm, data = stan_data,
                chains = 4, iter = 10000, warmup = 5000, cores= 4)




summary(fit)$summary %>%
  as.data.frame() %>%
  arrange(desc(Rhat)) 

mcmc_trace(as.array(fit), pars = "X_total[39]")
mcmc_trace(as.array(fit), pars = "C[2,2]")

mcmc_trace(as.array(fit), pars = "tau[2,30]")

(stan_data$H/7)/(stan_data$Y)

mcmc_trace(as.array(fit), pars = "tau[3,132]")

#20%

mcmc_trace(as.array(fit), pars = "tau[1,132]")

X


# Extract samples
X_samples <- rstan::extract(fit, pars = "X")$X  # array of dim [iterations, I, J]
p_samples <- rstan::extract(fit, pars = "p")$p  # array of dim [iterations, I, J]
X_tot_samples <- rstan::extract(fit, pars = "X_total")$X  # array of dim [iterations, I, J]

# Get mean across iterations
X_mean <- apply(X_samples, c(2, 3), mean)  # matrix of [I, J] means
X_tot_mean = apply(X_tot_samples, 2, median)  # matrix of [I, J] means
p_med <- apply(p_samples, c(2,3,4), median)

# Get quantiles
X_lwr <- apply(X_samples, c(2, 3), quantile, 0.025)
X_upr <- apply(X_samples, c(2, 3), quantile, 0.975)
p_upr <- apply(p_samples, c(2,3,4), quantile, 0.975)

X_tot_lwr <- apply(X_tot_samples, 2, quantile, 0.025)
X_tot_upr <- apply(X_tot_samples, 2, quantile, 0.025)
p_lwr <- apply(p_samples, c(2,3,4), quantile, 0.025)


X_mean[1,] %>% plot(type = "l")
X_mean[2,] %>% plot(type = "l")
X_mean[3,] %>% plot(type = "l")
X_mean[4,] %>% plot(type = "l")

mdlfit <- case_merged_age %>% 
  filter(age != "Unknown")
mdlfit$X_mean = c(X_mean)
mdlfit$X_lwr = c(X_lwr)
mdlfit$X_upr = c(X_upr)

mdlfit$prob_admission_med <- c(p_med[,,1])
mdlfit$prob_admission_lwr <- c(p_lwr[,,1])
mdlfit$prob_admission_upr <- c(p_upr[,,1])

mdlfit$prob_rep_med <- c(p_med[,,2])
mdlfit$prob_rep_lwr <- c(p_lwr[,,2])
mdlfit$prob_rep_upr <- c(p_upr[,,2])

mdlfit$prob_notrep_med <- c(p_med[,,3])
mdlfit$prob_notrep_lwr <- c(p_lwr[,,3])
mdlfit$prob_notrep_upr <- c(p_upr[,,3])

mdlfit$X_tot_mean = rep(c(X_tot_mean),each =4)
mdlfit$X_tot_lwr = rep(c(X_tot_lwr),each =4)
mdlfit$X_tot_upr = rep(c(X_tot_upr),each =4)

ggplot(mdlfit, aes(earliest_week_end_date, X_mean, col = age))+ 
  geom_ribbon(aes(ymin = X_lwr, ymax=X_upr, fill = age), alpha = 0.3)+
  geom_line()
mdlfit %>% 
  group_by(earliest_week_end_date) %>% 
  mutate(X_mean = sum(X_mean)) %>% 
  slice(1) %>% 
  ggplot(aes(earliest_week_end_date, X_mean))+ 
  geom_line()

mdlfit %>% 
  group_by(earliest_week_end_date) %>% 
  ggplot(aes(earliest_week_end_date, X_tot_mean))+ 
  geom_ribbon(aes(ymin = X_tot_lwr, ymax=X_tot_upr), alpha = 0.3)+
  geom_line()

mdlfit %>% 
  ggplot(aes(earliest_week_end_date, prob_admission_med, col = age))+ 
  geom_ribbon(aes(ymin = prob_admission_lwr, ymax=prob_admission_upr, fill = age), alpha = 0.3)+ 
  geom_line()

mdlfit %>% 
  ggplot(aes(earliest_week_end_date, prob_rep_med, col = age))+ 
  geom_ribbon(aes(ymin = prob_rep_lwr, ymax=prob_rep_upr, fill = age), alpha = 0.3)+ 
  geom_line()

mdlfit %>% 
  ggplot(aes(earliest_week_end_date, prob_notrep_med, col = age))+ 
  geom_ribbon(aes(ymin = prob_notrep_lwr, ymax=prob_notrep_upr, fill = age), alpha = 0.3)+ 
  geom_line()


