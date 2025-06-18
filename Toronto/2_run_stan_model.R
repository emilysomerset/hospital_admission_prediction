library(rstan) # rstan_2.32.7 
library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(forcats)
library(bayesplot)
library(ggplot2)

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


# sm <- stan_model("./Toronto/Stan models/model2_justcases.stan")
# 
# set.seed(2917)
# ss <- sample(1:3000,5)
# 
# mdl_fits <- lapply(as.list(1:5),function(i){
# ratio = data_foranalysis_full$analysis_d[[ss[i]]] %>% 
#   arrange(-desc(earliest_week_end_date)) %>% 
#   filter(!is.na(ratio_cases)) %$% ratio_v_u_fixed
# ######################################### Prepare data
# stan_data <- list(
#   J = nrow(case_age_matrix2),
#   I = ncol(case_age_matrix2[,-1]),
#   # Infection model  
#   Y = as.matrix(case_age_matrix2[,-1]) %>% t(),
#   ratio = ratio,
#   mean_X0 = init*2,
#   sd_X0 = c(200,50,20,3), 
#   phi_p = -log(0.5) / 1
# ) 
# 
# 
# 
# ### Fit model
# fit <- sampling(sm, data = stan_data,
#                 chains = 4, iter = 2000, warmup = 1000, seed = 123)
# 
# fit
# })
# 
# summary(fit)$summary %>%
#   as.data.frame() %>%
#   arrange(desc(Rhat))
# 
# mcmc_trace(as.array(fit), pars = "C_raw[1]")
# 
# 
# # Extract samples
# X_samples <- rstan::extract(fit, pars = "X")$X  # array of dim [iterations, I, J]
# X_tot_samples <- rstan::extract(fit, pars = "X_total")$X  # array of dim [iterations, I, J]
# 
# # Get mean across iterations
# X_mean <- apply(X_samples, c(2, 3), mean)  # matrix of [I, J] means
# X_tot_mean = apply(X_tot_samples, 2, median)  # matrix of [I, J] means
# 
# # Get quantiles
# X_lwr <- apply(X_samples, c(2, 3), quantile, 0.025)
# X_upr <- apply(X_samples, c(2, 3), quantile, 0.975)
# 
# X_tot_lwr <- apply(X_tot_samples, 2, quantile, 0.025)
# X_tot_upr <- apply(X_tot_samples, 2, quantile, 0.975)
# 
# 
# X_mean[1,] %>% plot(type = "l")
# X_mean[2,] %>% plot(type = "l")
# X_mean[3,] %>% plot(type = "l")
# X_mean[4,] %>% plot(type = "l")
# 
# mdlfit <- case_merged_age %>% 
#   filter(age != "Unknown")
# mdlfit$X_mean = c(X_mean)
# mdlfit$X_lwr = c(X_lwr)
# mdlfit$X_upr = c(X_upr)
# 
# mdlfit$X_tot_mean = rep(c(X_tot_mean),each =4)
# mdlfit$X_tot_lwr = rep(c(X_tot_lwr),each =4)
# mdlfit$X_tot_upr = rep(c(X_tot_upr),each =4)
# 
# ggplot(mdlfit, aes(earliest_week_end_date, X_mean, col = age))+ 
#   geom_ribbon(aes(ymin = X_lwr, ymax=X_upr, fill = age), alpha = 0.3)+
#   geom_line()+
#   scale_y_continuous(trans = "log10")
# 
# mdlfit %>% 
#   group_by(earliest_week_end_date) %>% 
#   mutate(X_mean = sum(X_mean)) %>% 
#   slice(1) %>% 
#   ggplot(aes(earliest_week_end_date, X_mean))+ 
#   geom_line()
# 
# mdlfit %>% 
#   group_by(earliest_week_end_date) %>% 
#   ggplot(aes(earliest_week_end_date, X_tot_mean))+ 
#   geom_ribbon(aes(ymin = X_tot_lwr, ymax=X_tot_upr), alpha = 0.3)+
#   geom_line()


#############################################################


stan_data <- list(
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


sm <- stan_model("./Toronto/Stan models/model4_cases_admissions_reparam_work.stan")

# X_val <- stan_data$Y * 2
# A_val <- pmax(round(stan_data$H / 10), 1)
# 
# init_fun <- function() {
#   list(
#     X = stan_data$Y*2,  # X > A, some positive values
#     A = round(stan_data$H/10),                # small positive numbers for admissions
#     logit_pi = rep(-1, stan_data$J),                                              # around 0 -> pi ~ 0.5
#     psi_p = 0,
#     C_raw = rep(0.8, length.out = (stan_data$I * (stan_data$I + 1) / 2)),
#     tau_logit = matrix(-1, nrow = stan_data$I, ncol = stan_data$J),
#     psi_tau = 0,
#     rho = 0,
#     theta = 0.5
#   )
# }

# I = stan_data$I
# J= stan_data$J
# init_fun <- function() {
#   list(
#     X_raw = matrix(rnorm(I * J, mean = 0.2, sd = 0.2), I, J),
#     logit_pi = rnorm(J, 0, 0.5),
#     psi_p = 0.5,
#     C_raw = abs(rnorm(I * (I + 1) / 2, mean = 1.5, sd = 0.5)),
#     
#     tau_logit = matrix(rnorm(I * J, -2, 0.5), I, J),
#     psi_tau = 0.5,
#     rho = 0
#   )
# }
# 
# init_fun <- function() list(
#   tau_logit = matrix(-2.3, stan_data$I, stan_data$J),   # gives p1 ≈ 0.05
#   logit_pi = rep(-0.105, stan_data$J)                   # gives p2 ≈ 0.45
# )


fit <- sampling(sm, data = stan_data,
                chains = 4, iter = 10000, warmup = 5000, cores= 4)




sm <- stan_model("./Toronto/Stan models/model4_cases_admissions_nonnormalapprox.stan")

init_fun <- function() {
  list(
    X = pmax(matrix(rpois(stan_data$I * stan_data$J, lambda = 2 * max(stan_data$Y)), stan_data$I, stan_data$J), stan_data$Y + 1),
    logit_pi_init = 0,
    z_logit_pi = rep(0, stan_data$J - 1),
    psi_p = 0.1,
    C_raw = rep(1, stan_data$I * (stan_data$I + 1) / 2),
    tau_logit = matrix(-2, stan_data$I, stan_data$J),
    psi_tau = 0.1
  )
}

fit <- sampling(sm, data = stan_data,
                chains = 4, iter = 10000, warmup = 5000, cores= 4,
                init = init_fun)



sm <- stan_model("./Toronto/Stan models/model4_cases_hosp_nonnormalapprox.stan")

stan_data <- list(
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
  H = as.matrix(hosp_age_matrix2[,-1]) %>% t(),
  phi_tau = -log(0.5)/1,
  K=3,
  A_tot = case_merged_age %>% group_by(earliest_week_end_date) %>% slice(1) %$% hospitalized_cases
) 


init_fun <- function() {
  list(
    X = pmax(matrix(rpois(stan_data$I * stan_data$J, lambda = 2 * max(stan_data$Y)), stan_data$I, stan_data$J), stan_data$Y + 1),
    logit_pi_init = 0,
    z_logit_pi = rep(0, stan_data$J - 1),
    psi_p = 0.1,
    C_raw = rep(1, stan_data$I * (stan_data$I + 1) / 2),
    tau_logit = matrix(-2, stan_data$I, stan_data$J),
    psi_tau = 0.1,
    A = as.matrix(adm_age_matrix2[,-1]) %>% t()
  )
}

fit <- sampling(sm, data = stan_data,
                chains = 4, iter = 10000, warmup = 5000, cores= 4, init=init_fun)

summary(fit)$summary %>%
  as.data.frame() %>%
  arrange(desc(Rhat)) 


summary(fit)$summary %>%
  as.data.frame() %>% 
  filter(grepl("theta", rownames(.)))

mcmc_trace(as.array(fit), pars = c("X[1,44]","X[2,44]"))+ 
  ggplot2::facet_wrap(~parameter, scales = "free_y")
mcmc_trace(as.array(fit), pars = "p[2,29,1]")
mcmc_trace(as.array(fit), pars = "X_total[48]")
stan_data$Y[2,29]

mcmc_scatter(as.array(fit),pars = c("p[1,44,1]","A[1,44]"), transformations = "log")

mcmc_trace(as.array(fit), pars = "tau[2,30]")

(stan_data$H/7)/(stan_data$Y)

mcmc_trace(as.array(fit), pars = "tau[3,132]")

#20%

mcmc_trace(as.array(fit), pars = "tau[1,132]")

X


# Extract samples
psi_p_samples <- rstan::extract(fit, pars = "psi_p")$psi_p  # array of dim [iterations, I, J]
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


