process_results <- function(data_foranalysis, 
                            MI_models,
                            adj = TRUE, 
                            mean_adj = NULL, 
                            sd_adj = NULL,
                            full_timeseries = FALSE){

ncounts = length(data_foranalysis$tmbdat[[1]]$case_counts)
lag = 0

set.seed(29172)
ss <- 1:500
len = ncounts+1 - lag 
logit_p_samps1 <- lapply(MI_models, function(samps){samps$samps[1:(len),ss ]})
Z_samps1 <- lapply(MI_models, function(samps){samps$samps[(len  + 1):((len  +len)),ss]})
thetasamps <- lapply(MI_models, function(samps){samps$thetasamples[[1]][ss]})

lensim = ncol(Z_samps1[[1]])*length(Z_samps1)

### combine
logit_p_samps1 <- Reduce(cbind, logit_p_samps1)
thetasamps <- Reduce(c, thetasamps)

p_samps = exp(logit_p_samps1)/(1+  exp(logit_p_samps1))

analysis_d <- data_foranalysis$analysis_d2 %>% 
  mutate(row_id = 1:nrow(.)) %>% 
  filter(row_id >= (index_start[1]-1)) 

if (full_timeseries ==TRUE){
  analysis_d <- analysis_d %>% 
    filter(!(is.na(ratio) & (row_id > (index_start[1]))))%>% 
    dplyr::select(-c("variable","value","exp_v","sum_v","shift_sum_v","ratio","wastewater_count_index"))
}else{  analysis_d <- analysis_d%>% 
  dplyr::select(-c("variable","value","exp_v","sum_v","shift_sum_v","ratio","wastewater_count_index")) }

analysis_d$p_med = as.numeric(apply(p_samps, MARGIN=1,median))
analysis_d$p_upr = as.numeric(apply(p_samps, MARGIN=1,quantile,0.975))
analysis_d$p_lwr = as.numeric(apply(p_samps, MARGIN=1,quantile, 0.025))


# gg1

#### Latent Z value
Z_samps1 <- Reduce(cbind, Z_samps1)
Z_samps1 <- abs(Z_samps1)
Z_samps1_cumsum <- apply(Z_samps1, 2, cumsum)
Z_samps1_cumsum_delta <- apply(Z_samps1_cumsum, 2, function(x){x-x[1]})

analysis_d$z_med = as.numeric(apply(Z_samps1, MARGIN=1,median))
analysis_d$z_upr = as.numeric(apply(Z_samps1, MARGIN=1,quantile,0.975))
analysis_d$z_lwr = as.numeric(apply(Z_samps1, MARGIN=1,quantile,0.025))

set.seed(2)
y_samps = rnorm(len*lensim, c(p_samps)*c(Z_samps1),sqrt(c(p_samps)*(1-c(p_samps))*c(Z_samps1)))
y_samps <- matrix(y_samps, ncol = lensim, byrow=FALSE)

analysis_d$y_med = as.numeric(apply(y_samps, MARGIN=1,median))
analysis_d$y_upr = as.numeric(apply(y_samps, MARGIN=1,quantile,0.975))
analysis_d$y_lwr = as.numeric(apply(y_samps, MARGIN=1,quantile,0.025))

analysis_d$z_cumsum_noadj_med = as.numeric(apply(Z_samps1_cumsum, MARGIN=1,median)) 
analysis_d$z_cumsum_noadj_upr = as.numeric(apply(Z_samps1_cumsum, MARGIN=1,quantile,0.975)) 
analysis_d$z_cumsum_noadj_lwr = as.numeric(apply(Z_samps1_cumsum, MARGIN=1,quantile,0.025)) 

analysis_d$z_cumsum_noadj_delta_med = as.numeric(apply(Z_samps1_cumsum_delta, MARGIN=1,median)) 
analysis_d$z_cumsum_noadj_delta_upr = as.numeric(apply(Z_samps1_cumsum_delta, MARGIN=1,quantile,0.975)) 
analysis_d$z_cumsum_noadj_delta_lwr = as.numeric(apply(Z_samps1_cumsum_delta, MARGIN=1,quantile,0.025)) 

if (adj == TRUE){
phase = rnorm(lensim, mean = mean_adj, sd = sd_adj)
analysis_d$z_cumsum_med = as.numeric(apply((Z_samps1_cumsum + phase), MARGIN=1,median)) 
analysis_d$z_cumsum_upr = as.numeric(apply((Z_samps1_cumsum + phase), MARGIN=1,quantile,0.975)) 
analysis_d$z_cumsum_lwr = as.numeric(apply((Z_samps1_cumsum + phase), MARGIN=1,quantile,0.025)) }

analysis_d <- analysis_d %>% 
  mutate(number_of_cases_cumsum = (cumsum(y)+ analysis_d$Y0[1]),
         number_of_cases_cumsum_delta = number_of_cases_cumsum- number_of_cases_cumsum[1])

return(analysis_d)

}