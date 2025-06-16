
prep_data <- function(results,
                      case_data,
                      y_var = "number_of_cases",
                      date_forfilter_strict_upr = NULL,
                      date_forfilter_strict_lwr = NULL, 
                      AR = FALSE,
                      weight_ratio = FALSE, 
                      weights_datadic = NULL){

# library(OSplines) # OSplines_0.1.1
library(aghq) # aghq_0.4.1
library(readr) # readr_2.1.5
library(dplyr) # dplyr_1.1.4
library(magrittr) # magrittr_2.0.3 
library(lubridate) # lubridate_1.9.3 
library(Matrix) # Matrix_1.7-0 
library(TMB) #TMB_1.9.14 

tmp <- results$v

tmp <- tmp %>% 
  mutate(exp_v = exp(value))

## 
tmp <-tmp %>% 
  group_by(variable) %>% 
  mutate(sum_v = c(rep(NA,6), zoo::rollapply(exp_v,7,sum))) %>% 
  # filter(!is.na(sum_v)) %>% 
  mutate(shift_sum_v = c(rep(NA,7), sum_v[1:(length(value)-7)])) %>% 
  mutate(ratio = sum_v/shift_sum_v )

tmp <- tmp %>% 
  filter(!is.na(ratio))

if (AR == TRUE){
  
tmp_ar <- lapply(list(results$v_u, results$v_u_fixed,results$v_fixed), 
       function(x){
         tmp2 <- x %>%
           mutate(value = exp(value)) %>%
           left_join(weights_datadic, by = "site_id") %>% 
           mutate(value = value*weight) %>%
           group_by(variable, sample_date) %>%
           summarise(value = sum(value)) %>%
           mutate(value = log(value))
         
         tmp2 <- tmp2 %>% 
           mutate(exp_value = exp(value))
         
         ## cases
         tmp2 <-tmp2 %>% 
           group_by(variable) %>% 
           mutate(sum_value = c(rep(NA,6), zoo::rollapply(exp_value,7,sum))) %>% 
           mutate(shift_sum_value = c(rep(NA,7), sum_value[1:(length(value)-7)])) %>% 
           mutate(ratio_value = sum_value/shift_sum_value )
         
         tmp2 <- tmp2 %>% 
           filter(!is.na(ratio_value))
       })
}
 

# if (weight_ratio == TRUE){
#   tmp3 <- split(results$v_u, results$v_u$site_id)
#   
#   tmp3 <- lapply(tmp3, function(x){x %>% mutate(exp_v_u = exp(value))})
# 
#   
#   tmp3 <- lapply(tmp3, function(x){
#    x%>% 
#     group_by(variable) %>% 
#     mutate(sum_v_u = c(rep(NA,6), zoo::rollapply(exp_v_u,7,sum))) %>% 
#     # filter(!is.na(sum_v)) %>% 
#     mutate(shift_sum_v_u = c(rep(NA,7), sum_v_u[1:(length(value)-7)])) %>% 
#     mutate(ratio_v_u = sum_v_u/shift_sum_v_u ) %>% 
#     filter(!is.na(ratio_v_u))
#   })
# 
#   tmp3 <- Reduce(rbind, tmp3) %>% 
#     left_join(weights_datadic, by = "site_id") %>% 
#     mutate(ratio_v_u = weight*ratio_v_u) %>% 
#     group_by(sample_date, variable) %>% 
#     summarise(ratio_v_u_weighted = sum(ratio_v_u))
# }

case_data <- case_data %>% 
  rename(y = y_var) %>% 
  arrange(-desc(earliest_week_end_date)) %>% 
  mutate(ratio_cases = c(NA, y[2:nrow(.)]/y[1:(nrow(.)-1)]))

if (!is.null(date_forfilter_strict_upr)){case_data = case_data %>% filter(earliest_week_end_date < date_forfilter_strict_upr)}
if (!is.null(date_forfilter_strict_lwr)){case_data = case_data %>% filter(earliest_week_end_date > date_forfilter_strict_lwr)}

analysis_d <- case_data %>% 
  left_join(tmp %>% 
              mutate(wastewater_count_index = 0:(length(value)-1)) %>% 
              mutate(sample_date = sample_date), 
            by = c( "earliest_week_end_date"="sample_date")) %>% 
  filter(!is.na(ratio))

if (AR == TRUE){
  analysis_d <- analysis_d %>% 
    left_join(tmp_ar[[1]] %>% dplyr::select(sample_date, ratio_value, variable) %>% 
                rename(ratio_v_u = ratio_value), 
              by = c( "earliest_week_end_date"="sample_date","variable"))
  
  analysis_d <- analysis_d %>% 
    left_join(tmp_ar[[2]] %>% dplyr::select(sample_date, ratio_value, variable) %>% 
                rename(ratio_v_u_fixed = ratio_value), 
              by = c( "earliest_week_end_date"="sample_date","variable"))
  
  analysis_d <- analysis_d %>% 
    left_join(tmp_ar[[3]] %>% dplyr::select(sample_date, ratio_value, variable) %>% 
                rename(ratio_v_fixed = ratio_value), 
              by = c( "earliest_week_end_date"="sample_date","variable"))
}

# if (weight_ratio == TRUE){
#   analysis_d <- analysis_d %>% 
#     left_join(tmp3 %>% dplyr::select(sample_date, ratio_v_u_weighted, variable), 
#               by = c( "earliest_week_end_date"="sample_date","variable"))
# }

analysis_d2 <- case_data %>% 
  left_join(tmp %>% 
              mutate(wastewater_count_index = 0:(length(value)-1)) %>% 
              mutate(sample_date = sample_date), 
            by = c( "earliest_week_end_date"="sample_date")) %>% 
  group_by(earliest_week_end_date) %>% 
  slice(1) %>% 
  ungroup() 

if (AR == TRUE){
  analysis_d2 <- analysis_d2 %>% 
    left_join(tmp_ar[[1]] %>% dplyr::select(sample_date, ratio_value, variable) %>% 
                rename(ratio_v_u = ratio_value), 
              by = c( "earliest_week_end_date"="sample_date","variable"))
  
  analysis_d2 <- analysis_d2 %>% 
    left_join(tmp_ar[[2]] %>% dplyr::select(sample_date, ratio_value, variable) %>% 
                rename(ratio_v_u_fixed = ratio_value), 
              by = c( "earliest_week_end_date"="sample_date","variable"))
  
  analysis_d2 <- analysis_d2 %>% 
    left_join(tmp_ar[[3]] %>% dplyr::select(sample_date, ratio_value, variable) %>% 
                rename(ratio_v_fixed = ratio_value), 
              by = c( "earliest_week_end_date"="sample_date","variable"))
}

# if (weight_ratio == TRUE){
#   analysis_d2 <- analysis_d2 %>% 
#     left_join(tmp3 %>% dplyr::select(sample_date, ratio_v_u_weighted, variable), 
#               by = c( "earliest_week_end_date"="sample_date","variable"))
# }

analysis_d <- split(analysis_d, analysis_d$variable)

tmbdat <- lapply(analysis_d, function(dat){
  
  list(
    case_counts = dat$y,
    ratio = dat$ratio,
    ratio_v_u = dat$ratio_v_u, #will sometimes be null
    ratio_v_u_fixed = dat$ratio_v_u_fixed, #will sometimes be null
    ratio_v_fixed = dat$ratio_v_fixed #will sometimes be null
  )
  
  
})

return(list(tmbdat=tmbdat,analysis_d = analysis_d,analysis_d2 = analysis_d2))

}
