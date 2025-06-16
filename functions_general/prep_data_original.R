
prep_data <- function(outcome_column_name,
                      denom_column_name,
                      site_id, 
                      sample_date,
                      data,
                      polyOrder,
                      pred_also,
                      nknots = NULL){

# library(OSplines) # OSplines_0.1.1
library(aghq) # aghq_0.4.1
library(readr) # readr_2.1.5
library(dplyr) # dplyr_1.1.4
library(magrittr) # magrittr_2.0.3 
library(lubridate) # lubridate_1.9.3 
library(Matrix) # Matrix_1.7-0 
library(TMB) #TMB_1.9.14 

# Needed function
compute_weights_precision <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}

df <- data %>% 
  rename(y = outcome_column_name, 
         denom = denom_column_name, 
         sample_date = sample_date, 
         site_id = site_id) %>% 
  dplyr::select(y, denom, sample_date, site_id, censored_y)

if (pred_also == TRUE){
  rr = df %$% sample_date %>% range()}

if (pred_also == FALSE){
  rr = df %>% filter(!is.na(y)) %$% sample_date %>% range()
}

df <- expand.grid(sample_date = seq(rr[1],rr[2],1),
                  site_id = unique(df$site_id)) %>% 
  left_join(df, by = c("sample_date","site_id")) %>% 
  dplyr::select(sample_date, site_id,  y, denom, censored_y) 

df =df %>%  
  arrange(site_id,sample_date) %>% 
  ungroup() 

df_full <- df %>% 
  mutate(obs = !is.na(y)) %>% 
  mutate(nindex = 1:nrow(.))

df <- df %>% 
  filter(!is.na(y)) %>% 
  mutate(nindex2 = 1:nrow(.))

df$t <- df$sample_date %>% as.numeric()
df$t <- df$t - min(df$t)
df_full$t <- df_full$sample_date %>% as.numeric()
df_full$t <- df_full$t - min(df_full$t)

stationsizes = df_full %>%  
  group_by(site_id) %>%  
  summarise(ll = length(site_id)) %$% ll

n1 = df_full %>%  
  mutate(n1 = 0:(nrow(.)-1)) %>%  
  group_by(site_id) %>%  
  slice(1) %$% n1 

n1 <- c(n1, (nrow(df_full)))

# Setup the design and precision matrix 
polyOrder = polyOrder
if (is.null(nknots)){knots <- seq(0,(max(df$t)), length.out = floor((max(df$t))/7))}else{
  knots <- seq(0,(max(df$t)), length.out = nknots)
}

B <- local_poly_helper(knots = knots, 
                       refined_x = df$t,
                       p = polyOrder)

P <- compute_weights_precision(x=knots)

X_global = global_poly_helper(x=df$t, p= polyOrder)
Xfstat <- model.matrix(~site_id, 
                       data = df %>% mutate(site_id = factor(site_id)), 
                       contrasts.arg = list(site_id = "contr.sum"))[,-1]
Xfstat <- as.matrix(Xfstat)
X = cbind(X_global, Xfstat)

daily <- model.matrix(~t -1, 
                      data = df_full %>% mutate(t = factor(t)))


obs <- model.matrix(~nindex -1, data = df_full %>% mutate(nindex=factor(nindex)) %>% filter(obs))
station <- model.matrix(~factor(site_id) -1, data = df)
y_ind_obs <- df %>% filter(censored_y==FALSE)%$%nindex2 -1
y_ind_cens <- df %>% filter(censored_y==TRUE)%$%nindex2 - 1
observed_y <- model.matrix(~nindex2 -1, data = df %>% mutate(nindex2=factor(nindex2)) %>% filter(censored_y==FALSE))
censored_y <- model.matrix(~nindex2 -1, data = df %>% mutate(nindex2=factor(nindex2)) %>% filter(censored_y==TRUE))

tmbdat <- list(
  # Design matrix
  Xfstat = as(Xfstat, 'dgTMatrix'),
  X_global = as(X_global, 'dgTMatrix'),
  X = as(X, 'dgTMatrix'),
  B = as(B, 'dgTMatrix'),
  P = as(P,'dgTMatrix'),
  daily = as(daily, 'dgTMatrix'),
  obs = as(obs, 'dgTMatrix'),
  logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
  # Response
  y = df$y,
  # PC Prior params
  n1=n1,
  stationsizes = stationsizes,
  denom = df$denom,
  y_ind_obs=y_ind_obs,
  y_ind_cens = y_ind_cens,
  cens_dir = rep(1, length(y_ind_cens)),
  station = as(station, 'dgTMatrix'),
  knots = knots
)

return(list(tmbdat=tmbdat,df_full = df_full))

  }
