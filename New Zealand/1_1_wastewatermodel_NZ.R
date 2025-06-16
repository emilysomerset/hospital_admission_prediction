rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(readr) #readr_2.1.5
library(janitor) #janitor_2.2.1 
library(lubridate) #lubridate_1.9.4

ww_data_all <- read_csv("./section_4/data/ww_data_all.csv")

compile(file="./section_4/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored.cpp")
try(dyn.unload(dynlib("./section_4/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored")),silent = TRUE)
dyn.load(dynlib("./section_4/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored"))

compute_weights_precision <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}
local_poly_helper <- function (knots, refined_x, p = 2) {
  if (min(knots) >= 0) {
    dif <- diff(knots)
    nn <- length(refined_x)
    n <- length(knots)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x[j] <= knots[i]) {
          D[j, i] <- 0
        }
        else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= 
                 knots[i]) {
          D[j, i] <- (1/factorial(p)) * (refined_x[j] - 
                                           knots[i])^p
        }
        else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - 
                                          knots[i + 1])^(p - k))/(factorial(k) * factorial(p - 
                                                                                             k)))
        }
      }
    }
  }
  else if (max(knots) <= 0) {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D[j, i] <- 0
        }
        else if (refined_x_neg[j] <= knots_neg[i + 1] & 
                 refined_x_neg[j] >= knots_neg[i]) {
          D[j, i] <- (1/factorial(p)) * (refined_x_neg[j] - 
                                           knots_neg[i])^p
        }
        else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - 
                                          knots_neg[i + 1])^(p - k))/(factorial(k) * 
                                                                        factorial(p - k)))
        }
      }
    }
  }
  else {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D1 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D1[j, i] <- 0
        }
        else if (refined_x_neg[j] <= knots_neg[i + 1] & 
                 refined_x_neg[j] >= knots_neg[i]) {
          D1[j, i] <- (1/factorial(p)) * (refined_x_neg[j] - 
                                            knots_neg[i])^p
        }
        else {
          k <- 1:p
          D1[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - 
                                           knots_neg[i + 1])^(p - k))/(factorial(k) * 
                                                                         factorial(p - k)))
        }
      }
    }
    refined_x_pos <- refined_x
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    dif <- diff(knots_pos)
    nn <- length(refined_x_pos)
    n <- length(knots_pos)
    D2 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_pos[j] <= knots_pos[i]) {
          D2[j, i] <- 0
        }
        else if (refined_x_pos[j] <= knots_pos[i + 1] & 
                 refined_x_pos[j] >= knots_pos[i]) {
          D2[j, i] <- (1/factorial(p)) * (refined_x_pos[j] - 
                                            knots_pos[i])^p
        }
        else {
          k <- 1:p
          D2[j, i] <- sum((dif[i]^k) * ((refined_x_pos[j] - 
                                           knots_pos[i + 1])^(p - k))/(factorial(k) * 
                                                                         factorial(p - k)))
        }
      }
    }
    D <- cbind(D1, D2)
  }
  D
}
global_poly_helper <- function (x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}
prior_conversion <- function (d, prior, p) {
  Cp <- (d^((2 * p) - 1))/(((2 * p) - 1) * (factorial(p - 1)^2))
  prior_q <- list(a = prior$a, u = (prior$u * (1/sqrt(Cp))))
  prior_q
}

raw_d <- ww_data_all %>% 
  rename(sample_date = Collected,
         site_id = SampleLocation) 

# clean column names
work_d <- raw_d %>% 
  janitor::clean_names() %>% 
  group_by(sample_date) %>% 
  mutate(sample_date = ymd(sample_date)) %>% 
  ungroup()

work_d$sample_date %>% range()

work_d <- work_d %>% 
  filter(!is.na(sars_gcl)) %>% #not sure if these were measured or just data error
  mutate(censored_y = ifelse(sars_gcl <= 500, TRUE, FALSE)) %>% 
  mutate(sars_gcl = ifelse(censored_y, 500, sars_gcl))

## Merge site populations
sites <- read_csv("./section_4/data/sites.csv")
sites <- sites %>% 
  dplyr::select(SampleLocation,DisplayName,SampleType,Region, Population) %>% 
  rename("site_id"= SampleLocation,
         "pop"=Population) %>% 
  janitor::clean_names() 

work_d <- work_d %>% 
  left_join(sites, by = "site_id")

work_d$region %>% unique() 


region_pop <- data.frame(region = c("Auckland","Bay of Plenty","Canterbury","Gisborne","Hawke's Bay",
                                    "Manawatu-Whanganui","Marlborough","Nelson","Northland","Otago",
                                    'Southland','Taranaki','Tasman','Waikato','Wellington','West Coast'),
                         pop_region_2023 = c(1.716*1000000,347700, 666300,38200,184800,
                                             260900, 52200, 55600,203900, 254600,
                                             103900,127300,58700,513800,550600,32900))


work_d_full <- work_d 

work_d_comp <- work_d %>% 
  filter(sample_date >= "2021-11-12" & sample_date <= "2023-03-31")




## Run all before. 
source("./functions_general/prep_data_original.R")

# You might not have enough vector memory space to run this. 
data_foranalysis <- prep_data(outcome_column_name = "sars_gcl",
                              denom_column_name = "denom",
                              site_id = "site_id",
                              sample_date = "sample_date", 
                              data = work_d_comp %>% mutate(denom = 1), 
                              polyOrder = 3,
                              pred_also = FALSE)

data_foranalysis$df_full$site_id %>% unique() %>% length() #134 stations


tmbdat <- data_foranalysis$tmbdat

prior_IWP <- prior_conversion(d=20, prior =list(u=log(2),alpha = 0.5),p=3)

# Set other priors
tmbdat$u1 = prior_IWP$u
tmbdat$alpha1 = 0.5
tmbdat$u2 = 3
tmbdat$alpha2 = 0.5
tmbdat$betaprec = 0.01
tmbdat$lambda_phi = -log(0.5)/50
tmbdat$lambda_tau = -log(0.5)/0.5
tmbdat$lambda_cov = -log(0.5)/0.5

set.seed(2)
init_daily <- rnorm(ncol(tmbdat$daily), 0, 0.1);
init_W <- rnorm(ncol(tmbdat$obs),0,0.1)
tmbparams <- list(
  W = c(rep(0, (ncol(tmbdat$X)+ ncol(tmbdat$B))), init_daily, init_W, 0, diff(init_W)), # W = c(U,beta,Z); 
  theta1 = 10, # -2log(sigma)
  theta2 = 3,
  cov_log = 0,
  theta3 = 2.2,
  theta4 = -4.8
)


ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored",
  silent = TRUE
)


aghq_k = 3
start.time <- Sys.time()
mdl1 <- aghq::marginal_laplace_tmb(ff,k=aghq_k,startingvalue = c(10,0,0,0,0))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
samps1 <- aghq::sample_marginal(mdl1, M = 3000) 
marginals <- mdl1$marginals  

df_full =data_foranalysis$df_full

save(file="./section_4/wastewater_model_newzealand_censored_comp_fullAR.RData", list = c("samps1","marginals","tmbdat","df_full"))


load("./section_4/wastewater_model_newzealand_censored_comp_fullAR.RData")
source("./functions_general/process_results_posteriorsamps.R")

# You might not have enough vector memory space to run this.
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
results <- process_results(df_full =df_full,
                           tmbdat = tmbdat, 
                           samps1 = samps1, 
                           polyOrder=3, 
                           AR=TRUE)


save(file="./section_4/results_model_NZ.RData", list = "results")
