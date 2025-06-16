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
compute_post_fun <- function (samps, global_samps = NULL, knots, refined_x, p, degree = 0) {
  if (p <= degree) {
    return(message("Error: The degree of derivative to compute is not defined. Should consider higher order smoothing model or lower order of the derivative degree."))
  }
  if (is.null(global_samps)) {
    global_samps = matrix(0, nrow = p, ncol = ncol(samps))
  }
  if (nrow(global_samps) != p | nrow(samps) != (length(knots) - 
                                                1)) {
    return(message("Error: Incorrect dimension of global_samps or samps. Check whether the choice of p or the choice of knots are consistent with the fitted model."))
  }
  if (ncol(samps) != ncol(global_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and global polynomials."))
  }
  X = global_poly_helper(refined_x, p = p)
  X <- as.matrix(X[, 1:(p - degree)])
  for (i in 1:ncol(X)) {
    X[, i] <- (factorial(i + degree - 1)/factorial(i - 1)) * 
      X[, i]
  }
  B = as(local_poly_helper(knots, refined_x = refined_x, p = (p - 
                                                                degree)), "dgTMatrix")
  fitted_samps_deriv <- X %*% global_samps[(1 + degree):p, 
  ] + B %*% samps
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}

process_results <- function(df_full, tmbdat, samps1, polyOrder,  id_group, id_group_name = NULL, weights_df){
  
  P = as.matrix(tmbdat$P)
  X = as.matrix(tmbdat$X)
  Xfstat = as.matrix(tmbdat$Xfstat)
  daily = as.matrix(tmbdat$daily)
  obs = as.matrix(tmbdat$obs)
  knots = tmbdat$knots
  
  coefsamps1 <- samps1$samps[1:ncol(P),]
  global_samps1 <- samps1$samps[(ncol(P) + 1):(ncol(P) + ncol(X)-ncol(Xfstat)),]
  Xf_samps1 <- samps1$samps[(ncol(P) + ncol(X)-ncol(Xfstat) + 1):(ncol(P) + ncol(X)),]
  daily_samps1 <- samps1$samps[(ncol(P) + ncol(X) + 1):(ncol(P) + ncol(X) + ncol(daily)),]
  u_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+1):(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)),]
  u_deriv_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)+1):(nrow(samps1$samps)),]
  
  v <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                         knots = knots, 
                         refined_x = df_full$t,
                         p = polyOrder, degree = 0)
  
  vderiv <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                              knots = knots, 
                              refined_x = df_full$t,
                              p = polyOrder, degree = 1)
  
  Xfstat_full <-model.matrix(~site_id, 
                             data = df_full %>% mutate(site_id = factor(site_id)), 
                             contrasts.arg = list(site_id = "contr.sum"))[,-1]
  
  # norm <- Xfstat_full
  # norm[,2] <- -1
  # norm[,1]<- -1
  
  if (is.null(dim(Xf_samps1))){Xf_samps1 <- t(Xf_samps1)}
  
  v_u_fixed <- v[,-1]+ u_samps1 + Xfstat_full%*%Xf_samps1
  v_u <- v[,-1]+ u_samps1 
  v_u_deriv <-vderiv[,-1]+ u_deriv_samps1
  v_u_fixed_deriv <-vderiv[,-1]+ u_deriv_samps1
  u_fixed <- u_samps1 + Xfstat_full%*%Xf_samps1
  v_fixed <- v[,-1] + Xfstat_full%*%Xf_samps1
  v_fixed_deriv <- vderiv[,-1]
  
  ## Nominal AR2 + ospline+ fixed effects
  df_full$exp_v_u_fixed <- as.numeric(apply(exp(v_u_fixed), MARGIN=1,median))
  df_full$exp_v_u_fixed_lwr <- as.numeric(apply(exp(v_u_fixed), MARGIN=1,quantile, 0.025))
  df_full$exp_v_u_fixed_upr <- as.numeric(apply(exp(v_u_fixed), MARGIN=1,quantile, 0.975))
  
  df_full$exp_v_u <- as.numeric(apply(exp(v_u), MARGIN=1,median))
  df_full$exp_v_u_lwr <- as.numeric(apply(exp(v_u), MARGIN=1,quantile, 0.025))
  df_full$exp_v_u_upr <- as.numeric(apply(exp(v_u), MARGIN=1,quantile, 0.975))
  
  df_full$exp_v_u_fixed_deriv <- as.numeric(apply(exp(v_u_fixed)*(v_u_deriv), MARGIN=1,median))
  df_full$exp_v_u_fixed_deriv_lwr <- as.numeric(apply(exp(v_u_fixed)*(v_u_deriv), MARGIN=1,quantile, 0.025))
  df_full$exp_v_u_fixed_deriv_upr <- as.numeric(apply(exp(v_u_fixed)*(v_u_deriv), MARGIN=1,quantile, 0.975))
  
  
  ## Nominal AR2 + fixed effects
  df_full$exp_u_fixed <- as.numeric(apply(exp(u_fixed), MARGIN=1,median))
  df_full$exp_u_fixed_lwr <- as.numeric(apply(exp(u_fixed), MARGIN=1,quantile, 0.025))
  df_full$exp_u_fixed_upr <- as.numeric(apply(exp(u_fixed), MARGIN=1,quantile, 0.975))
  
  df_full$exp_u_fixed_deriv <- as.numeric(apply(exp(u_fixed)*(u_deriv_samps1), MARGIN=1,median))
  df_full$exp_u_fixed_deriv_lwr <- as.numeric(apply(exp(u_fixed)*(u_deriv_samps1), MARGIN=1,quantile, 0.025))
  df_full$exp_u_fixed_deriv_upr <- as.numeric(apply(exp(u_fixed)*(u_deriv_samps1), MARGIN=1,quantile, 0.975))
  
  ## Log AR2 + ospline+ fixed effects
  df_full$v_u_fixed <- as.numeric(apply(v_u_fixed, MARGIN=1,median))
  df_full$v_u_fixed_lwr <- as.numeric(apply(v_u_fixed, MARGIN=1,quantile, 0.025))
  df_full$v_u_fixed_upr <- as.numeric(apply(v_u_fixed, MARGIN=1,quantile, 0.975))
  
  df_full$v_u_fixed_deriv <- as.numeric(apply(v_u_deriv, MARGIN=1,median))
  df_full$v_u_fixed_deriv_lwr <- as.numeric(apply(v_u_deriv, MARGIN=1,quantile, 0.025))
  df_full$v_u_fixed_deriv_upr <- as.numeric(apply(v_u_deriv, MARGIN=1,quantile, 0.975))
  
  ## Log AR2 + fixed effects
  df_full$u_fixed <- as.numeric(apply(u_fixed, MARGIN=1,median))
  df_full$u_fixed_lwr <- as.numeric(apply(u_fixed, MARGIN=1,quantile, 0.025))
  df_full$u_fixed_upr <- as.numeric(apply(u_fixed, MARGIN=1,quantile, 0.975))
  
  df_full$u_fixed_deriv <- as.numeric(apply(u_deriv_samps1, MARGIN=1,median))
  df_full$u_fixed_deriv_lwr <- as.numeric(apply(u_deriv_samps1, MARGIN=1,quantile, 0.025))
  df_full$u_fixed_deriv_upr <- as.numeric(apply(u_deriv_samps1, MARGIN=1,quantile, 0.975))
  
  ## Log Ospline + fixed effects
  df_full$v_fixed <- as.numeric(apply(v_fixed, MARGIN=1,median))
  df_full$v_fixed_upr <- as.numeric(apply(v_fixed, MARGIN=1,quantile,p=0.975))
  df_full$v_fixed_lwr <- as.numeric(apply(v_fixed, MARGIN=1,quantile,p=0.025))
  
  ## Ospline + fixed effects
  df_full$exp_v_fixed <- as.numeric(apply(exp(v_fixed), MARGIN=1,median))
  df_full$exp_v_fixed_upr <- as.numeric(apply(exp(v_fixed), MARGIN=1,quantile,p=0.975))
  df_full$exp_v_fixed_lwr <- as.numeric(apply(exp(v_fixed), MARGIN=1,quantile,p=0.025))
  
  df_full$exp_v_fixed_deriv <- as.numeric(apply(exp(v_fixed)*vderiv[,-1], MARGIN=1,median))
  df_full$exp_v_fixed_deriv_upr <- as.numeric(apply(exp(v_fixed)*vderiv[,-1], MARGIN=1,quantile, p=0.975))
  df_full$exp_v_fixed_deriv_lwr<- as.numeric(apply(exp(v_fixed)*vderiv[,-1], MARGIN=1,quantile, p=0.025))
  
  ## Log Ospline 
  df_full$v <- as.numeric(apply(v[,-1], MARGIN=1,median))
  df_full$v_upr <- as.numeric(apply(v[,-1], MARGIN=1,quantile,p=0.975))
  df_full$v_lwr <- as.numeric(apply(v[,-1], MARGIN=1,quantile,p=0.025))
  
  df_full$v_deriv <- as.numeric(apply(vderiv[,-1], MARGIN=1,median))
  df_full$v_deriv_upr <- as.numeric(apply(vderiv[,-1], MARGIN=1,quantile, p=0.975))
  df_full$v_deriv_lwr<- as.numeric(apply(vderiv[,-1], MARGIN=1,quantile, p=0.025))
  
  ## Ospline
  df_full$exp_v <- as.numeric(apply(exp(v[,-1]), MARGIN=1,median))
  df_full$exp_v_upr <- as.numeric(apply(exp(v[,-1]), MARGIN=1,quantile, p = 0.975))
  df_full$exp_v_lwr <- as.numeric(apply(exp(v[,-1]), MARGIN=1,quantile, p = 0.025))
  
  df_full$exp_v_deriv <- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,median))
  df_full$exp_v_deriv_upr <- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,quantile, 0.975))
  df_full$exp_v_deriv_lwr<- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,quantile, 0.025))
  
  df_full$inst_repro = as.numeric(apply(exp(vderiv[,-1]), MARGIN=1,median))
  df_full$inst_repro_lwr = as.numeric(apply(exp(vderiv[,-1]), MARGIN=1,quantile, 0.025))
  df_full$inst_repro_upr = as.numeric(apply(exp(vderiv[,-1]), MARGIN=1,quantile, 0.975))
  
  df_full$inst_repro_v_u = as.numeric(apply(exp(v_u_fixed_deriv), MARGIN=1,median))
  df_full$inst_repro_v_u_lwr = as.numeric(apply(exp(v_u_fixed_deriv), MARGIN=1,quantile, 0.025))
  df_full$inst_repro_v_u_upr = as.numeric(apply(exp(v_u_fixed_deriv), MARGIN=1,quantile, 0.975))
  
if ( (1 %in% id_group) | (2 %in% id_group) ){
  post_samps_df <- df_full %>% 
    dplyr::select('sample_date','site_id',id_group_name) %>% 
    cbind(as.data.frame(v_u_fixed)) %>% 
    melt(id.vars = 1:(2+length(id_group_name)))
  
  post_samps_df_deriv <- df_full %>% 
    dplyr::select('sample_date','site_id',id_group_name) %>% 
    cbind(as.data.frame(v_u_fixed_deriv)) %>% 
    melt(id.vars = 1:(2+length(id_group_name)))
  
  post_samps_df_v_fixed <- df_full %>% 
    dplyr::select('sample_date','site_id',id_group_name) %>% 
    cbind(as.data.frame(v_fixed)) %>% 
    melt(id.vars = 1:(2+length(id_group_name)))
  
  post_samps_df_v_fixed_deriv <- df_full %>% 
    dplyr::select('sample_date','site_id',id_group_name) %>% 
    cbind(as.data.frame(v_fixed_deriv)) %>% 
    melt(id.vars = 1:(2+length(id_group_name)))
  
  post_samps_df_noint <- df_full %>% 
    dplyr::select('sample_date','site_id',id_group_name) %>% 
    cbind(as.data.frame(v_u)) %>% 
    melt(id.vars = 1:(2+length(id_group_name)))
  
  output_postsamps = post_samps_df %>% 
    rename("v_u_fixed" = value) %>% 
    cbind(post_samps_df_deriv%>% dplyr::select(value) %>% rename("v_u_fixed_deriv"= value) ) %>% 
    cbind(post_samps_df_v_fixed%>% dplyr::select(value) %>% rename("v_fixed"= value) ) %>% 
    cbind(post_samps_df_v_fixed_deriv%>% dplyr::select(value) %>% rename("v_fixed_deriv"= value) ) %>% 
    cbind(post_samps_df_noint %>% dplyr::select(value) %>% rename("v_u"= value) )}
  
if (1 %in% id_group){
  tmp<- output_postsamps  %>% 
    left_join(weights_df, by = "site_id") %>% 
    group_by(variable, sample_date) %>% 
    mutate(denom_weights = sum(weights)) %>% 
    summarise(ave_exps_v_u_fixed = sum(exp(v_u_fixed)*weights/denom_weights),
              ave_exps_v_u_fixed_deriv = sum(v_u_fixed_deriv*exp(v_u_fixed)*weights/denom_weights),
              ave_exps_v_fixed = sum(exp(v_fixed)*weights/denom_weights),
              ave_exps_v_fixed_deriv = sum(v_fixed_deriv*exp(v_fixed)*weights/denom_weights)) %>% 
    mutate(inst_repro = exp(ave_exps_v_u_fixed_deriv/ave_exps_v_u_fixed),
           cum_inst_repro = cumprod(inst_repro)) %>% 
    ungroup() %>% 
    group_by(sample_date) %>% 
    summarise(ave_exp_v_u_fixed_med = median(ave_exps_v_u_fixed),
              ave_exp_v_u_fixed_upr = quantile(ave_exps_v_u_fixed, 0.975),
              ave_exp_v_u_fixed_lwr = quantile(ave_exps_v_u_fixed, 0.025),
              ave_exp_v_u_fixed_deriv_med = median(ave_exps_v_u_fixed_deriv),
              ave_exp_v_u_fixed_deriv_upr = quantile(ave_exps_v_u_fixed_deriv, 0.975),
              ave_exp_v_u_fixed_deriv_lwr = quantile(ave_exps_v_u_fixed_deriv, 0.025),
              ave_exp_v_fixed_med = median(ave_exps_v_fixed),
              ave_exp_v_fixed_upr = quantile(ave_exps_v_fixed, 0.975),
              ave_exp_v_fixed_lwr = quantile(ave_exps_v_fixed, 0.025),
              ave_exp_v_fixed_deriv_med = median(ave_exps_v_fixed_deriv),
              ave_exp_v_fixed_deriv_upr = quantile(ave_exps_v_fixed_deriv, 0.975),
              ave_exp_v_fixed_deriv_lwr = quantile(ave_exps_v_fixed_deriv, 0.025),
              inst_repro_med = median(inst_repro),
              inst_repro_upr = quantile(inst_repro, 0.975),
              inst_repro_lwr = quantile(inst_repro, 0.025),
              prod_inst_repro_med = median(cum_inst_repro),
              prod_inst_repro_upr = quantile(cum_inst_repro, 0.975),
              prod_inst_repro_lwr = quantile(cum_inst_repro, 0.025))}else{tmp = NULL}
  
  if (2 %in% id_group){
    tmp2<- output_postsamps  %>% 
      rename(grouping_var = id_group_name) %>% 
      group_by(variable, sample_date, grouping_var) %>% 
      summarise(ave_exps = mean(exp(v_u_fixed)),
                ave_exps_deriv = mean(v_u_fixed_deriv*exp(v_u_fixed)),
                ave_exps_noint = mean(exp(v_u))) %>% 
      ungroup() %>% 
      group_by(sample_date, grouping_var) %>% 
      summarise(ave_exp_v_u_fixed = median(ave_exps),
                ave_exp_v_u_fixed_upr = quantile(ave_exps, 0.975),
                ave_exp_v_u_fixed_lwr = quantile(ave_exps, 0.025),
                ave_exp_v_u = median(ave_exps_noint),
                ave_exp_v_u_upr = quantile(ave_exps_noint, 0.975),
                ave_exp_v_u_lwr = quantile(ave_exps_noint, 0.025),
                ave_exp_v_u_fixed_deriv = median(ave_exps_deriv),
                ave_exp_v_u_fixed_deriv_upr = quantile(ave_exps_deriv, 0.975),
                ave_exp_v_u_fixed_deriv_lwr = quantile(ave_exps_deriv, 0.025),
                post_prob_ave_exp_v_u_fixed_deriv = length(which(ave_exps_deriv>0))/length(ave_exps_deriv))}else{tmp2 = NULL}
  
  return(list(df_full=df_full, station_ave_df = tmp,custom_ave_df = tmp2))
}