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


process_results <- function(df_full, tmbdat, samps1, polyOrder, AR = FALSE,version =1 ){
  
  P = as.matrix(tmbdat$P)
  X = as.matrix(tmbdat$X)
  Xfstat = as.matrix(tmbdat$Xfstat)
  daily = as.matrix(tmbdat$daily)
  if (version == 2){
  obs = as.matrix(tmbdat$obs)}
  knots = tmbdat$knots
  
  coefsamps1 <- samps1$samps[1:ncol(P),]
  global_samps1 <- samps1$samps[(ncol(P) + 1):(ncol(P) + ncol(X)-ncol(Xfstat)),]
  Xf_samps1 <- samps1$samps[(ncol(P) + ncol(X)-ncol(Xfstat) + 1):(ncol(P) + ncol(X)),]
  daily_samps1 <- samps1$samps[(ncol(P) + ncol(X) + 1):(ncol(P) + ncol(X) + ncol(daily)),]
  if (AR == TRUE){
  if (version == 2){
  u_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+1):(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)),]
  u_deriv_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)+1):(nrow(samps1$samps)),]}else{
    u_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+1):(ncol(P) + ncol(X) + ncol(daily)+nrow(daily)),]
    u_deriv_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+nrow(daily)+1):(nrow(samps1$samps)),]  
  }
    }
  
  v <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                         knots = knots, 
                         refined_x = unique(df_full$t),
                         p = polyOrder, degree = 0)
  
  v_full <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                        knots = knots, 
                        refined_x = df_full$t,
                        p = polyOrder, degree = 0)
  
  Xfstat_full <-model.matrix(~site_id, 
                             data = df_full %>% mutate(site_id = factor(site_id)), 
                             contrasts.arg = list(site_id = "contr.sum"))[,-1]
  
  v <- v[,-1] 
  v_full <- v_full[,-1] 
  v_fixed <- v_full +  Xfstat_full%*%Xf_samps1
  
  if (AR == TRUE){
  v_u <- v_full + u_samps1;
  v_u_fixed <- v_full + u_samps1 + Xfstat_full%*%Xf_samps1}
  
  v <- v %>% 
    mutate(sample_date = df_full$sample_date[1:nrow(v)])
  
  v <- v %>% 
    melt(id.vars = c("sample_date"))
  
  v_fixed <- v_fixed %>% 
    mutate(sample_date = df_full$sample_date, 
           site_id = df_full$site_id) %>% 
    melt(id.vars = c("sample_date", "site_id"))
  
  if (AR == TRUE){
    v_u <- v_u %>% 
      mutate(sample_date = df_full$sample_date, 
             site_id = df_full$site_id) %>% 
      melt(id.vars = c("sample_date", "site_id"))
    
    v_u_fixed <- v_u_fixed %>% 
      mutate(sample_date = df_full$sample_date, 
             site_id = df_full$site_id) %>% 
      melt(id.vars = c("sample_date", "site_id"))
  } else {v_u = NULL; v_u_fixed = NULL}
  
  
  return(list(v = v, v_fixed=v_fixed, v_u = v_u, v_u_fixed = v_u_fixed))
}