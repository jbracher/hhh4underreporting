#' @useDynLib hhh4underreporting
#' @importFrom Rcpp sourceCpp

# get stationary mean in periodic case; there is also an Rcpp version,
# this is only to check results.
stat_mean_seas <- function(nu, phi, kappa){
  if(length(nu) != length(phi) | length(nu) != length(kappa)){
    stop("The provided nu, phi and kappa need to be of the same length.")
  }

  # length of period:
  L <- length(nu)
  # set up matrix of contributions of different nu to different means:
  contributions <- matrix(NA, ncol = L, nrow = L)
  # fill this matrix row-wise:
  for(i in 1:L){
    # order the nu and phi like this:
    # i, i - 1, ..., 1, L, L - 1, ..., i + 1
    # i.e. increasing by their distance from i when moving backwards
    # vector for re-ordering:
    order_back_in_time <- if(i == L){
      L:1
    }else{c(i:1, L:(i + 1))}
    # re-order nu and phi:
    nu_back_in_time <- nu[order_back_in_time]
    phi_back_in_time <- phi[order_back_in_time]
    kappa_back_in_time <- kappa[order_back_in_time]
    # now the weights (cumulative products of phi, see manuscript)
    # can be calculated easily:
    weights <- cumprod(c(1, (phi_back_in_time + kappa_back_in_time)[c(1:(L - 1))]))
    # fill line in contributions:
    contributions[i, ] <- weights * nu_back_in_time /
      (1 - prod(phi + kappa))
  }
  # sum up the contributions (row-wise):
  rowSums(contributions)
}

# get stationary variance in periodic case; there is also an Rcpp version,
# this is only to check results.
stat_var_seas <- function(nu, phi, kappa, psi){

  if(length(nu) != length(phi) | length(nu) != length(kappa) | length(nu) != length(psi)){
    stop("The provided nu, phi and kappa need to be of the same length.")
  }

  L <- length(nu)
  psi_shifted <- psi[c(L, 1:(L - 1))]
  mu <- stat_mean_seas(nu, phi, kappa)
  v <- psi*mu^2 + mu
  v_shifted <- v[c(L, 1:(L - 1))]
  h <- (phi + kappa)^2 + phi^2*psi_shifted
  summands <- phi^2*v_shifted

  # get all summands entering tinto the different var(lambda_t):
  contributions <- matrix(NA, ncol = L, nrow = L)
  for(i in 1:L){
    order_back_in_time <- if(i == L){
      L:1
    }else{c(i:1, L:(i + 1))}

    weights <- cumprod(c(1, h[order_back_in_time]))
    contributions[i, ] <- weights[1:L]*summands[order_back_in_time]/(1 - prod(weights))
  }
  # obtain var(lambda_t) as rowSums:
  stat_var_lambda <- rowSums(contributions)
  # move to var(X_t:)
  stat_var_X <- (1 + psi)*stat_var_lambda + v
  return(list(stat_var_lambda = stat_var_lambda,
              stat_var_X = stat_var_X))
}

# wrapper: get stationary properties:
compute_sop_seas <- function(nu, phi, kappa, psi, q){
  # mean:
  mu_X <- stat_mean_seas_cpp(nu, phi, kappa)
  # variance:
  sigma2_X <- stat_var_seas_cpp(nu, phi, kappa, psi)$stat_var_X
  # autocovariance for d = 1 (corresponds to g in manuscript)
  inds_shifted <- c(length(phi), 1:(length(phi) - 1))
  cov1_X <- phi*sigma2_X[inds_shifted] +
    kappa*(sigma2_X[inds_shifted] - mu_X[inds_shifted] - psi[inds_shifted]*mu_X[inds_shifted]^2)/
    (1 + psi[inds_shifted])
  # decay of the autocovariance function (h in the manuscript):
  decay_X <- phi + kappa

  # and for the thinned version:
  mu_X_tilde <- q*mu_X
  sigma2_X_tilde <- q^2*sigma2_X + q*(1 - q)*mu_X
  cov1_X_tilde <- q^2*cov1_X
  decay_X_tilde <- decay_X

  return(list(mu_X = mu_X, sigma2_X = sigma2_X,
              cov1_X = cov1_X, v_X = sigma2_X, decay_cov_X = decay_X,
              mu_X_tilde = mu_X_tilde, v_X_tilde = sigma2_X_tilde,
              cov1_X_tilde = cov1_X_tilde, decay_cov_X_tilde = decay_X_tilde))
}


reparam_seas <- function(nu, phi, kappa, psi, q){
  L <- length(nu)

  # if q == 1 nothing needs to be done:
  if(q == 1) return(list(nu = nu, phi = phi, kappa = kappa, psi = psi, q = 1))

  # the second-order properties we want to achieve:
  target_sop <- compute_sop_seas(nu, phi, kappa, psi, q)

  # now find a completely observed process Y with these properties:
  nu_Y <- q*nu # known for theoretical reasons
  phi_plus_kappa_Y <- target_sop$decay_cov_X_tilde # this, too.
  mu_Y <- target_sop$mu_X_tilde # by definition
  v_Y <- target_sop$v_X_tilde # by definition
  cov1_Y <- target_sop$cov1_X_tilde # by definition

  # for the remaining one the old values are ususally good starting values:
  phi_Y <- phi
  kappa_Y <- phi_plus_kappa_Y - phi_Y
  psi_Y <- psi

  # given a starting value of psi[L] we can get var_lambda_Y[L]
  # (the variance of the conditional mean process for Y)
  v_lambda_Y <- numeric(L)
  v_lambda_Y[L] <- (v_Y[L] - mu_Y[L] - psi_Y[L]*mu_Y[L]^2)/
    (1 + psi_Y[L])

  # now we can iteratively roll through the weeks and correct the values of phi and psi
  for(k in 1:3){
    for(i in 1:L){
      im1 <- ifelse(i == 1, L, i - 1)
      # correct phi_Y[i]
      phi_Y[i] <- (cov1_Y[i] - phi_plus_kappa_Y[i]*v_lambda_Y[im1] -
                        phi_plus_kappa_Y[i]*mu_Y[im1]^2 -
                        nu_Y[i]*mu_Y[im1] + mu_Y[i]*mu_Y[im1])/
        (v_Y[im1] - v_lambda_Y[im1])
      kappa_Y[i] <- phi_plus_kappa_Y[i] - phi_Y[i]
      # update v_lambda_Y:
      v_lambda_Y[i] <- phi_Y[i]^2*v_Y[im1] +
        (2*phi_Y[i]*kappa_Y[i] + kappa_Y[i]^2)*v_lambda_Y[im1]
      # now correct psi_Y[i]
      psi_Y[i] <- (v_Y[i] - v_lambda_Y[i] - mu_Y[i])/
        (mu_Y[i]^2 + v_lambda_Y[i])
    }
  }
  return(list(nu = nu_Y, phi = phi_Y, kappa = kappa_Y, psi = psi_Y, q = 1))
}

# function for obtaining nu_star (compare manuscript), necessary for geometric-lag formulation
nu_to_nu_star_seas <- function(nu, kappa, max_lag = 5){
  nu_prolonged <- c(tail(nu, max_lag), nu)
  nu_matr <- matrix(nrow = length(nu), ncol = max_lag + 1)
  for(i in 1:length(nu)){
    nu_matr[i, ] <- rev(nu_prolonged[seq(from = i, length.out = max_lag + 1)])
  }
  weight_matrix <- get_weight_matrix_seas_cpp(phi = rep(1, length(nu)), kappa = kappa, max_lag = max_lag + 1)
  nu_star <- rowSums(weight_matrix*nu_matr)
  return(nu_star)
}

# auxililary function to get lagged variables:
back_seas <- function(vect, to, lgt){
  L <- length(vect)
  vect <- rep(vect, ceiling(lgt / L) + 1)
  # order backwards:
  vect_back <- if(to == L){
    rev(tail(vect, lgt))
  }else{
    rev(head(tail(vect, L - to + lgt), lgt))
  }
  return(vect_back)
}

# auxiliary function to get weight matrix for lags:
get_weight_matrix_seas <- function(phi, kappa, max_lag){
  if(length(phi) != length(kappa)) stop("phi and kappa need to be the same length.")
  L <- length(phi)
  wgts <- matrix(ncol = L, nrow = max_lag)
  for(i in 1:L){
    # fill in phi values
    wgts[, i] <- back_seas(phi, to = i, lgt = max_lag)
    # fill in kappa values, one more for each row which we move down
    for(j in 2:max_lag){
      wgts[j:max_lag, i] <- wgts[j:max_lag, i]*back_seas(kappa, to = i, lgt = max_lag - j + 1)
    }
  }
  t(wgts)
}

#' Fits a sesaonal endemic-epidemic model with underreporting using an approximative maximum
#' likelihood scheme. The likelihood approximation is based on an approximation of the process by
#' a second-order equivalent process with complete reporting. The reporting probability cannot be
#' estimated from the data (in most cases these contain no information on it) and thus needs to be
#' specified in advance.
#'
#' @param observed a time series of counts (numeric vector)
#' @param L the number of observations per period (integer; e.g. 52 for weekly data and yearly seasonality)
#' @param q the assumed reporting probability
#' @param seas_phi should seasonality also be accounted for in the autoregressive parameter $phi$?
#' or only in the immigration parameter $nu$?
#' @param include_kappa should the parameter $kappa$ and thus lagged versions of the conditional expectation $lambda$ be included?
#' @param initial the initial value of the parameter vector passed to optim
#' (note: the function tries different starting values in any case)
#' @param max_lag in evaluation of likelihood only lags up to max_lag are taken into account
#' @return the return object from \code{optim} providing the maximum likelihood estimates
#' (mostly on the log scale).
#' @export
fit_hhh4u_seasonal <- function(observed, L, q, seas_phi = FALSE, include_kappa = TRUE,
                               initial = c(alpha_nu = 4, gamma_nu = 0, delta_nu = 0,
                                              alpha_phi = -1, gamma_phi = 0, delta_phi = 0,
                                              alpha_kappa = -1, log_psi = -3),
                         max_lag = 10, iter_optim = 3, ...){

  # adapt initial values depending on which aparameters are included in the model:
  if(seas_phi == FALSE){
    initial <- initial[!(names(initial) %in% c("gamma_phi", "delta_phi"))]
  }
  if(include_kappa == FALSE){
    initial <- initial[names(initial) != "alpha_kappa"]
  }

  # wrapper around nllik_seas, to be passed to optim
  lik_vect <- function(pars){
    alpha_nu <- pars["alpha_nu"]
    gamma_nu <- pars["gamma_nu"]
    delta_nu <- pars["delta_nu"]
    alpha_phi <- pars["alpha_phi"]
    # use gamma_phi and delta_phi only if seas_phi == TRUE:
    gamma_phi <- ifelse(seas_phi, pars["gamma_phi"], 0)
    delta_phi <- ifelse(seas_phi, pars["delta_phi"], 0)
    alpha_kappa <- ifelse(include_kappa, pars["alpha_kappa"], -10)
    psi <- exp(pars["log_psi"])
    nllik_seas(observed = observed,
             alpha_nu = alpha_nu, gamma_nu = gamma_nu, delta_nu = delta_nu,
             alpha_phi = alpha_phi, gamma_phi = gamma_phi, delta_phi = delta_phi,
             alpha_kappa = alpha_kappa, psi = psi, q = q, L = L, max_lag = max_lag)
  }

  # optimize:
  opt <- optim(par = initial, fn = lik_vect,...)
  # repeat optimizatin (improves convergence)
  for(i in 1:iter_optim){
    opt <- optim(par = opt$par, fn = lik_vect,...)
  }

  return(opt)
}

#' Evaluate the approximate negative log-likelihood of a seasonal underreported endemic-epidemic model
#'
#' The likelihood approximation is based on an approximation of the process by
#' a second-order equivalent process with complete reporting.
#'
#' @param observed a time series of counts (numeric vector)
#' @param alpha_nu,gamma_nu,delta_nu the endemic model parameters (scalar)
#' @param alpha_phi,gamma_phi,delta_phi,alpha_kappa the autoregressive model parameters (scalars)
#' @param psi overdispersion parameter (scalar)
#' @param q the assumed reporting probability
#' @param L the period length (integer; e.g. 52 for weekly data and yearly seasonality)
#' @param max_lag in evaluation of likelihood only lags up to max_lag are taken into account
#' @param return_contributions shall the log-likelihood contributions of each time point be
#' returned (as vector)?
#'
#' @return The negative log-likelihood as scalar or (if \code{return_contributions == TRUE}) the vector of
#' negative log-likelihood contributions.
nllik_seas <- function(observed, alpha_nu, gamma_nu, delta_nu,
                     alpha_phi, gamma_phi = 0, delta_phi = 0,
                     alpha_kappa, psi, q, L, max_lag = 5, return_contributions = FALSE){
  lgt <- length(observed)
  # get model matrix:
  # mod_matr <- matrix(nrow = length(observed), ncol = max_lag)
  # for(i in 1:max_lag){
  #   mod_matr[, i] <- c(rep(NA, i), head(observed, lgt - i))
  # }
  mod_matr <- get_mod_matr_cpp(observed = observed, max_lag = max_lag)

  # extract parameter values over one season:
  vect_t <- seq(from = 1, length.out = L)
  nu <- exp(alpha_nu + gamma_nu*sin(2*pi*vect_t/L) + delta_nu*cos(2*pi*vect_t/L))
  phi <- exp(alpha_phi + gamma_phi*sin(2*pi*vect_t/L) + delta_phi*cos(2*pi*vect_t/L))
  kappa <- rep(exp(alpha_kappa), L)
  psi <- rep(psi, L)

  # get second order properties:
  sop <- compute_sop_seas(nu = nu, phi = phi, kappa = kappa, psi = psi, q = q)

  if(any(unlist(sop) < 0)) return(-Inf)

  # get corresponding parameters for unthinned process:
  pars_Y <- reparam_seas_cpp(nu = nu, phi = phi, kappa = kappa, psi = psi, q = q)
  nu_Y <- pars_Y$nu
  phi_Y <- pars_Y$phi
  kappa_Y <- pars_Y$kappa
  psi_Y <- pars_Y$psi

  # nu_Y needs to be transformed to move to the observation-driven formulation
  # of our process
  nu_Y_tilde <- nu_to_nu_star_seas_cpp(nu = nu_Y, kappa = kappa_Y, max_lag = max_lag)

  # get weight matrix:
  weight_matrix <- get_weight_matrix_seas_cpp(phi = phi_Y, kappa = kappa_Y, max_lag = max_lag)
  # paste as appropriate.
  weight_matrix <- weight_matrix[rep(1:L, length.out = lgt), ]

  # get likelihood:
  lambda <- rep(nu_Y_tilde, length.out = lgt) + rowSums(weight_matrix*mod_matr)
  llik <- -1*dnbinom(observed, mu = lambda, size = rep(1/psi_Y, length.out = lgt),
                        log = TRUE)
  if(return_contributions) return(llik) else return(sum(llik[-(1:max_lag)]))
}

#' Simulate from seasonal underreported model
#'
#' Simulate from a seasonal endemic-epidemic model with underreporting. Includes a burn-in period
#' to reach stationary phase.
#'
#' @param alpha_nu,gamma_nu,delta_nu,alpha_phi,gamma_phi,delta_phi,kappa,psi the model parameters
#' @param q the reporting probability
#' @param L the length of one cycle (integer)
#' @param n_seas the number of seasons to simulate (integer)
#' @param start initial value of both $X$ and $lambda$ (at beginning of burn in period; can usually be ignored)
#' @param burn_in number of seasons to discard to reach stationary phase
#' @return A named list with elements \code{"X"} and \code{"Y"} containing the unthinned and thinned simulated time series.
#' @export
simulate_hhh4u_seasonal <- function(alpha_nu, gamma_nu, delta_nu,
                                    alpha_phi, gamma_phi, delta_phi,
                                    kappa, psi, q = 1, L = 52, n_seas = 10,
                             start = 10, burn_in = 10){

  lgt <- n_seas*L
  lgt_total <- lgt + burn_in*L

  sin_t <- sin(2*pi*1:lgt/L)
  cos_t <- cos(2*pi*1:lgt/L)

  nu <- exp(alpha_nu + gamma_nu*sin_t + delta_nu*cos_t)
  phi <- exp(alpha_phi + gamma_phi*sin_t + delta_phi*cos_t)

  lambda <- X <- rep(NA, lgt)
  lambda[1] <- X[1] <- start
  for(i in 2:lgt){
    ind <- if(i %% L == 0){
      L
    }else{
      i%%L
    }
    lambda[i] <- nu[ind] + phi[ind]*X[i - 1] + kappa*lambda[i - 1]
    X[i] <- rnbinom(1, mu = lambda[i], size = 1/psi)
  }
  X <- tail(X, lgt)
  lambda <- tail(lambda, lgt)
  X_tilde <- rbinom(lgt, X, q)
  list(X = X, lambda = lambda, X_tilde = X_tilde)
}

# helper function to obtain second-order properties of simulated time series
emp_sop_seas <- function(observed, L){
  # bring in matrix form where rows correspond to calendar weeks
  # and columns to seasons.
  obs_matrix <- matrix(observed, nrow = L)

  means_emp <- rowMeans(Y_by_stratum)
  sq_means_emp <- rowMeans(Y_by_stratum^2)
  var_emp <- apply(Y_by_stratum, 1, var)

  ar1_emp <- numeric(L)
  ar1_emp[1] <- cov(Y_by_stratum[1, -1], Y_by_stratum[L, -ncol(Y_by_stratum)])

  for(i in 2:L){
    ar1_emp[i] <- cov(Y_by_stratum[i, ], Y_by_stratum[i - 1, ])
  }

  ar2_emp <- numeric(L)
  ar2_emp[1] <- cov(obs_matrix[1, -1], obs_matrix[L - 1, -ncol(obs_matrix)])
  ar2_emp[2] <- cov(obs_matrix[2, -1], obs_matrix[L, -ncol(obs_matrix)])

  for(i in 3:L){
    ar2_emp[i] <- cov(obs_matrix[i, ], obs_matrix[i - 2, ])
  }

  return(list(means_emp = means_emp,
              var_emp = var_emp,
              ar1_emp = ar1_emp,
              ar2_emp = ar2_emp))
}
