#' @useDynLib hhh4underreporting
#' @importFrom Rcpp sourceCpp

# get stationary mean in periodic case:
stat_mean_seas <- function(nu, phi, kappa){
  if(length(nu) != length(phi) | length(nu) != length(kappa)){
    stop("The provided nu, phi and kappa need to be of the same length.")
  }

  # length of period:
  L <- length(nu)
  # set up matrix of contributions of different g-values to different means:
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

# get stationary variance in periodic case:
stat_var_seas <- function(nu, phi, kappa, psi){

  if(length(nu) != length(phi) | length(nu) != length(kappa) | length(nu) != length(psi)){
    stop("The provided nu, phi and kappa need to be of the same length.")
  }

  L <- length(nu)
  a <- psi; b <- 1; c <- 0
  a_shifted <- a[c(L, 1:(L - 1))]
  mu <- stat_mean_seas(nu, phi, kappa)
  v <- a*mu^2 + b*mu + c
  v_shifted <- v[c(L, 1:(L - 1))]
  h <- (phi + kappa)^2 + phi^2*a_shifted
  summands <- phi^2*v_shifted

  contributions <- matrix(NA, ncol = L, nrow = L)
  for(i in 1:L){
    order_back_in_time <- if(i == L){
      L:1
    }else{c(i:1, L:(i + 1))}

    weights <- cumprod(c(1, h[order_back_in_time]))
    contributions[i, ] <- weights[1:L]*summands[order_back_in_time]/(1 - prod(weights))
  }
  stat_var_lambda <- rowSums(contributions)
  stat_var_X <- (1 + a)*stat_var_lambda + v
  return(list(stat_var_lambda = stat_var_lambda,
              stat_var_X = stat_var_X,
              h = h))
}

# wrapper: get stationary properties:
compute_sop_seas <- function(nu, phi, kappa, psi, p){
  mu_X <- stat_mean_seas(nu, phi, kappa)
  sigma2_X <- stat_var_seas(nu, phi, kappa, psi)$stat_var_X
  inds_shifted <- c(length(phi), 1:(length(phi) - 1))
  cov1_X <- phi*sigma2_X[inds_shifted] +
    kappa*(sigma2_X[inds_shifted] - mu_X[inds_shifted] - psi[inds_shifted]*mu_X[inds_shifted]^2)/
    (1 + psi[inds_shifted])
  decay_X <- phi + kappa

  mu_Y <- p*mu_X
  sigma2_Y <- p^2*sigma2_X + p*(1 - p)*mu_X # hier weitermachen!!! beide Manuskripte ansehen.
  cov1_Y <- p^2*cov1_X
  decay_Y <- decay_X

  return(list(mu_X = mu_X, sigma2_X = sigma2_X,
              cov1_X = cov1_X, v_X = sigma2_X, decay_cov_X = decay_X,
              mu_Y = mu_Y, v_Y = sigma2_Y,
              cov1_Y = cov1_Y, decay_cov_Y = decay_Y))
}


reparam_seas <- function(nu, phi, kappa, psi, p){
  L <- length(nu)

  if(p == 1) return(list(nu = nu, phi = phi, kappa = kappa, psi = psi, p = 1))

  target_sop <- compute_sop_seas(nu, phi, kappa, psi, p)

  # now find a completely observed process with the same properties:
  nu_star <- p*nu # known for theoretical reasons
  phi_plus_kappa_star <- target_sop$decay_cov_Y # this, too.
  mu_X_star <- target_sop$mu_Y # by definition
  v_X_star <- target_sop$v_Y # by definition
  cov1_X_star <- target_sop$cov1_Y # by definition

  # for the remaining one the old values are ususally good starting values:
  phi_star <- phi
  kappa_star <- phi_plus_kappa_star - phi_star
  psi_star <- psi

  # given a starting value of psi[L] we can get var_lambda_star[L]
  v_lambda_star <- numeric(L)
  v_lambda_star[L] <- (v_X_star[L] - mu_X_star[L] - psi_star[L]*mu_X_star[L]^2)/
    (1 + psi_star[L])

  # now we can iteratively roll through the weeks and correct the values of phi and kappa
  for(k in 1:3){
    for(i in 1:L){
      im1 <- ifelse(i == 1, L, i - 1)
      # correct phi_star[i]
      phi_star[i] <- (cov1_X_star[i] - phi_plus_kappa_star[i]*v_lambda_star[im1] -
                        phi_plus_kappa_star[i]*mu_X_star[im1]^2 -
                        nu_star[i]*mu_X_star[im1] + mu_X_star[i]*mu_X_star[im1])/
        (v_X_star[im1] - v_lambda_star[im1])
      kappa_star[i] <- phi_plus_kappa_star[i] - phi_star[i]
      # update v_lambda_star:
      v_lambda_star[i] <- phi_star[i]^2*v_X_star[im1] +
        (2*phi_star[i]*kappa_star[i] + kappa_star[i]^2)*v_lambda_star[im1]
      # now correct psi_star[i]
      psi_star[i] <- (v_X_star[i] - v_lambda_star[i] - mu_X_star[i])/
        (mu_X_star[i]^2 + v_lambda_star[i])
    }
  }
  return(list(nu = nu_star, phi = phi_star, kappa = kappa_star, psi = psi_star, p = 1))
}

# function for obtaining "shifted" nu necessary for observation-driven formulation
nu_to_nu_tilde_seas <- function(nu, kappa, max_lag = 5){
  nu_prolonged <- c(tail(nu, max_lag), nu)
  nu_matr <- matrix(nrow = length(nu), ncol = max_lag + 1)
  for(i in 1:length(nu)){
    nu_matr[i, ] <- rev(nu_prolonged[seq(from = i, length.out = max_lag + 1)])
  }
  weight_matrix <- get_weight_matrix_seas(phi = rep(1, length(nu)), kappa = kappa, max_lag = max_lag + 1)
  nu_transformed <- rowSums(weight_matrix*nu_matr)
  return(nu_transformed)
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

# function for likelihood inference
#' @export
fit_lik_seas <- function(Y, L, p, seas_phi = FALSE, initial = c(alpha_nu = 4, gamma_nu = 0, delta_nu = 0,
                                              alpha_phi = -1, gamma_phi = 0, delta_phi = 0,
                                              alpha_kappa = -1, log_psi = -3),
                         max_lag = 10, ...){

  if(seas_phi == FALSE){
    initial <- initial[c("alpha_nu", "gamma_nu", "delta_nu",
                         "alpha_phi", "alpha_kappa", "log_psi")]
  }

  lik_vect <- function(pars){
    alpha_nu <- pars["alpha_nu"]
    gamma_nu <- pars["gamma_nu"]
    delta_nu <- pars["delta_nu"]
    alpha_phi <- pars["alpha_phi"]
    # use gamma_phi and delta_phi only if seas_phi == TRUE:
    gamma_phi <- ifelse(seas_phi, pars["gamma_phi"], 0)
    delta_phi <- ifelse(seas_phi, pars["delta_phi"], 0)
    alpha_kappa <- pars["alpha_kappa"]
    psi <- exp(pars["log_psi"])
    lik_seas(Y = Y,
             alpha_nu = alpha_nu, gamma_nu = gamma_nu, delta_nu = delta_nu,
             alpha_phi = alpha_phi, gamma_phi = gamma_phi, delta_phi = delta_phi,
             alpha_kappa = alpha_kappa, psi = psi, p = p, L = L, max_lag = max_lag)
  }
  initials <- list(initial,
                   initial*c(1.5, rep(1, length(initial) - 1)),
                   initial*c(0.5, rep(1, length(initial) - 1)))
  opt <- optim(par = initials[[1]], fn = lik_vect,...)
  for(i in 2:3){
    opt_temp <- optim(par = initials[[i]], fn = lik_vect, ...)
    if(opt_temp$value < opt$value) opt <- opt_temp
  }
  return(opt)
}

# function to evaluate likelihood
lik_seas <- function(Y, alpha_nu, gamma_nu, delta_nu,
                     alpha_phi, gamma_phi = 0, delta_phi = 0,
                     alpha_kappa, psi, p, L, max_lag = 5){
  lgt <- length(Y)
  # get model matrix:
  mod_matr <- matrix(nrow = length(Y), ncol = max_lag)
  for(i in 1:max_lag){
    mod_matr[, i] <- c(rep(NA, i), head(Y, lgt - i))
  }

  # extract parameter values over one season:
  vect_t <- seq(from = 0, length.out = L)
  nu <- exp(alpha_nu + gamma_nu*sin(2*pi*vect_t/52) + delta_nu*cos(2*pi*vect_t/52))
  phi <- exp(alpha_phi + gamma_phi*sin(2*pi*vect_t/52) + delta_phi*cos(2*pi*vect_t/52))
  kappa <- rep(exp(alpha_kappa), L)
  psi <- rep(psi, L)

  # get second order properties:
  sop <- compute_sop_seas(nu = nu, phi = phi, kappa = kappa, psi = psi, p = p)

  if(any(unlist(sop) < 0)) return(-Inf)

  # get corresponding parameters for unthinned process:
  pars_star <- reparam_seas(nu = nu, phi = phi, kappa = kappa, psi = psi, p = p)
  nu_star <- pars_star$nu
  phi_star <- pars_star$phi
  kappa_star <- pars_star$kappa
  psi_star <- pars_star$psi

  # nu_star needs to be transformed to move to the observation-driven formulation
  # of our process
  nu_star_tilde <- nu_to_nu_tilde_seas(nu = nu_star, kappa = kappa_star, max_lag = max_lag)

  # get weight matrix:
  weight_matrix <- get_weight_matrix_seas(phi = phi_star, kappa = kappa_star, max_lag = max_lag)
  # paste as appropriate.
  weight_matrix <- weight_matrix[rep(1:L, length.out = lgt), ]

  # get likelihood:
  lambda <- rep(nu_star_tilde, length.out = lgt) + rowSums(weight_matrix*mod_matr)
  llik <- - sum(dnbinom(Y, mu = lambda, size = rep(1/psi_star, length.out = lgt),
                        log = TRUE)[-(1:max_lag)])
  return(llik)
}

# function to simulate from seasonal model:
generate_ar_seas <- function(nu, phi, kappa, psi, p = 1, n_seas = 10,
                             start = 10, burn_in = 10){
  L <- length(nu)
  lgt <- L*(n_seas + burn_in)

  if(length(nu) != length(phi) |
     length(nu) != length(psi) |
     length(nu) != length(kappa)){
    stop("nu and lambda, kappa and psi need to have the same length.")
  }

  L <- length(nu)
  lambda <- X <- rep(NA, lgt)
  lambda[1] <- X[1] <- start
  for(i in 2:lgt){
    ind <- if(i %% L == 0){
      L
    }else{
      i%%L
    }
    lambda[i] <- nu[ind] + phi[ind]*X[i - 1] + kappa[ind]*lambda[i - 1]
    X[i] <- rnbinom(1, mu = lambda[i], size = 1/psi[ind])
  }
  X <- tail(X, lgt - burn_in*L)
  lambda <- tail(lambda, lgt - burn_in*L)
  Y <- rbinom(lgt - burn_in*L, X, p)
  list(X = X, lambda = lambda, Y = Y)
}

# function to get some empirical second order properties
emp_sop_seas <- function(Y, L){
  Y_by_stratum <- matrix(Y, nrow = L)

  Y_means_emp <- rowMeans(Y_by_stratum)
  Y_sq_means_emp <- rowMeans(Y_by_stratum^2)
  Y_var_emp <- apply(Y_by_stratum, 1, var)

  Y_ar1_emp <- numeric(L)
  Y_ar1_emp[1] <- cov(Y_by_stratum[1, -1], Y_by_stratum[L, -ncol(Y_by_stratum)])

  for(i in 2:L){
    Y_ar1_emp[i] <- cov(Y_by_stratum[i, ], Y_by_stratum[i - 1, ])
  }

  Y_ar2_emp <- numeric(L)
  Y_ar2_emp[1] <- cov(Y_by_stratum[1, -1], Y_by_stratum[L - 1, -ncol(Y_by_stratum)])
  Y_ar2_emp[2] <- cov(Y_by_stratum[2, -1], Y_by_stratum[L, -ncol(Y_by_stratum)])

  for(i in 3:L){
    Y_ar2_emp[i] <- cov(Y_by_stratum[i, ], Y_by_stratum[i - 2, ])
  }

  return(list(Y_means_emp = Y_means_emp,
              Y_var_emp = Y_var_emp,
              Y_ar1_emp = Y_ar1_emp,
              Y_ar2_emp = Y_ar2_emp))
}
