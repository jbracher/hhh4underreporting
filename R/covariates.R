#' Simulate from underreported model with one covariate
#'
#' Simulate from an endemic-epidemic model with a covariate entering
#' into the endemic component. Requires specification of the initial value for $lambda$.
#'
#' @param lambda1 the initial mean
#' @param alpha_nu,beta_nu,phi,kappa,psi the parameters of the latent model
#' @param q the reporting probability
#' @param z a covariate entering into the endemic parameter $nu$
#' @return A named list with elements \code{"X"} and \code{"Y"} containing the unthinned and thinned simulated time series.
#' @examples
#' # define a covariate:
#' z <- sin(2*pi*1:200/20)
#' # define the model parameters:
#' lambda1 <- 10
#' alpha_nu <- 3
#' beta_nu <- 1
#' phi <- 0.5
#' kappa <- 0.2
#' psi <- 0.1
#' q <- 0.5
#' sim <- generate_ar_cov(lambda1 = lambda1,
#'                        alpha_nu = alpha_nu,
#'                        beta_nu = beta_nu,
#'                        phi = phi,
#'                        kappa = kappa,
#'                        psi = psi,
#'                        q = q,
#'                        z = z)
#' @export
simulate_hhh4u_covariate <- function(lambda1, alpha_nu, beta_nu, phi, kappa, psi, q = 1, z){
  lgt = length(z)

  # compute nu:
  nu <- exp(alpha_nu + beta_nu*z)

  lambda <- X <- Y <- numeric(lgt)
  lambda[1] <- lambda1
  X[1] <- rnbinom(1, mu = lambda[1], size = 1/psi)

  for(t in 2:lgt){
    lambda[t] <- nu[t] + phi*X[t - 1] + kappa*lambda[t - 1]
    X[t] <- rnbinom(1, mu = lambda[t], size = 1/psi)
  }

  X_tilde <- rbinom(n = lgt, size = X, prob = q)
  return(list(X = X, X_tilde = X_tilde))
}

cond_mean_cov <- function(m1, nu, phi, kappa){
  lgt <- length(nu)
  mu <- numeric(lgt)
  mu[1] <- m1
  for(i in 2:lgt){
    mu[i] <- nu[i] + (phi[i] + kappa[i])*mu[i - 1]
  }
  return(mu)
}

compute_sop_cov <- function(m1, vl1, nu, phi, kappa, psi, q, compute_Sigma = FALSE){
  lgt <- length(nu)
  # get means:
  mu_X <- cond_mean_cov(m1 = m1, nu = nu,  phi = phi, kappa = kappa)
  # compute variances of lambda as they will be required:
  v_lambda <- v_X <- cov1_X <- numeric(lgt)
  v_lambda[1] <- vl1
  v_X[1] <- (1 + psi[1])*vl1 + m1 + psi[1]*m1^2
  for(i in 2:lgt){
    v_lambda[i] <- phi[i]^2*v_X[i - 1] + (2*phi[i]*kappa[i] + kappa[i]^2)*v_lambda[i - 1]
    v_X[i] <- (1 + psi[i])*v_lambda[i] + mu_X[i] + psi[i]*mu_X[i]^2
    cov1_X[i] <- kappa[i]*v_lambda[i - 1] +
      phi[i]*v_X[i - 1] +
      (phi[i] + kappa[i])*mu_X[i - 1]^2 +
      nu[i]*mu_X[i - 1] -
      mu_X[i - 1]*mu_X[i]
  }

  if(compute_Sigma){
    # move to matrix
    cov_X <- matrix(nrow = lgt, ncol = lgt)
    # fill the diagonal:
    diag(cov_X) <- v_X
    # the first off-diagonal is tricky:
    for(i in 1:(lgt - 1)){
      cov_X[i, i + 1] <- kappa[i + 1]*v_lambda[i] +
        phi[i + 1]*v_X[i] +
        (phi[i + 1] + kappa[i + 1])*mu_X[i]^2 +
        nu[i + 1]*mu_X[i] -
        mu_X[i]*mu_X[i + 1]
    }
    # the others are simple
    for(ro in 1:(lgt - 2)){ # loop over rows
      for(co in (ro + 2):lgt){
        cov_X[ro, co] <- (phi[co] + kappa[co])*cov_X[ro, co - 1]
      }
    }
    cov_X_tilde <- q^2*cov_X + q*(1 - q)*diag(mu_X)
  }else{
    cov_X <- cov_X_tilde <- NULL
  }

  # thinning:
  mu_X_tilde <- q*mu_X
  v_X_tilde <- q^2*v_X + q*(1 - q)*mu_X
  cov1_X_tilde <- q^2*cov1_X

  return(list(mu_X = mu_X, v_X = v_X, cov1_X = cov1_X,
              decay_cov_X = phi + kappa,  cov_X = cov_X, v_lambda = v_lambda,
              mu_X_tilde = mu_X_tilde, v_X_tilde = v_X_tilde, cov1_X_tilde = cov1_X_tilde,
              decay_cov_X_tilde = phi + kappa, cov_X_tilde = cov_X_tilde))
}


reparam_cov <- function(m1, vl1, nu, phi, kappa, psi, q){
  lgt <- length(nu)
  # compute target second-order properties:
  target_sop <- compute_sop_cov(m1 = m1, vl1 = vl1, nu = nu, phi = phi, kappa = kappa, psi = psi, q = q)

  # now find a completely observed process with the same properties:
  nu_Y <- p*nu # known for theoretical reasons (?)
  phi_plus_kappa_Y <- phi + kappa # this, too.

  mu_Y <- target_sop$mu_Y
  phi_Y <- kappa_Y <- psi_Y <- v_lambda_Y <- numeric(lgt)

  v_lambda_Y[1] <- q^2*vl1
  psi_Y[1] <- psi[1]

  v_Y <- target_sop$v_Y # by definition
  cov1_Y <- target_sop$cov1_Y # by definition
  phi_plus_kappa_Y <- target_sop$decay_cov_Y

  # element [1, 1] is correct by construction (thinning property of NB), the others can be adapted step by step.
  # done in loop:
  for(i in 2:lgt){
    # correct phi_Y[i]
    phi_Y[i] <- (cov1_Y[i] - phi_plus_kappa_Y[i]*v_lambda_Y[i - 1] -
                      phi_plus_kappa_Y[i]*mu_Y[i - 1]^2 -
                      nu_Y[i]*mu_Y[i - 1] + mu_Y[i]*mu_Y[i - 1])/
      (v_Y[i - 1] - v_lambda_Y[i - 1])
    kappa_Y[i] <- phi_plus_kappa_Y[i] - phi_Y[i]
    # update v_lambda_Y:
    v_lambda_Y[i] <- phi_Y[i]^2*v_Y[i - 1] +
      (2*phi_Y[i]*kappa_Y[i] + kappa_Y[i]^2)*v_lambda_Y[i - 1]
    # now correct psi_Y[i]
    psi_Y[i] <- (v_Y[i] - v_lambda_Y[i] - mu_Y[i])/
      (mu_Y[i]^2 + v_lambda_Y[i])
  }

  # current_sop <- compute_sop_cov(m1 = mu_Y[1], vl1 = v_lambda_Y[1], nu = nu_Y,
  #                              phi = phi_Y, kappa = kappa_Y, psi = psi_Y, p = 1)
  # if(all.equal(target_sop$v_Y, current_sop$v_Y)) print("worked")

  return(list(m1 = mu_Y[1], vl1 = v_lambda_Y[1], nu = nu_Y, phi = phi_Y, kappa = kappa_Y,
              psi = psi_Y, q = 1))
}

# auxiliary function to get weight matrix for lags:
get_weight_matrix_cov <- function(phi, kappa, max_lag){
  if(length(phi) != length(kappa)) stop("phi and kappa need to be the same length.")
  lgt <- length(phi)

  wgts <- matrix(ncol = lgt, nrow = max_lag)
  for(i in (max_lag + 1):lgt){
    # fill in phi values
    wgts[, i] <- phi[rev(seq(to = i, length.out = max_lag))]
    # fill in kappa values, one more for each row which we move down
    for(j in 2:max_lag){
      wgts[j:max_lag, i] <- wgts[j:max_lag, i]*kappa[rev(seq(to = i, length.out = max_lag - j + 1))]
    }
  }
  t(wgts)
}

# function for obtaining "shifted" nu necessary for observation-driven formulation
nu_to_nu_star_cov <- function(nu, kappa, max_lag){
  lgt <- length(nu)
  nu_transformed <- nu
  for(i in 1:max_lag){
    nu_transformed <- nu + kappa*c(NA, head(nu_transformed, lgt - 1))
  }
  nu_transformed
}

#' Evaluate the approximate log-likelihood of an underreported endemic-epidemic
#' model with time-varying parameters
#'
#' The likelihood approximation is based on an approximation of the process by
#' a second-order equivalent process with complete reporting. Note that this is
#' an R version while a reimplementation in Rcpp is used in the optimization.
#'
#' @param Y a time series of counts (numeric vector)
#' @param m1 the initial mean, i.e. $E(lambda_1)$
#' @param vl1 the initial variance of lambda, i.e. $Var(lambda_1)$
#' @param nu,phi,kappa,psi the time-varying model parameters (vectors of same length)
#' @param psi overdispersion parameter (scalar)
#' @param q the assumed reporting probability
#' @param max_lag in evaluation of likelihood only lags up to max_lag are taken into account
#' @param return_contributions shall the log-likelihood contributions of each time point be
#' returned (as vector)?
#'
#' @return The log-likelihood as scalar or (if \code{return_contributions == TRUE}) the vector of
#' log-likelihood contributions.
#'
#' @export
llik_cov <- function(Y, m1, vl1, nu, phi, kappa, psi, q, max_lag = 5, return_contributions = FALSE){
  lgt <- length(Y)
  # get model matrix:
  mod_matr <- matrix(nrow = length(Y), ncol = max_lag)
  for(i in 1:max_lag){
    mod_matr[, i] <- c(rep(NA, i), head(Y, lgt - i))
  }

  # get corresponding parameters for unthinned process:
  pars_star <- reparam_cov(m1 = m1, vl1 = vl1, nu = nu, phi = phi,
                         kappa = kappa, psi = psi, q = q)
  m1_Y <- pars_star$m1
  vl1_Y <- pars_star$vl1
  nu_Y <- pars_star$nu
  phi_Y <- pars_star$phi
  kappa_Y <- pars_star$kappa
  psi_Y <- pars_star$psi

  nu_Y_star <- nu_to_nu_star_cov(nu = nu_Y, kappa = kappa_Y, max_lag = max_lag)

  # get weight matrix:
  weight_matrix <- get_weight_matrix_cov(phi = phi_Y, kappa = kappa_Y, max_lag = max_lag)

  # get likelihood:
  lambda <- rep(nu_Y_star, length.out = lgt) + rowSums(weight_matrix*mod_matr)
  llik <- dnbinom(Y, mu = lambda, size = 1/psi_Y,
                  log = TRUE)
  if(return_contributions) return(llik) else return(sum(llik[-(1:max_lag)]))
  return(llik)
}

# version where the cpp bits are held together by R code.
llik_cov_r_cpp <- function(Y, m1, vl1, nu, phi, kappa, psi, q, max_lag = 5){
  lgt <- length(Y)
  # get model matrix:
  mod_matr <- get_mod_matr_cpp(Y = Y, max_lag = max_lag)

  # get corresponding parameters for unthinned process:
  pars_star <- reparam_cov_cpp(m1 = m1, vl1 = vl1, nu = nu, phi = phi,
                             kappa = kappa, psi = psi, q = q)
  m1_Y <- pars_star$m1
  vl1_Y <- pars_star$vl1
  nu_Y <- pars_star$nu
  phi_Y <- pars_star$phi
  kappa_Y <- pars_star$kappa
  psi_Y <- pars_star$psi

  nu_Y_star <- nu_to_nu_star_cov_cpp(nu = nu_Y, kappa = kappa_Y)

  # get weight matrix:
  weight_matrix <- get_weight_matrix_cov_cpp(phi = phi_Y, kappa = kappa_Y, max_lag = max_lag)

  # get likelihood:
  lambda <- rep(nu_Y_star, length.out = lgt) + rowSums(weight_matrix*mod_matr)
  llik <- sum(dnbinom(Y[-(1:max_lag)], mu = lambda[-(1:max_lag)], size = 1/psi_Y[-(1:max_lag)],
                      log = TRUE))
  return(llik)
}

#' Fit an endemic-epidemic model with underreporting and one covariate
#'
#' Fits an endemic-epidemic model with underreporting and one covariate in the endemic component
#' using an approximative maximum likelihood scheme. The likelihood approximation is based on an
#' approximation of the process by a second-order equivalent process with complete reporting.
#' The reporting probability cannot be estimated from the data (in most cases these contain no
#' information on it) and thus needs to be specified in advance.
#'
#' @param Y a time series of counts (numeric vector)
#' @param q the assumed reporting probability
#' @param m1 the initial mean, i.e. $E(lambda_1)$
#' @param vl1 the initial variance of lambda, i.e. $Var(lambda_1)$
#' @param covariate the values of the covariate entering into the endemic component (numeric vector)
#' @param initial the initial value of the parameter vector passed to optim
#' (note: the function tries different starting values in any case)
#' @param max_lag in evaluation of likelihood only lags up to max_lag are taken into account
#' @return the return object from \code{optim} providing the maximum likelihood estimates
#' (mostly on the log scale).
#' @export
fit_hhh4u_covariate <- function(Y, q, m1, vl1, covariate,
                      initial = c(alpha_nu = 4, beta_nu = 0,
                                  alpha_phi = -1,
                                  alpha_kappa = -1, log_psi = -3),
                      max_lag = 10, iter_optim = 3, ...){

  nllik_vect <- function(pars){
    lgt <- length(Y)
    nu <- exp(pars["alpha_nu"] + pars["beta_nu"]*covariate)
    phi <- rep(exp(pars["alpha_phi"]), lgt)
    kappa <- rep(exp(pars["alpha_kappa"]), lgt)
    psi <- rep(exp(pars["log_psi"]), lgt)
    nllik_cov_cpp(Y = Y, m1 = m1, vl1 = vl1,
                nu = nu, phi = phi, kappa = kappa, psi = psi,
                q = q, max_lag = max_lag)
  }

  opt <- optim(par = initial, fn = nllik_vect, ...)
  for(i in 1:iter_optim){
    opt <- optim(par = opt$par, fn = nllik_vect, ...)
  }
  return(opt)
}

