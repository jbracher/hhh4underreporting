#' Simulate from time-homogeneous underreported model
#'
#' Simulate from a time-homogeneous endemic-epidemic model with underreporting. Includes a burn-in period
#' to reach stationary phase.
#'
#' @param nu,phi,kappa,psi the model parameters (scalars)
#' @param q the reporting probability
#' @param lgt the length of the generated time series (integer)
#' @param start initial value of both $X$ and $lambda$ (at beginning of burn in period)
#' @param burn_in number of iterations to discard to reach stationary phase
#' @return A named list with elements \code{"X"}, \code{"Y"} and \code{"lambda"} containing
#' the unthinned and thinned simulated time series as well as the conditional mean process.
#' @examples
#' sim <- simulate_hhh4u(nu = 3, phi = 0.4, kappa = 0.2, psi = 0.1, q = 0.5, lgt = 100)
#' @export
simulate_hhh4u <- function(nu, phi, kappa, psi, q = 1, lgt = 100,
                             start = 10, burn_in = 1000){
  simulate_hhh4u_seasonal(nu = nu, phi = phi, kappa = kappa, psi = psi, q = q,
                   n_seas = lgt, start = start, burn_in = burn_in)
}

# function to compute second-order properties
compute_sop <- function(nu, phi, kappa, psi, q, par_list = NULL){
  if(!is.null(par_list)){
    nu <- par_list$nu; phi <- par_list$phi; kappa <- par_list$kappa
    psi <- par_list$psi; q <- par_list$q
    names(nu) <- names(phi) <- names(kappa) <- names(psi) <- names(q) <- NULL
  }

  soc <- list()
  X <- X_tilde <- list()

  X$mu <- nu/(1 - phi - kappa)
  X$sigma2 <- (1 - (phi + kappa)^2 + phi^2)/
    (1 - (phi + kappa)^2 - psi*phi^2) *(X$mu + X$mu^2*psi)
  X$g <- phi*(1 - kappa*(phi + kappa))/(1 - (phi + kappa)^2 + phi^2)
  X$h <- phi + kappa

  X_tilde$mu <- q*X$mu
  X_tilde$sigma2 <- q^2*X$sigma2 + q*(1 - q)*X$mu
  X_tilde$g <- X$sigma2/(X$sigma2 + (1 - q)/q*X$mu)*X$g
  X_tilde$h <- phi + kappa

  return(list(X = X,
              X_tilde = X_tilde))
}

# function to recover nu, phi, kappa, psi for a given q:
recover_pars <- function(mu, sigma2, g, h, q = 1, sop_list = NULL){
  if(!is.null(sop_list)){
    mu <- sop_list$mu; sigma2 <- sop_list$sigma2; g <- sop_list$g; h <- sop_list$h
  }

  pars <- NULL
  pars$nu <- mu*(1 - h)/q
  c <- (sigma2 - (1 - q)*mu)/sigma2
  pars$phi <- (sqrt(c^2*(1 - h^2)^2 + 4*(c*h - g)*g*(1 - h^2)) - c*(1 - h^2))/
    (2*(c*h - g))
  pars$kappa <- h - pars$phi

  pars$psi <- ((sigma2 - (1 - q)*mu)*(1 - h^2) - q*mu*(1 - h^2 + pars$phi^2))/
    (pars$phi^2*(sigma2 - (1 - q)*mu) + mu^2*(1 - h^2 + pars$phi^2))

  pars$q <- q
  return(pars)
}

reparam <- function(nu, phi, kappa, psi, q){
  sop <- compute_sop(nu = nu, phi = phi, kappa = kappa, psi = psi, q = q)$X_tilde
  recover_pars(sop_list = sop)
}

#' Fitting time-homogeneous underreported endemic-epidemic model using maximum likelihood
#'
#' Fits a time-homogeneous endemic-epidemic model with underreporting using an approximative maximum
#' likelihood scheme. The likelihood approximation is based on an approximation of the process by
#' a second-order equivalent process with complete reporting. The reporting probability cannot be
#' estimated from the data (in most cases these contain no information on it) and thus needs to be
#' specified in advance.
#'
#' @param observed a time series of counts (numeric vector)
#' @param q the assumed reporting probability
#' @param initial the initial value of the parameter vector passed to optim
#' (note: the function tries different starting values in any case)
#' @param max_lag in evaluation of likelihood only lags up to max_lag are taken into account
#' @return the return object from \code{optim} providing the maximum likelihood estimates
#' on the log scale.
#' @export
fit_hhh4u <- function(observed, q, include_kappa = TRUE,
                    initial = c(log_nu = 2, log_phi = -1, log_kappa = -2, log_psi = -3),
                    max_lag = 5, iter_optim = 3, hessian = FALSE, ...){

  # remove kappa from initial values if desired:
  if(include_kappa == FALSE){
    initial <- initial[names(initial) != "alpha_kappa"]
  }

  nllik_vect <- function(pars){
    # extract parameter values:
    nu <- exp(pars["log_nu"])
    phi <- exp(pars["log_phi"])
    kappa <- ifelse(include_kappa, exp(pars["log_kappa"]), -10)
    psi <- exp(pars["log_psi"])
    nllik(observed = observed, nu = nu, phi = phi, kappa = kappa, psi = psi, q = q, max_lag = max_lag)
  }

  opt <- optim(par = initial, fn = nllik_vect, hessian = hessian, ...)
  # re-run optimization starting from previous optimum (improves convergence)
  for(i in 1:iter_optim){
    opt <- optim(par = opt$par, fn = nllik_vect, hessian = hessian, ...)
  }
  return(opt)
}

#' Evaluate the approximate log-likelihood of a time-homogeneous underreported endemic-epidemic model
#'
#' The likelihood approximation is based on an approximation of the process by
#' a second-order equivalent process with complete reporting.
#'
#' @param observed a time series of counts (numeric vector)
#' @param nu,phi,kappa,psi the model parameters (scalars)
#' @param q the assumed reporting probability
#' @param max_lag in evaluation of likelihood only lags up to max_lag are taken into account
#' @param return_contributions shall the log-likelihood contributions of each time point be
#' returned (as vector)?
#' @return The log-likelihood as scalar or (if \code{return_contributions == TRUE}) the vector of
#' log-likelihood contributions.
#'
#' @export
nllik <- function(observed, nu, phi, kappa, psi, q, max_lag = 5, return_contributions = FALSE){
  # get second order properties:
  sop <- compute_sop(nu = nu, phi = phi, kappa = kappa, psi = psi, q = q)$X_tilde

  if(any(unlist(sop) < 0)) return(-Inf)

  # get corresponding parameters for unthinned process:
  pars_Y <- recover_pars(q = 1, sop_list = sop)
  nu_Y <- pars_Y$nu
  phi_Y <- pars_Y$phi
  kappa_Y <- pars_Y$kappa
  psi_Y <- pars_Y$psi
  nu_Y_star <- nu_Y/(1 - kappa_Y) # move to geometric-lag display
  q <- 1
  # get likelihood:
  lgt <- length(observed)
  mod_matr <- matrix(nrow = length(observed), ncol = max_lag)
  for(i in 1: max_lag){
    mod_matr[, i] <- c(rep(NA, i), head(observed, lgt - i))
  }
  lin_pred <- nu_Y_star + mod_matr %*% matrix(phi_Y*kappa_Y^(1:max_lag - 1), ncol = 1)
  llik <- - dnbinom(observed, mu = as.vector(lin_pred), size = 1/psi_Y, log = TRUE)

  if(return_contributions) return(llik) else return(sum(llik[-(1:max_lag)]))
}

# helper function to obtain second-order properties of simulated time series
emp_sop <- function(observed){
  lgt <- length(observed)
  list(
    mean = mean(observed),
    var = var(observed),
    ar1 = cor(head(observed, lgt - 1), tail(observed, lgt - 1)),
    ar2 = cor(head(observed, lgt - 2), tail(observed, lgt - 2))
  )
}
