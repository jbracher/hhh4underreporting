# function to simulate from time-homogeneous model (just a wrapper around the seasonal case):
generate_ar <- function(nu, phi, kappa, psi, lgt = 100,
                             start = 10, burn_in = 1000){
  generate_ar_seas(nu = nu, phi = phi, kappa = kappa, psi = psi,
                   n_seas = lgt, start = start, burn_in = 1000)
}

# function to recover second-order properties
compute_sop <- function(nu, phi, kappa, psi, p, par_list = NULL){
  if(!is.null(par_list)){
    nu <- par_list$nu; phi <- par_list$phi; kappa <- par_list$kappa; psi <- par_list$psi; p <- par_list$p
    names(nu) <- names(phi) <- names(kappa) <- names(psi) <- names(p) <- NULL
  }

  soc <- list()
  unthinned <- thinned <- list()

  unthinned$mu <- nu/(1 - phi - kappa)
  unthinned$sigma2 <- (1 - (phi + kappa)^2 + phi^2)/
    (1 - (phi + kappa)^2 - psi*phi^2) *(unthinned$mu + unthinned$mu^2*psi)
  unthinned$g <- phi*(1 - kappa*(phi + kappa))/(1 - (phi + kappa)^2 + phi^2)
  unthinned$h <- phi + kappa

  thinned$mu <- p*unthinned$mu
  thinned$sigma2 <- p^2*unthinned$sigma2 + p*(1 - p)*unthinned$mu
  thinned$g <- unthinned$sigma2/(unthinned$sigma2 + (1 - p)/p*unthinned$mu)*unthinned$g
  thinned$h <- phi + kappa

  return(list(unthinned = unthinned,
              thinned = thinned))
}

# function to recover nu, pi, kappa, psi for a given p:
recover_pars <- function(mu, sigma2, g, h, p = 1, sop_list = NULL){
  if(!is.null(sop_list)){
    mu <- sop_list$mu; sigma2 <- sop_list$sigma2; g <- sop_list$g; h <- sop_list$h
  }

  pars <- NULL
  pars$nu <- mu*(1 - h)/p
  gamma <- (sigma2 - (1 - p)*mu)/sigma2
  pars$phi <- (gamma*(1 - h^2) - sqrt(gamma^2*(1 - h^2)^2 - 4*(g - gamma*h)*g*(1 - h^2)))/
    (2*(g - gamma*h))
  pars$kappa <- h - pars$phi

  sigma2_star <- (sigma2 - (1 - p)*mu)/p^2
  mu_star <- mu/p

  pars$psi <- (sigma2_star*(1 - h^2) - mu_star*(1 - h^2 + pars$phi^2))/
    (pars$phi^2*sigma2_star + mu_star^2*(1 - h^2 + pars$phi^2))

  pars$p <- p
  return(pars)
}

reparam <- function(nu, phi, kappa, psi, p){
  sop <- compute_sop(nu = nu, phi = phi, kappa = kappa, psi = psi, p = p)$thinned
  recover_pars(sop_list = sop)
}

# function to maximize likelihood:
fit_lik <- function(Y, p, initial = c(log_nu = 4, log_phi = -2, log_kappa = -3, log_psi = 3), max_lag = 5, hessian = FALSE, ...){
  lik_vect <- function(pars){
    # extract parameter values:
    nu <- exp(pars["log_nu"])
    phi <- exp(pars["log_phi"])
    kappa <- exp(pars["log_kappa"])
    psi <- exp(pars["log_psi"])
    lik(Y = Y, nu = nu, phi = phi, kappa = kappa, psi = psi, p = p, max_lag = max_lag)
  }
  optim(par = initial, fn = lik_vect, hessian = hessian, ...)
}

lik <- function(Y, nu, phi, kappa, psi, p, max_lag = 5){
  # get second order properties:
  sop <- compute_sop(nu = nu, phi = phi, kappa = kappa, psi = psi, p = p)$thinned

  if(any(unlist(sop) < 0)) return(-Inf)

  # get corresponding parameters for unthinned process:
  pars_star <- recover_pars(p = 1, sop_list = sop)
  nu_star <- pars_star$nu
  phi_star <- pars_star$phi
  kappa_star <- pars_star$kappa
  psi_star <- pars_star$psi
  nu_star_tilde <- nu_star/(1 - kappa_star) # move to geometric-lag display
  p <- 1
  # get likelihood:
  lgt <- length(Y)
  mod_matr <- matrix(nrow = length(Y), ncol = max_lag)
  for(i in 1: max_lag){
    mod_matr[, i] <- c(rep(NA, i), head(Y, lgt - i))
  }
  lin_pred <- nu_star_tilde + mod_matr %*% matrix(phi_star*kappa_star^(1:max_lag - 1), ncol = 1)
  llik <- - sum(dnbinom(Y, mu = as.vector(lin_pred), size = 1/psi, log = TRUE), na.rm = TRUE)
  return(llik)
}

emp_sop <- function(Y){
  lgt <- length(Y)
  list(
    Y_mean_emp = mean(Y),
    Y_var_emp = var(Y),
    Y_ar1_emp = cor(head(Y, lgt - 1), tail(Y, lgt - 1)),
    Y_ar2_emp = cor(head(Y, lgt - 2), tail(Y, lgt - 2))
  )
}
