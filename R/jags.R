fit_jags <- function(Y, p, n.chains = 5, n.iter = 10000,
                     n.burnin = 1000, n.thin = 10, ...){
  # model definition
  fct_jags <- function ()
  {
    nu_geo ~ dexp(0.001)
    phi ~ dexp(0.01)
    kappa ~ dexp(0.01)
    psi ~ dexp(0.001)
    nu = nu_geo*(1 - kappa)

    for(t in 1:5){
      lambda[t] ~ dexp(0.001)
      omega[t] ~ dgamma(1/psi, 1/psi)
      X[t] ~ dpois(omega[t]*lambda[t])
      Y[t] ~ dbinom(p, X[t])
    }

    for (t in 6:n_timepoints) {
      lambda[t] = nu_geo + phi*X[t - 1] + phi*kappa*X[t - 2] + phi*kappa^2*X[t - 3] +
        phi*kappa^3*X[t - 4] + phi*kappa^4*X[t - 5]
      omega[t] ~ dgamma(1/psi, 1/psi)
      X[t] ~ dpois(omega[t]*lambda[t])
      Y[t] ~ dbinom(p, X[t])
    }
  }

  pars_to_save <-  c("nu", "phi", "kappa", "psi", "p")
  dat <- list(n_timepoints = length(Y), Y = Y, p = p)
  inits0 <- list(list(X = round(Y/p), nu_geo = 10, phi = 0.4, kappa = 0.2, psi = 0.1))
  inits <- rep(inits0, n.chains)

  R2jags::jags(data = dat,
               parameters.to.save = pars_to_save,
               inits = inits,
               n.chains = n.chains, n.iter = n.iter,
               n.burnin = n.burnin, n.thin = n.thin,
               model.file = fct_jags, ...)
}
