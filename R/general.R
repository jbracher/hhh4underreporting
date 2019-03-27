#' Simulate from underreported model with time-varying parameters
#'
#' Simulate from a time-homogeneous endemic-epidemic model with underreporting.
#' Can also be used through the wrapper \code{simulate_hhh4u}.
#'
#' @param lambda1 the initial value of lambda
#' @param nu,phi,kappa,psi the model parameters (vectors of equal length
#' corresponding to the length of the time series; scalars are recycled.)
#' @param q the reporting probability
#' @return A named list with elements \code{"X"}, \code{"X_tilde"} and \code{"lambda"} containing
#' the unthinned and thinned simulated time series as well as the conditional mean process.
#' @examples
#' sim <- simulate_hhh4u(lambda1 = 1, nu = c(1:5, 4:1), phi = 0.4, kappa = 0.2, psi = 0.1, q = 0.5)
#' @export
simulate_hhh4u <- function(lambda1, nu, phi, kappa, psi, q = 1,
                           return_sts = TRUE, args_sts = list()){
  lgt = length(nu)

  if(!length(phi) %in% c(1, lgt) |
     !length(kappa) %in% c(1, lgt) |
     !length(psi) %in% c(1, lgt)){
    stop("Parameter vectors nu, phi, kappa, psi all need to have the same length. Alternatively phi, kappa and psi can also be scalar (will be recycled).")
  }

  if(length(q) != 1) stop("Reporting probability q needs to be time-constant (just one scalar).")

  if(length(nu) == 1) nu <- rep(nu, lgt)
  if(length(phi) == 1) phi <- rep(phi, lgt)
  if(length(kappa) == 1) kappa <- rep(kappa, lgt)
  if(length(psi) == 1) psi <- rep(psi, lgt)

  lambda <- X <- Y <- numeric(lgt)
  lambda[1] <- lambda1
  X[1] <- rnbinom(1, mu = lambda[1], size = 1/psi[1])

  for(t in 2:lgt){
    lambda[t] <- nu[t] + phi[t]*X[t - 1] + kappa[t]*lambda[t - 1]
    X[t] <- rnbinom(1, mu = lambda[t], size = 1/psi[t])
  }

  X_tilde <- rbinom(n = lgt, size = X, prob = q)

  # return either sts object with observed counts or a named list with the
  # observed and unobserved processes
  if(return_sts){
    args_sts$observed <- X_tilde
    return(do.call(sts, args_sts))
  }else{
    return(list(X = X, X_tilde = X_tilde, lambda = lambda))
  }
}

#' Simulate from underreported model with time-varying parameters
#'
#' Simulate from a time-homogeneous endemic-epidemic model with underreporting.
#' Wrapper around \code{simulate_hhh4u0} which takes vector-valued parameters.
#'
#' @param lambda1 the initial value of lambda
#' @param beta_nu parameters for the linear form entering into \eqn{log(nu)}
#' @param beta_phi parameters for the linear form entering into \eqn{log(phi)}
#' @param kappa,psi Time-constant model parameters
#' @param q the reporting probability
#' @param z_nu Matrix containing the covariates entering into \eqn{log(nu)} (in columns)
#' @param z_phi Matrix containing the covariates entering into \eqn{log(phi)} (in columns)
#' @return A named list with elements \code{"X"}, \code{"X_tilde"} and \code{"lambda"} containing
#' the unthinned and thinned simulated time series as well as the conditional mean process.
#' @examples
#' beta_nu <- c(1, 0.005)
#' z_nu <- cbind(1, 1:100)
#' sim <- simulate_hhh4u2(lambda1 = 1, beta_nu = beta_nu, beta_phi = -0.6,
#' kappa = 0.2, psi = 0.1, q = 0.5, z_nu = z_nu)
#' @export
simulate_hhh4u2 <- function(lambda1, beta_nu, beta_phi, kappa, psi, q = 1, lgt,
                            z_nu = matrix(rep(1, lgt)),
                            z_phi = matrix(rep(1, lgt))){
  nu <- exp(z_nu%*%beta_nu)
  phi <- exp(z_phi%*%beta_phi)
  simulate_hhh4u(lambda1 = lambda1,
                 nu = nu, phi = phi, kappa = kappa,
                 psi = psi, q = q)
}

#' Marginal mean of an underreported model
#'
#' Compute the time-varying marginal means of an underreported model.
#'
#' @param lambda1 the initial value of lambda
#' @param nu,phi,kappa the time-varying parameters of the model (\code{kappa} is also allowed to be time-varying here);
#' vectors of equal length
#' @param q the reporting probability; defaults to 1.
#' @return a vector containing the marginal means
marg_mean_tv <- function(lambda1, nu, phi, kappa, q = 1){
  lgt <- length(nu)
  mu_X <- numeric(lgt)
  mu_X[1] <- lambda1
  for(i in 2:lgt){
    mu_X[i] <- nu[i] + (phi[i] + kappa[i])*mu_X[i - 1]
  }
  mu_X_tilde <- q*mu_X
  return(mu_X_tilde)
}

#' Marginal second-order properties of an underreported model
#'
#' Compute the time-varying mean and covariance structure of an underreported model.
#'
#' @param lambda1 the initial value of lambda
#' @param nu,phi,kappa,psi the time-varying parameters of the model
#'  (\code{kappa} and \code{psi} are also allowed to be time-varying here);
#' vectors of equal length
#' @param q the reporting probability
#' @param compute_Sigma logical: shall a full covariance matrix be computed?
#' @return a named list containing the means, variances, first-order autocovariances.
#' decay parameters of the autocovariance functions and, if desired, full autocovariance
#' matrix for both \eqn{X} and \eqn{\tilde{X}}.
compute_sop_tv <- function(lambda1, nu, phi, kappa, psi, q, compute_Sigma = FALSE){
  lgt <- length(nu)
  # get means:
  mu_X <- marg_mean_tv(lambda1 = lambda1, nu = nu,  phi = phi, kappa = kappa)
  # compute variances of lambda as they will be required:
  v_lambda <- v_X <- cov1_X <- numeric(lgt)
  v_lambda[1] <- 0
  v_X[1] <- lambda1 + psi[1]*lambda1^2

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

#' Approximate an underreported model by a fully observed one
#'
#' This function returns the parameters of a fully observed model with the same second-order properties
#' as the specified underreported model.
#'
#' @param nu,phi,kappa,psi the time-varying parameters of the model
#'  (\code{kappa} and \code{psi} are also allowed to be time-varying here);
#' vectors of equal length
#' @param q the reporting probability
reparam_tv <- function(lambda1, nu, phi, kappa, psi, q){
  lgt <- length(nu)
  # compute target second-order properties:
  target_sop <- compute_sop_tv(lambda1 = lambda1, nu = nu,
                               phi = phi, kappa = kappa,
                               psi = psi, q = q)

  # now find a completely observed process with the same properties:
  nu_Y <- q*nu # known for theoretical reasons (?)
  phi_plus_kappa_Y <- phi + kappa # this, too.

  mu_Y <- target_sop$mu_X_tilde
  phi_Y <- kappa_Y <- psi_Y <- v_lambda_Y <- numeric(lgt)

  v_lambda_Y[1] <- 0
  psi_Y[1] <- psi[1]

  v_Y <- target_sop$v_X_tilde # by definition
  cov1_Y <- target_sop$cov1_X_tilde # by definition
  phi_plus_kappa_Y <- target_sop$decay_cov_X_tilde

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

  return(list(lambda1 = mu_Y[1], nu = nu_Y, phi = phi_Y, kappa = kappa_Y,
              psi = psi_Y, q = 1))
}

#' Compute the fitted values of a fully observed process
#'
#' Compute the fitted values of a fully observed process, which can also serve as an
#' approximation to those of an underreported process which is approximated by the given
#' fully observed process.
#'
#' @param observed The observed sequence of counts
#' @param lambda1_Y initial value for \eqn{lambda}
#' @param nu_Y,phi_Y,kappa_Y the time-varying parameters of the model
#'  (\code{kappa} is also allowed to be time-varying here);
#' vectors of equal length
compute_fitted_values_tv <- function(observed, lambda1_Y, nu_Y, phi_Y, kappa_Y){
  lgt <- length(observed)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1_Y
  for(t in 2:lgt){
    lambda[t] <- nu_Y[t] + phi_Y[t]*observed[t - 1] + kappa_Y[t]*lambda[t - 1]
  }
  return(lambda)
}

#' Evaluate the (approximate) likelihood of an underreported model
#'
#' Compute the fitted values of a fully observed process, which can also serve as an
#' approximation to those of an underreported process which is approximated by the given
#' fully observed process.
#'
#' @param observed The observed sequence of counts
#' @param lambda1 initial value for \eqn{lambda}
#' @param nu,phi,kappa,psi the time-varying parameters of the model
#'  (\code{kappa} is also allowed to be time-varying here);
#' vectors of equal length
#' @param q the assumed reporting probability
#' @param return_contributions logical: should the likelihood contributions
#' of the different weeks be returned?
#' @param log logical: return log-likelihood?
llik_tv <- function(observed, lambda1, nu, phi, kappa, psi, q,
                    return_contributions = FALSE, log = TRUE){
  lgt <- length(observed)

  # get corresponding parameters for unthinned process:
  pars_Y <- reparam_tv(lambda1 = lambda1, nu = nu, phi = phi,
                       kappa = kappa, psi = psi, q = q)

  lambda <- compute_fitted_values_tv(observed = observed, lambda1_Y = pars_Y$lambda1,
                                     nu = pars_Y$nu, phi = pars_Y$phi, kappa = pars_Y$kappa)


  llik <- dnbinom(observed, mu = lambda, size = 1/pars_Y$psi,
                  log = log)
  if(return_contributions){
    return(llik)
  }else{
    if(log) return(sum(llik)) else return(prod(llik))
  }
}

#' A slower version of hhh4u using only R code
hhh4u_R <- function(stsObj,
                    control = list(ar = list(f = ~ -1, use_kappa = TRUE),
                                   end = list(f = ~ 1),
                                   family = c("Poisson", "NegBin1"),
                                   q = NULL,
                                   subset = 1:nrow(stsObj),
                                   data = list(t = stsObj@epoch - min(stsObj@epoch)),
                                   start = NULL,
                                   optimizer = list(maxit = 1000),
                                   return_se = TRUE)
){

  control <- setControl(control = control, stsObj = stsObj)
  observed <- stsObj@observed[control$subset]

  nllik_vect <- function(pars){
    lgt <- length(observed)
    lambda1 = exp(pars["log_lambda1"])

    beta_nu <- pars[grepl("end.", names(pars))]
    beta_phi <- pars[grepl("ar.", names(pars))]

    nu <- exp(control$matr_nu %*% beta_nu)
    phi <- exp(control$matr_phi %*% beta_phi)
    kappa <- rep(exp(pars["log_kappa"]), lgt)
    psi <- if(control$family[1] == "NegBin1"){
      rep(exp(pars["log_psi"]), lgt)
    }else{
      rep(0.0001, lgt)
    }
    q <- control$q

    -llik_tv(observed = observed, lambda1 = lambda1, nu = nu,
             phi = phi, kappa = kappa, psi = psi,
             q = q)
  }

  opt <- optim(par = control$start, fn = nllik_vect, hessian = control$return_se,
               control = control$optimizer)
  for(i in 1:3){
    opt <- optim(par = opt$par, fn = nllik_vect, control = control$optimizer,
                 hessian = control$return_se)
  }

  # extract time-varying parameter values:
  lambda1 = exp(opt$par["log_lambda1"])
  beta_nu <- opt$par[grepl("end.", names(opt$par))]
  beta_phi <- opt$par[grepl("ar.", names(opt$par))]
  nu <- exp(control$matr_nu %*% beta_nu)
  phi <- exp(control$matr_phi %*% beta_phi)
  kappa <- rep(exp(opt$par["log_kappa"]), length(observed))
  psi <- rep(exp(opt$par["log_psi"]), length(observed))
  q <- control$q
  # get corresponding parameters for unthinned process:
  pars_Y <- reparam_tv(lambda1 = lambda1, nu = nu, phi = phi,
                       kappa = kappa, psi = psi, q = q)
  # compute fitted values:
  lambda <- compute_fitted_values_tv(observed = observed, lambda1_Y = pars_Y$lambda1,
                                     nu = pars_Y$nu, phi = pars_Y$phi, kappa = pars_Y$kappa)

  # create an informative return object inspired by hhh4:
  ret <- list()

  # parameter estimates and standard errors:
  ret$coefficients <- opt$par
  ret$se <- ret$cov <- NULL
  if(control$return_se){
    ret$cov <- solve(opt$hessian)
    ret$se <- diag(ret$cov)
  }
  ret$par_long <- list(lambda1 = lambda1, nu = nu, phi = phi,
                       kappa = kappa, psi = psi, q = q)
  ret$par_long_approximation <- pars_Y

  # fitted values:
  ret$fitted.values <- lambda

  # other:
  ret$dim <- length(ret$par)
  ret$loglikelihood <- -opt$value
  ret$convergence <- (opt$convergence == 0)
  ret$control <- control
  ret$stsObj <- stsObj
  ret$nobs <- length(control$subset)
  ret$optim <- opt

  class(ret) <- "hhh4u"
  return(ret)
}


#' Fit an underreported model to an observed time series
#'
#' The reporting probability needs to be specified!
#'
#' @param stsObj An sts object containing the observed time series. Only univariate time series are supported.
#' @param control a control list, similar to the one used in \code{hhh4} from the \code{surveillance} package.
#' The following aspects can be specified:
#' \itemize{
#' \item \code{ar}: List specifying the autoregressive component. Possible elements:
#' \itemize{
#' \item \code{f}: a formula
#' \item \code{use_kappa}: should a feedback mechanism be used?
#' }
#' \item \code{end}: List specifying the endemic component. Possible elements:
#' \itemize{
#' \item \code{f}: a formula
#' }
#' \item \code{family} The distibutional family, either \code{"Poisson"} or \code{"NegBin1"} for a negative binomial.
#' \item \code{q} The assumed reporting probability
#' \item \code{subset} The subset of the data to which the model shall be fitted; needs to contained within \code{1:nrow(stsObj@observed) }
#' \item \code{start} Initial value passed to \code{optim}; needs to have the correct length and naming.
#' \item \code{optimize} Additional arguments passed to \code{optim}
#' \item \code{return_se} Should the standard errors be returned? In case of numerical problems it can be reasonable to deactivate this.
#' }
#' @return An object of class \code{hhh4u}, which is a named list with the following elements:
#' \itemize{
#' \item \code{coefficients}: The parameter estimates
#' \item \code{se}: The estimated standard errors of the coefficients
#' \item \code{cov}: The estimated autocovariance matrix of the parameter estimates
#' \item \code{par_long}: a list containing the time-varying parameters in vector form
#' \item \code{par_long_approximation}: a list containing the time-varying parameters of the approximating fully observed process
#' \item \code{fitted.values}: the fitted values
#' \item \code{dim}: the model dimension (i.e. number of parameters)
#' \item \code{loglikelihood:} the log-likelihood of the model
#' \item \code{convergence}: code for convergence status returned by \code{optim}; 0 indicates successful convergence.
#' \item \code{control}: the control list used for the fit
#' \item \code{stsObj}: the \code{sts} object used for the fit
#' \item \code{nobs}: the number of observations used for the fit
#' \item \code{optim}: the return object of the call to optim
#' }
#' @export
hhh4u <- function(stsObj,
                  control = list(ar = list(f = ~ 1, use_kappa = TRUE),
                                 end = list(f = ~ 1),
                                 family = c("Poisson", "NegBin1"),
                                 q = NULL,
                                 subset = 1:nrow(stsObj),
                                 data = list(t = stsObj@epoch - min(stsObj@epoch)),
                                 start = NULL,
                                 optimizer = list(maxit = 5000),
                                 return_se = TRUE)
){

  control <- setControl(control = control, stsObj = stsObj)
  observed <- stsObj@observed[control$subset]

  nllik_vect <- function(pars){
    lgt <- length(observed)
    lambda1 = exp(pars["log_lambda1"])

    beta_nu <- pars[grepl("end.", names(pars))]
    beta_phi <- pars[grepl("ar.", names(pars))]

    nu <- exp(control$matr_nu %*% beta_nu)
    phi <- exp(control$matr_phi %*% beta_phi)
    kappa <- rep(exp(pars["log_kappa"]), lgt)
    psi <- if(control$family[1] == "NegBin1"){
      rep(exp(pars["log_psi"]), lgt)
    }else{
      rep(0.0001, lgt)
    }
    q <- control$q

    nllik_tv_cpp(observed = observed, lambda1 = lambda1, nu = nu,
                 phi = phi, kappa = kappa, psi = psi,
                 q = q)
  }

  opt <- optim(par = control$start, fn = nllik_vect, hessian = control$return_se,
               control = control$optimizer)
  for(i in 1:3){
    opt <- optim(par = opt$par, fn = nllik_vect, control = control$optimizer,
                 hessian = control$return_se)
  }

  # extract time-varying parameter values:
  lambda1 = exp(opt$par["log_lambda1"])
  beta_nu <- opt$par[grepl("end.", names(opt$par))]
  beta_phi <- opt$par[grepl("ar.", names(opt$par))]
  nu <- exp(control$matr_nu %*% beta_nu)
  phi <- exp(control$matr_phi %*% beta_phi)
  kappa <- rep(exp(opt$par["log_kappa"]), length(observed))
  psi <- rep(exp(opt$par["log_psi"]), length(observed))
  q <- control$q
  # get corresponding parameters for unthinned process:
  pars_Y <- reparam_tv(lambda1 = lambda1, nu = nu, phi = phi,
                       kappa = kappa, psi = psi, q = q)
  # compute fitted values:
  lambda <- compute_fitted_values_tv(observed = observed, lambda1_Y = pars_Y$lambda1,
                                     nu = pars_Y$nu, phi = pars_Y$phi, kappa = pars_Y$kappa)

  # create an informative return object inspired by hhh4:
  ret <- list()

  # parameter estimates and standard errors:
  ret$coefficients <- opt$par
  ret$q <- control$q
  ret$se <- ret$cov <- NULL
  if(control$return_se){
    ret$cov <- solve(opt$hessian)
    ret$se <- diag(ret$cov)
  }
  ret$par_long <- list(lambda1 = lambda1, nu = nu, phi = phi,
                       kappa = kappa, psi = psi, q = q)
  ret$par_long_approximation <- pars_Y

  # fitted values:
  ret$fitted.values <- lambda
  ret$observed <- observed
  ret$timepoints <- seq(from = stsObj@start[1] + stsObj@start[2]/stsObj@freq,
                        length.out = nrow(stsObj@observed), by = 1/stsObj@freq)[control$subset]

  # other:
  ret$dim <- length(ret$par)
  ret$loglikelihood <- -opt$value
  ret$convergence <- (opt$convergence == 0)
  ret$control <- control
  ret$stsObj <- stsObj
  ret$nobs <- length(control$subset)
  ret$nTime <- length(control$subset)
  ret$optim <- opt
  if(!ret$convergence) warning("Optimizer did not converge.")
  ret$call <- match.call()

  class(ret) <- "hhh4u"
  return(ret)
}

#' Check and complete the control list
#'
#' Similar to \code{surveillance:::setControl and surveillance:::makeControl} this function
#' checks the arguments in the control list and adds default settings where necessary.
#' @param control the control argument as provided by the user
#' @param stsObj the stsObj to be used for the fit
#' @return a complete control list, with default values added where nothing was specified.
setControl <- function(control, stsObj){
  stopifnot(is.list(control))
  nTime <- nrow(stsObj)
  nUnit <- ncol(stsObj)
  if(nUnit > 1) stop("hhh4u is only available for univariate analyses.")
  if (nTime <= 2)
    stop("too few observations")
  defaultControl <- eval(formals(hhh4u)$control)
  control <- modifyList(defaultControl, control)

  # check provided arguments:
  if(!is.list(control$ar)){
    stop("control$ar must be a list.")
  }
  if(!is.list(control$end)){
    stop("control$end must be a list.")
  }
  if(is.null(control$q)){
    stop("A value for q needs to be provided in the control argument.")
  }
  if(!inherits(control$end$f, "formula")){
    stop("'control$end$f' must be a formula")
  }
  if(!inherits(control$ar$f, "formula")){
    stop("'control$ar$f' must be a formula")
  }
  if(control$ar$f == ~-1 | control$end$f == ~-1){
    stop("The formula arguments must not be ~-1 (both components must be included).")
  }
  if(!(control$family[1] %in% c("Poisson", "NegBin1"))){
    stop("control$family needs to be either 'Poisson' or 'NegBin1'.")
  }
  if(!is.logical(control$ar$use_kappa)){
    stop("control$ar$use_kappa needs to be logical.")
  }
  if(!is.numeric(control$subset) |
     any(diff(control$subset) != 1) |
     any(!(control$subset %in% seq_along(stsObj@observed)))){
    stop("control$subset needs to be a vector of consecutive integers between 1 and nrow(stsObj@observed).")
  }
  if(!is.logical(control$return_se)){
    stop("return_se needs to be logical.")
  }

  # check data:
  if(!is.list(control$data)){
    stop("control$data must be a named list of covariates")
  }
  for(i in seq_along(control$data)){
    if(length(control$data[[i]]) != nrow(stsObj@observed)){
      stop("data needs to contain covariate vectors of same length as stsObj@observed.")
    }
  }
  # construct the model matrices
  data.fr <- as.data.frame(control$data)[control$subset, , drop = FALSE]
  control$matr_nu <- model.matrix(control$end$f, data = data.fr)
  control$matr_phi <- model.matrix(control$ar$f, data = data.fr)

  # address starting values:
  start_lambda1 <- c("log_lambda1" = 1)
  start_end <- rep(0, ncol(control$matr_nu))
  names(start_end) <- paste0("end.", colnames(control$matr_nu))
  start_ar <- rep(0, ncol(control$matr_phi))
  names(start_ar) <- paste0("ar.", colnames(control$matr_phi))
  start_ar["ar.(Intercept)"] <- -1
  start_kappa <- if(control$ar$use_kappa) c(log_kappa = -1) else NULL
  start_psi <- if(control$family[1] == "NegBin1") c(log_psi = -1) else NULL
  default_start <- c(start_lambda1, start_end, start_ar,
                     start_kappa, start_psi)
  if(is.null(control$start)){
    control$start <- default_start
  }else{
    if(length(control$start) != length(default_start)) stop("control$start is not of the right length.")
  }

  return(control)
}

#' Method to return fitted values
#' @export
fitted.hhh4u <- function(fit){
  fit$fitted.values
}

#' Method to return observed values
#' @export
observed.hhh4u <- function(fit){
  fit$observed
}

#' Method for a simple plot
#' @export
plot.hhh4u <- function(fit, type  = c("fitted"), ...){

  if(type == "fitted"){
    plot(fit$timepoints, fit$fitted.values, type = "l", ylim = c(0, max(fit$observed)),
         xlab = "Time", ylab = "observed cases", ...)
    points(fit$timepoints, fit$observed, cex = 0.6)
  }
}

#' Method to obtain summary
#' @export
summary.hhh4u <- function(fit){
  ret <- fit[c("call", "convergence", "dim", "loglikelihood",
                    "nTime", "coefficients", "se")]
  class(ret) <- "summary.hhh4u"
  return(ret)
}

#' @export
logLik.hhh4u <- function(fit){
  fit$loglikelihood
}

#' @export
AIC.hhh4u <- function(fit){
  2*(fit$dim - fit$loglikelihood)
}

#' Methods to get BIC
#' @export
BIC.hhh4u <- function(fit){
  2*(fit$dim*log(fit$nTime) - fit$loglikelihood)
}

#' Printing method for hhh4u
#' @export
print.hhh4u <- function(object){
  if (!object$convergence) {
    cat("Note: optimizer did not converge.\n")
  }
  cat("Underreported endemic-epidemic model fitted under the assumption q =", object$q, "\n")
  cat("Coefficients:\n")
  print(object$coefficients)
  cat("Standard errors:\n")
  print(object$se)
}

#' Method to print summary
#' @export
print.summary.hhh4u <- function(object){
  print.hhh4u(object)
  cat("\nLog-likelihood:", round(object$loglikelihood, 3), "\n")
  cat("Number of time points:", object$nTime, "\n")
  cat("AIC:", round(AIC.hhh4u(object), 3), "\n")
  cat("BIC:", round(BIC.hhh4u(object), 3), "\n")
}

#' Simulation method
#' @export
simulate.hhh4u <- function(fit){
  do.call(simulate_hhh4u, fit$par_long)
}
