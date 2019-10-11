# Demonstrating that debiasing estimates from a fit with reporting probability 1
# is not a good idea

# Johannes Bracher, johannes.bracher@uzh.ch


library(hhh4underreporting)

setwd("/home/johannes/Documents/underreporting/Article_hhh4/simulations/")

n_sim <- 1000 # number of simulation iterations
lgt <- 8*52 # length of time series

# set true parameters:
nu <- 20
phi <- 0.4
kappa <- 0.2
psi <- 0.1
q <- 0.2

# function to choose starting values based on moment estimators:
get_starting_values <- function(sts, q = 1){
  # extract counts:
  vect <- as.vector(sts@observed)
  # get empirical second-order properties:
  ac <- acf(vect, lag.max = 2, plot = FALSE)$acf[2:4]
  emp_sop <- list(mu = mean(vect), sigma2 = var(vect),
                  g = ac[1], h = min(ac[2]/ac[1], 0.9))
  # transform to parameter starting values:
  start0 <- hhh4underreporting:::recover_pars_homog(sop_list = emp_sop, q = q)

  # catch cases where moment estimators ill-defined:
  if(is.nan(start0$nu) | start0$nu < 0) start0$nu <- 5
  if(is.nan(start0$phi) | start0$phi < 0) start0$phi <- max(ac[1], 0.1)
  if(is.nan(start0$kappa) | start0$kappa < 0) start0$kappa <- 0.05
  if(is.nan(start0$psi) | start0$psi < 0) start0$psi <- 0.01

  # format so that hhh4u can use them:
  start <- c("log_lambda1" = log(vect[1] + 1),
             "end.(Intercept)" = log(start0$nu),
             "ar.(Intercept)" = log(start0$phi),
             "log_kappa" = log(start0$kappa),
             "log_psi" = log(start0$psi))

  return(start)
}

# matrices to store results:
estim <- estim_corr <- estim_naive <-
  matrix(ncol = 5, nrow = n_sim,
         dimnames = list(NULL, c("lambda1", "nu", "phi", "kappa", "psi")))


for(i in 1:n_sim){
  set.seed(i)

  # generate time series
  sts_temp <- sts(simulate_hhh4u_homog(nu, phi, kappa, psi, q, lgt = 8*52)$X_tilde)

  # starting values for fitting with correct q
  start_q <- get_starting_values(sts = sts_temp, q = q)
  # fit with correct q:
  fit_q <- hhh4u(sts_temp, control = list(q = q,
                                          family = "NegBin1",
                                          start = start_q,
                                          return_se = FALSE))
  estim[i, ] <- exp(fit_q$coefficients)

  # starting values for naive fitting
  start_1 <- get_starting_values(sts = sts_temp, q = 1)
  # fit with q = 1
  fit_1 <- hhh4u(sts_temp, control = list(q = 1,
                                          family = "NegBin1",
                                          start = start_1,
                                          return_se = FALSE))
  estim_naive[i, ] <- exp(fit_1$coefficients)

  # correct estimates (de-biasing)
  par_list_naive_temp <- as.list(estim_naive[i, -1])
  estim_corr[i, ] <-
    unlist(hhh4underreporting:::recover_pars_homog(
      q = q, # correct reporting probability
      sop_list = hhh4underreporting:::compute_sop_homog( # second-order properties of fitted model
        par_list = par_list_naive_temp, q = 1)$X
    )
    )[c(NA, "nu", "phi", "kappa", "psi")]
  print(i)
}

results <- list(estim = estim, estim_naive = estim_naive, estim_corr = estim_corr)

save(nu, phi, kappa, psi, q, lgt, n_sim,
     results,
     file = paste0("approximation_vs_correction/results_approximation_vs_correction_", lgt, "_", n_sim, ".rda"))

# nice plot for manuscript:

bp <- function(vect, at, shift = 0, col = "white"){
  boxplot(vect, at = at + shift, add = TRUE, boxwex = 0.1, col = col, pch = 16, lty = 1,
          staplewex = 0, medlwd = 1, cex = 0.5)
}

all_bps <- function(samples, variable, vals_q, col = "white",
                    true_value = NA, expected_bias = rep(NA, length(vals_q)),
                    xlim = 0:1, ylim = 0:1,
                    xlab = expression(true~pi), ...){

  if(length(col) == 1) col <- rep(col, length(vals_q))

  plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, axes = FALSE, ...)
  axis(1, at = c(0, 0.1, 0.25, 0.5, 0.75, 1))
  abline(h = true_value, col = "chartreuse3")
  lines(seq_along(expected_bias)/length(expected_bias), expected_bias, col = "red")
  for(i in seq_along(vals_q)){
    bp(samples[[i]][, variable], at = vals_q[i], col = col[i])
  }
}

layout(matrix(c(1, 2, 3, 4,
                1, 2, 3, 4,
                1, 2, 3, 4,
                5, 5, 5, 5), ncol = 4, byrow = TRUE))

par(las = 1, mar = c(1, 4.5, 1, 1), font.main = 1, family = "serif")

# my method:
all_bps(samples = results, variable = "nu", vals_q = 1:3/4, xlim = c(0, 1), ylim = c(0, 30),
        ylab = expression(hat(nu)), xlab = "", col = c("white", "red", "yellow"),
        true_value = nu)
abline(h = 0)

all_bps(samples = results, variable = "phi", vals_q = 1:3/4, xlim = c(0, 1), ylim = c(-1, 1),
        ylab = expression(hat(phi)), xlab = "", col = c("white", "red", "yellow"),
        true_value = phi)
abline(h = 0)

all_bps(samples = results, variable = "kappa", vals_q = 1:3/4, xlim = c(0, 1), ylim = c(-1, 1),
        ylab = expression(hat(kappa)), xlab = "", col = c("white", "red", "yellow"),
        true_value = kappa)
abline(h = 0)

all_bps(samples = results, variable = "psi", vals_q = 1:3/4, xlim = c(0, 1), ylim = c(0, 0.2),
        ylab = expression(hat(psi)), xlab = "", col = c("white", "red", "yellow"),
        true_value = psi)
abline(h = 0)

par(las = 1, mar = c(0, 4.5, 1, 1), font.main = 1, family = "serif")
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes =  FALSE)
legend("top", legend = c("correction during\n optimization", "no correction", "correction after\n estimation"),
       fill = c("white", "red", "yellow"), col = "black", ncol = 3, bty = "n")

