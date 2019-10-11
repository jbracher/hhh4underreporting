# exploring the biases occurring  when underreporting is ignored or
# multiplication factors are used.

# Johannes Bracher, johannes.bracher@uzh.ch

library(hhh4underreporting)

setwd("/home/johannes/Documents/underreporting/Article_hhh4/simulations/")

# number of simulation iterations:
n_sim <- 1000

# true parameters:
nu <- 15
phi <- 0.4
kappa <- 0.3
psi <- 0.1
vals_q <- c(1, 2.5, 5, 7.5, 10)/10 # five different values for reporting probability
lgt <- 8*52 # set length of simulated time series
pars <- c("nu" = nu, "phi" = phi, "kappa" = kappa, "psi" = psi)

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

# prepare lists of matrices to store results:
estim0 <- matrix(ncol = 5, nrow = n_sim)
colnames(estim0) <- c("lambda1", "nu", "phi", "kappa", "psi")
estim <- Y <- list()
for(i in seq_along(vals_q)){
  estim[[i]] <- estim0
  Y[[i]] <- matrix(ncol = lgt, nrow = n_sim)
}
names(estim) <- names(Y) <- paste0("p", vals_q)
estim_naive <- estim_stretch <- estim

# set seeds. We set the seed in each iteration to be able to get back to certain settings
# in case convergence issues occur.
set.seed(123)
seeds <- sample(1:10000, n_sim)

# run:
for(i in 1:n_sim){

  set.seed(seeds[i]) # set seed

  # simulate latent proces X:
  X_temp <- simulate_hhh4u_homog(nu = nu, phi = phi, kappa = kappa, psi = psi, q = 1, lgt = lgt)$X

  # run over different reporting probabilities:
  for(j in seq_along(vals_q)){
    # generate underreported counts:
    q_temp <- vals_q[j]
    Y_temp <- rbinom(lgt, X_temp, q_temp)
    Y[[j]][i, ] <- Y_temp
    sts_temp <- sts(observed = Y_temp)

    # fit model using our method:
    start_fl <- get_starting_values(sts = sts_temp, q = q_temp)
    fl <- hhh4u(sts_temp, control = list(q = q_temp,
                                         family = "NegBin1",
                                         start = start_fl,
                                         return_se = FALSE))
    estim[[j]][i, ] <- exp(fl$coefficients) # keep parameter estimates

    # fit model ignoring underreporting:
    start_fl_naive <- get_starting_values(sts = sts_temp, q = 1)
    fl_naive <- hhh4u(sts_temp, control = list(q = 1,
                                               family = "NegBin1",
                                               start = start_fl_naive,
                                               return_se =  FALSE))
    estim_naive[[j]][i, ] <- exp(fl_naive$coefficients)

    # fit model using stretching approach:
    sts_stretched_temp <- sts(observed = round(Y_temp/q_temp))
    start_fl_stretched <- get_starting_values(sts = sts_stretched_temp, q = 1)
    fl_stretch <- hhh4u(sts_stretched_temp, control = list(q = 1,
                                                 family = "NegBin1",
                                                 start = start_fl_stretched,
                                                 return_se = FALSE))
    estim_stretch[[j]][i, ] <- exp(fl_stretch$coefficients)
  }
  print(i)
}

# save:
save(nu, phi, kappa, psi, vals_q, lgt, seeds, estim, estim_naive, estim_stretch,
     file = paste0("bias_variance/results_bias_variance_", lgt, "_", n_sim, "_general.rda"))
# load(paste0("bias_variance/results_bias_variance_", lgt, "_", n_sim, "_general.rda"))

# store as csv:
estim_frame <- data.frame(true_nu = rep(nu, n_sim),
                          true_phi = rep(phi, n_sim),
                          true_kappa = rep(kappa, n_sim),
                          true_psi = rep(psi, n_sim))

for(p in names(estim)){
  colnames(estim[[p]]) <- paste0("estim_", colnames(estim[[p]]), "_", p)
  colnames(estim_naive[[p]]) <- paste0("estim_naive_", colnames(estim_naive[[p]]), "_", p)
  colnames(estim_stretch[[p]]) <- paste0("estim_stretch_", colnames(estim_stretch[[p]]), "_", p)

  estim_frame <- cbind(estim_frame, estim[[p]], estim_naive[[p]], estim_stretch[[p]])
}

colnames(estim_frame)

write.csv(estim_frame, file = paste0("bias_variance/results_bias_variance_", lgt, "_", n_sim, "_general.csv"))


# get expectations for naive estimates:
vals_q_long <- 1:100/100
expec_naive <- matrix(nrow = 4, ncol = length(vals_q_long))
rownames(expec_naive) <- c("nu", "phi", "kappa", "psi")
colnames(expec_naive) <- vals_q_long
for(k in seq_along(vals_q_long)){
  true_sop <- hhh4underreporting:::compute_sop_homog(nu, phi, kappa, psi, vals_q_long[k])$X_tilde
  expec_naive[, k] <- unlist(hhh4underreporting:::recover_pars_homog(q = 1, sop_list = true_sop))[1:4]
}

# get expectations for estimates based on stretching with underreporting factors:
expec_stretched <- NA*expec_naive
for(k in seq_along(vals_q_long)){
  true_sop <- hhh4underreporting:::compute_sop_homog(nu, phi, kappa, psi, vals_q_long[k])$X_tilde
  stretched_sop <- true_sop
  stretched_sop$mu <- true_sop$mu/vals_q_long[k]
  stretched_sop$sigma2 <- true_sop$sigma2/vals_q_long[k]^2

  expec_stretched[, k] <- unlist(hhh4underreporting:::recover_pars_homog(q = 1, sop_list = stretched_sop))[1:4]
}

# add estimates of effective reproduction numbers:
for(i in seq_along(vals_q)){
  estim[[i]] <- cbind(estim[[i]],
                      Reff = estim[[i]][, "phi"]/(1 - estim[[i]][, "kappa"]))
  estim_naive[[i]] <- cbind(estim_naive[[i]],
                            Reff = estim_naive[[i]][, "phi"]/
                              (1 - estim_naive[[i]][, "kappa"]))
  estim_stretch[[i]] <- cbind(estim_stretch[[i]],
                              Reff = estim_stretch[[i]][, "phi"]/
                                (1 - estim_stretch[[i]][, "kappa"]))
}


# helper functions for small boxplots:
# plot one small boxplot
bp <- function(vect, at, shift = 0, col = "white"){
  boxplot(vect, at = at + shift, add = TRUE, boxwex = 0.1, col = col, pch = 16, lty = 1,
          staplewex = 0, medlwd = 1, cex = 0.5)
}
# plot several small boxplots
all_bps <- function(samples, variable, vals_q, col = "white",
                    true_value = NA, expected_bias = rep(NA, length(vals_q)),
                    xlim = 0:1, ylim = 0:1,
                    xlab = expression(true~pi), ...){
  plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, axes = FALSE, ...)
  axis(1, at = c(0, 0.1, 0.25, 0.5, 0.75, 1))
  abline(h = true_value, col = "chartreuse3")
  lines(seq_along(expected_bias)/length(expected_bias), expected_bias, col = "red")
  for(i in seq_along(vals_q)){
    bp(samples[[i]][, variable], at = vals_q[i], col = col)
  }
}

# plot:
par(las = 1, mfcol = c(5, 3), mar = c(4, 4.5, 1.8, 1), font.main = 1, family = "serif")

# my method:
all_bps(samples = estim, variable = "nu", vals_q = vals_q, ylim = c(0, 60),
        ylab = expression(hat(nu)), true_value = nu)
mtext("(a) Accounting for underreporting", side = 3, at = 0.5, cex = 0.75, line = 0.9)

all_bps(samples = estim, variable = "phi", vals_q = vals_q,
        ylab = expression(hat(phi1)), true_value = phi)

all_bps(samples = estim, variable = "kappa", vals_q = vals_q,
        ylab = expression(hat(kappa)), true_value = kappa)

all_bps(samples = estim, variable = "psi", vals_q = vals_q, ylim = c(0, 0.25),
        ylab = expression(hat(psi)), true_value = psi)

all_bps(samples = estim, variable = "Reff", vals_q = vals_q, ylim = c(0, 1),
        ylab = expression(hat(R)[eff]), true_value = phi/(1 - kappa))


# naive:
all_bps(samples = estim_naive, variable = "nu", vals_q = vals_q, ylim = c(0, 60),
        ylab = expression(hat(nu)), expected_bias = expec_naive["nu", ], true_value = nu)
mtext("(b) Ignoring underreporting", side = 3, at = 0.5, cex = 0.75, line = 0.9)

all_bps(samples = estim_naive, variable = "phi", vals_q = vals_q,
        ylab = expression(hat(phi1)), expected_bias = expec_naive["phi", ], true_value = phi)

all_bps(samples = estim_naive, variable = "kappa", vals_q = vals_q,
        ylab = expression(hat(kappa)), expected_bias = expec_naive["kappa", ], true_value = kappa)

all_bps(samples = estim_naive, variable = "psi", vals_q = vals_q, ylim = c(0, 0.25),
        ylab = expression(hat(psi)), expected_bias = expec_naive["psi", ], true_value = psi)

all_bps(samples = estim_naive, variable = "Reff", vals_q = vals_q, ylim = c(0, 1),
        ylab = expression(hat(R)[eff]),
        expected_bias = expec_naive["phi", ]/(1 - expec_naive["kappa", ]), true_value = phi/(1 - kappa))

# stretched:
all_bps(samples = estim_stretch, variable = "nu", vals_q = vals_q, ylim = c(0, 60),
        ylab = expression(hat(nu)), true_value = nu)
mtext("(c) Using multiplication factors", side = 3, at = 0.5, cex = 0.75, line = 0.9)

all_bps(samples = estim_stretch, variable = "phi", vals_q = vals_q,
        ylab = expression(hat(phi1)), expected_bias = expec_stretched["phi", ], true_value = phi)

all_bps(samples = estim_stretch, variable = "kappa", vals_q = vals_q,
        ylab = expression(hat(kappa)), expected_bias = expec_stretched["kappa", ], true_value = kappa)

all_bps(samples = estim_stretch, variable = "psi", vals_q = vals_q, ylim = c(0, 0.25),
        ylab = expression(hat(psi)), expected_bias = expec_stretched["psi", ], true_value = psi)

all_bps(samples = estim_stretch, variable = "Reff", vals_q = vals_q, ylim = c(0, 1),
        ylab = expression(hat(R)[eff]),
        expected_bias = expec_stretched["phi", ]/(1 - expec_stretched["kappa", ]), true_value = phi/(1 - kappa))

