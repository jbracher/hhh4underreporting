# Case study on rotavirus gastroenteritis in Berlin, Germany
# Johannes Bracher, johannes.bracher@uzh.ch

library(hhh4underreporting)

setwd("/home/johannes/Documents/underreporting/Article_hhh4/case_study")

# get data
data("rota_germany")
region <- "Berlin"
dat <- rota_germany[rota_germany$year %in% 2001:2008, region] # restrict to years 2001-2008

# for numerical stability we increase the first observation from 0 to 1
# otherwise problems with log(lambda1) goes to minus infinity and numerical problems occur
# when differentiating the log likelihood function (note that point estimates are unaffected
# by this minor change).
dat[1] <- 1
sts_dat <- sts(observed = dat, start = c(2001, 1))

plot(sts_dat)

# define range for reporting probability:
n_steps_q <- 50
values_q <- 1:n_steps_q/n_steps_q

# get yearly sums:
yearly_sums <- colSums(matrix(dat, nrow = 52))

#############################################################
# analyse with seasonality only in endemic component:
fits_seas_end <- list()
pars_seas_end <- matrix(ncol = 7, nrow = length(values_q))
colnames(pars_seas_end) <- c("log_lambda1", "alpha_nu", "gamma_nu", "delta_nu", "alpha_phi", "alpha_kappa", "log_psi")
for(i in 1:n_steps_q){
  # store fit
  fits_seas_end[[i]] <- hhh4u(sts_dat,
                             control = list(end = list(f = addSeason2formula(~1)),
                                            family = "NegBin1",
                                            q = values_q[i],
                                            return_se = FALSE))
  # store parameters
  pars_seas_end[i, ] <- fits_seas_end[[i]]$coefficients
  print(i)
}
# save(fits_seas_end, pars_seas_end, file = paste0("results/fits_rota_seas_end_", region, "_0105.rda"))
# load("/home/johannes/Documents/underreporting/real_data/rotavirus_hhh4latent/likelihood/fits_seas_end_rota.rda")

# look at parameter values:
par(mfrow = c(2, 2))

plot(values_q, pars_seas_end[, "alpha_nu"], type = "b")
plot(values_q, pars_seas_end[, "gamma_nu"], ylim = c(0, 1), type = "b")
lines(values_q, pars_seas_end[, "delta_nu"], type = "b")

plot(values_q, exp(pars_seas_end[, "alpha_phi"]), ylim = c(0, 1), type = "b")
lines(values_q, exp(pars_seas_end[, "alpha_kappa"]), type = "b")
abline(v = 0.2)
abline(h = 0.79)

plot(values_q, exp(pars_seas_end[, "log_psi"]), type = "b", ylim = c(0, 0.2))

# check convergence:
for(i in 1:n_steps_q) print(fits_seas_end[[i]]$convergence) # all successfully converged

#############################################################
# analyse with seasonality in both components:
fits_seas_both <- list()
pars_seas_both <- matrix(ncol = 9, nrow = length(values_q))
colnames(pars_seas_both) <- c("log_lambda1",
                              "alpha_nu", "gamma_nu", "delta_nu",
                              "alpha_phi", "gamma_phi", "delta_phi",
                              "alpha_kappa", "log_psi")


for(i in n_steps_q:1){
  print(paste("Starting ", i))

  # we here re-run the fitting function three times, always using the result from the
  # previous iteration as the starting value. This improves convergence.
  fits_seas_both[[i]] <- hhh4u(sts_dat,
                               control = list(end = list(f = addSeason2formula(~1)),
                                              ar = list(f = addSeason2formula(~1)),
                                              family = "NegBin1",
                                              q = values_q[i],
                                              return_se = FALSE))
  print(paste("First: loglik =", fits_seas_both[[i]]$loglikelihood))

  fits_seas_both[[i]] <- hhh4u(sts_dat,
                               control = list(end = list(f = addSeason2formula(~1)),
                                              ar = list(f = addSeason2formula(~1)),
                                              family = "NegBin1",
                                              q = values_q[i],
                                              return_se = FALSE,
                                              start = fits_seas_both[[i]]$coefficients))
  print(paste("Second: loglik =", fits_seas_both[[i]]$loglikelihood))

  fits_seas_both[[i]] <- hhh4u(sts_dat,
                               control = list(end = list(f = addSeason2formula(~1)),
                                              ar = list(f = addSeason2formula(~1)),
                                              family = "NegBin1",
                                              q = values_q[i],
                                              return_se = FALSE,
                                              start = fits_seas_both[[i]]$coefficients))
  print(paste("Third: loglik =", fits_seas_both[[i]]$loglikelihood))

  # store parameters:
  pars_seas_both[i, ] <- fits_seas_both[[i]]$coefficients
}

# save(fits_seas_both, pars_seas_both, file = paste0("results/fits_rota_seas_both_", region, "0108_general.rda"))
# load(paste0("results/fits_rota_seas_both_", region, "0108_general.rda"))
write.csv(pars_seas_both, paste0("results/fits_rota_seas_both_", region, "0108_general.csv"), row.names = FALSE)

# compute Reff in the form of intercept and sine/cosine waves (required for plot function)
Reff <- pars_seas_both[, 5:7]
Reff[, 1] <- Reff[, 1] - log(1 - exp(pars_seas_both[, "alpha_kappa"]))
colnames(Reff) <- c("alpha_Reff", "gamma_Reff", "delta_Reff")

# check convergence:
for(i in 1:n_steps_q) print(fits_seas_both[[i]]$convergence) # all successfully converged

#############################################################
# plot as in manuscript:

# helper function to plot parameters over the course of a season:
lines_seas <- function(vect, exp = FALSE, ...){
  x <- 1:52
  if(exp){
    y <- exp(vect[1] + vect[2]*sin(2*pi*x/52) + vect[3]*cos(2*pi*x/52))
  }else{
    y <- vect[1] + vect[2]*sin(2*pi*x/52) + vect[3]*cos(2*pi*x/52)
  }
  lines(x, y, ...)
}

# define colours:
library(viridis, quietly = TRUE)
cols <- rev(viridis(n_steps_q))
cols[2] <- "red"

# structure plot area:
layout(matrix(c(1, 2,
                1, 2,
                3, 5,
                4, 5), ncol = 2, byrow = TRUE))
order_plot <- c(1:n_steps_q, 2)
par(mar = c(5.1, 4.5, 0.2, 1.0), las = 1, cex = 0.9, font.main = 1, family = "serif")

# plot nu:
plot(NULL, ylim = c(0, 8), xlim = c(1, 52), xlab = "calendar week", ylab = expression(hat(nu)[t]), axes = FALSE)
axis(1, at = c(1, 13, 26, 39, 52))
box(); axis(2, at = log(c(1, 2, 5, 10, 20, 50, 100, 250, 500, 1000)), labels = c(1, 2, 5, 10, 20, 50, 100, 250, 500, 1000))
for(i in order_plot) lines_seas(pars_seas_both[i, 2:4], col = cols[i], exp = FALSE)

# plot phi:
plot(NULL, ylim = c(0, 1.15), xlim = c(1, 52), xlab = "calendar week",
     ylab = expression(hat(phi1)[t]), axes = FALSE)
axis(1, at = c(1, 13, 26, 39, 52))
axis(2); box()
for(i in order_plot) lines_seas(pars_seas_both[i, 5:7], col = cols[i], exp = TRUE)
vals_x <- seq(from = 3, to = 47, length.out = n_steps_q)
text(26, 1.1, expression(pi))
points(vals_x, rep(1.05, n_steps_q), col = cols, pch = 15)
text(vals_x[c(1, 12, 25, 37, 50)], 1, labels = c(".02", ".25", ".5", ".75", "1"), cex = 0.8)

# plot kappa:
plot(NULL, xlim = 0:1, ylim = c(0, 0.35), axes = FALSE,
     type = "l", xlab = expression(pi), ylab = expression(hat(kappa)))
points(values_q, exp(pars_seas_both[, "alpha_kappa"]), pch = 16, cex = 0.7, col = cols)
axis(1)
axis(2, at = 0:3/10)
box()

# plot psi:
plot(NULL, xlim = 0:1,  ylim = c(0.08, 0.11), axes = FALSE,
     type = "l", xlab = expression(pi), ylab = expression(hat(psi)))
axis(1)
axis(2, at = c(0.08, 0.09, 0.1, 0.11))
box()
points(values_q, exp(pars_seas_both[, "log_psi"]), pch = 16, cex = 0.7, col = cols)

# plot Reff
plot(NULL, ylim = c(0, 1.1), xlim = c(1, 52), xlab = "calendar week",
     ylab = expression(hat(R)["eff, t"]), axes = FALSE)
axis(1, at = c(1, 13, 26, 39, 52))
axis(2); box()
for(i in order_plot) lines_seas(Reff[i, ], col = cols[i], exp = TRUE)

#############################################################
# get numbers for table:

# model with reporting probability 0.043:
# fit model:
ctrl_0043 <- list(end = list(f = addSeason2formula(~1)),
                  ar = list(f = addSeason2formula(~1)),
                  family = "NegBin1",
                  return_se = TRUE,
                  q = 0.043)
fit_0043.0 <- hhh4u(sts_dat, ctrl_0043) # re-fit once to improve convergence
ctrl_0043$start <- fit_0043.0$coefficients
fit_0043 <- hhh4u(sts_dat, ctrl_0043)
# for comparison:
fit_0043$loglikelihood - fits_seas_both[[2]]$loglikelihood
fit_0043$coefficients / fits_seas_both[[2]]$coefficients
# get entries for table:
entries_tab_0043 <- paste0(format(fit_0043$coefficients, digits = 1),
                           " (", format(fit_0043$se, digits = 1, scientific = FALSE), ")")
names(entries_tab_0043) <- c("log_lambda1", "alpha_nu", "gamma_nu", "delta_nu",
                             "alpha_phi", "gamma_phi", "delta_phi",
                             "alpha_kappa", "log_psi")

# model with reporting probability 0.043, then 0.063:
# define time-varying reporting probability:
q_timevar <- c(rep(0.043, 4*52), seq(from = 0.043, to = 0.063, length.out = 52),
               rep(0.063, 3*52))
# fit model:
ctrl_0043.0063 <- list(end = list(f = addSeason2formula(~1)),
                       ar = list(f = addSeason2formula(~1)),
                       family = "NegBin1",
                       return_se = TRUE,
                       q = q_timevar)
fit_0043.0063.0 <- hhh4u(sts_dat, control = ctrl_0043.0063)
ctrl_0043.0063$start <- fit_0043.0063.0$coefficients
fit_0043.0063 <- hhh4u(sts_dat, ctrl_0043.0063)
# get entries for table:
entries_tab_0043.0063 <- paste0(format(fit_0043.0063$coefficients, digits = 1),
                                " (", format(fit_0043.0063$se, digits = 1, scientific = FALSE), ")")
names(entries_tab_0043.0063) <- names(entries_tab_0043)


# model with reporting probability 1:
# fit model:
ctrl_1 <- list(end = list(f = addSeason2formula(~1)),
               ar = list(f = addSeason2formula(~1)),
               family = "NegBin1",
               q = 1)
fit_1.0 <- hhh4u(sts_dat, control = ctrl_1)
ctrl_1$start <- fit_1.0$coefficients
fit_1 <- hhh4u(sts_dat, ctrl_1)
# get entries for table:
entries_tab_1 <- paste0(format(fit_1$coefficients, digits = 1),
                        " (", format(fit_1$se, digits = 1, scientific = FALSE), ")")
names(entries_tab_1) <- names(entries_tab_0043)

tab <- rbind(q0043 = entries_tab_0043,
             q0043.0063 = entries_tab_0043.0063,
             q1 = entries_tab_1)


##########################################
# plot of confidence bands:

library(mvtnorm)

# helper functions for sampling-based computation of confidence bands:
# for Reff
sample_quantiles_Reff <- function(fit, n_sim = 10000, probs = c(0.05, 0.95)){
  samples_par <- rmvnorm(n_sim, fit$coefficients, fit$cov)
  samples_Reff <- matrix(ncol = 52, nrow = n_sim)

  for(i in 1:n_sim){
    samples_Reff[i, ] <- exp(samples_par[i, "ar.(Intercept)"] +
                               samples_par[i, "ar.sin(2 * pi * t/52)"]*sin(2*pi*1:52/52) +
                               samples_par[i, "ar.cos(2 * pi * t/52)"]*cos(2*pi*1:52/52))/
      (1 - exp(samples_par[i, "log_kappa"]))
  }
  quantiles_Reff <- apply(samples_Reff, 2, quantile, probs = c(probs[1], 0.5, probs[2]))
  quantiles_Reff[2, ] <- exp(fit$coefficients["ar.(Intercept)"] +
                               fit$coefficients["ar.sin(2 * pi * t/52)"]*sin(2*pi*1:52/52) +
                               fit$coefficients["ar.cos(2 * pi * t/52)"]*cos(2*pi*1:52/52))/
    (1 - exp(fit$coefficients["log_kappa"]))
  return(quantiles_Reff)
}

# for nu:
sample_quantiles_nu <- function(fit, n_sim = 10000, probs = c(0.05, 0.95)){
  samples_par <- rmvnorm(n_sim, fit$coefficients, fit$cov)
  samples_nu <- matrix(ncol = 52, nrow = n_sim)

  for(i in 1:n_sim){
    samples_nu[i, ] <- exp(samples_par[i, "end.(Intercept)"] +
                             samples_par[i, "end.sin(2 * pi * t/52)"]*sin(2*pi*1:52/52) +
                             samples_par[i, "end.cos(2 * pi * t/52)"]*cos(2*pi*1:52/52))
  }

  quantiles_nu <- apply(samples_nu, 2, quantile, probs = c(probs[1], 0.5, probs[2]))
  # replace middle column by point estimates:
  quantiles_nu[2, ] <- exp(fit$coefficients["end.(Intercept)"] +
                             fit$coefficients["end.sin(2 * pi * t/52)"]*sin(2*pi*1:52/52) +
                             fit$coefficients["end.cos(2 * pi * t/52)"]*cos(2*pi*1:52/52))
  return(quantiles_nu)
}

# obtain confidence bands:
set.seed(123)

quantiles_Reff_0043 <- sample_quantiles_Reff(fit_0043)
quantiles_Reff_0043.0063 <- sample_quantiles_Reff(fit_0043.0063)
quantiles_Reff_1 <- sample_quantiles_Reff(fit_1)

quantiles_nu_0043 <- sample_quantiles_nu(fit_0043)
quantiles_nu_0043.0063 <- sample_quantiles_nu(fit_0043.0063)
quantiles_nu_1 <- sample_quantiles_nu(fit_1)

# helper functions for plotting:
plot_quantiles <- function(quantiles, add = FALSE, col = "black", lwd_point_est = 2,...){
  if(!add){
    plot(quantiles[2, ], ylim =c(0, 1.6), type = "l", col = col, lwd = lwd_point_est, ...)
  }else{
    lines(quantiles[2, ], col = col, lwd = lwd_point_est)
  }

  lines(quantiles[1, ], lty = 3, col = col)
  lines(quantiles[3, ], lty = 3, col = col)
}


plot_quantiles_nu <- function(quantiles, add = FALSE, col = "black", lwd_point_est = 2,...){
  if(!add){
    plot(log(quantiles[2, ]), type = "l", col = col, lwd = lwd_point_est, ...)
  }else{
    lines(log(quantiles[2, ]), col = col, lwd = lwd_point_est)
  }

  lines(log(quantiles[1, ]), lty = 3, col = col)
  lines(log(quantiles[3, ]), lty = 3, col = col)
}

#plot:
# nu:
par(mfrow = c(1, 2), las = 1, font.main = 1, family = "serif", mar = c(3, 4.5, 1, 1), cex = 0.8)

plot_quantiles_nu(quantiles_nu_0043, xlab = "calendar week",
                  ylab = expression(hat(nu)[t]), col = "red",
                  axes = FALSE, ylim = c(0, 7))
axis(1); box(); axis(2, at = log(c(1, 2, 5, 10, 20, 50, 100, 250, 500, 1000)), labels = c(1, 2, 5, 10, 20, 50, 100, 250, 500, 1000))
plot_quantiles_nu(quantiles_nu_1, add = TRUE, col = cols[length(cols)])
plot_quantiles_nu(quantiles_nu_0043.0063, add = TRUE, col = "darkgreen")

# Reff:
plot_quantiles(quantiles_Reff_0043, xlab = "calendar week",
               ylab = expression(hat(R)[eff, t]), col = "red",
               axes = FALSE)
axis(1, at = c(1, 13, 26, 39, 52))
axis(2); box()
plot_quantiles(quantiles_Reff_1, add = TRUE, col = cols[length(cols)])
plot_quantiles(quantiles_Reff_0043.0063, add = TRUE, col = "darkgreen")

legend("top", legend = c(expression(pi==1),
                         expression(pi == 0.043),
                         expression(paste(pi==0.043))),
       col = c(cols[length(cols)], "red", "darkgreen"), lty = 1, lwd = 2, ncol = 3, bty = "n")
text(48, 1.43, "then 0.063")
