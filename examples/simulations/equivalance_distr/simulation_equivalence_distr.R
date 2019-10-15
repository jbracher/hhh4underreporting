# Comparing properties of an underreported process and the approximating
# fully observed one

# Johannes Bracher, johannes.bracher@uzh.ch

setwd("/home/johannes/Documents/underreporting/Article_hhh4/simulations/equivalance_distr")

library(hhh4underreporting)

# Comparing the stationary and conditional distributions of second-order
# equivalent processes with different reporting probabilities

# specify model parameters
nu <- 30
phi <- 0.5
kappa <- 0.2
psi <- 0.1
q <- 0.1
pars_q <- list(nu = nu, phi = phi, kappa = kappa, psi = psi, q = q)

# parameters of second-order equivalent completely observed process:
pars1 <- hhh4underreporting:::reparam_homog(nu, phi, kappa, psi, q)
nu1 <- pars1$nu
phi1 <- pars1$phi
kappa1 <- pars1$kappa
psi1 <- pars1$psi


# length of simulated time series:
lgt <- 3000000

# simulate from both settings (with reporting probability q and 1):
set.seed(8)

# with reporting probability 1
ts1 <- simulate_hhh4u_homog(nu = nu1,
                            phi = phi1,
                            kappa = kappa1,
                            psi = psi1,
                            q = 1,
                            lgt = lgt)
X_tilde1 <- ts1$X_tilde
lambda1 <- ts1$lambda

# with reporting probability q
ts_q <- simulate_hhh4u_homog(nu = nu,
                             phi = phi,
                             kappa = kappa,
                             psi = psi,
                             q = q,
                             lgt = lgt)
X_tilde_q <- ts_q$X_tilde
lambda_q <- ts_q$lambda

# Rao-Blackwellize using the lambda vectors which are also returned by the simulation functions:
support <- 0:100

rbize <- function(lambda, psi, q, support = 0:100){
  dens <- matrix(ncol = length(support), nrow = length(lambda))
  for(i in 1:nrow(dens)){
    dens[i, ] <- dnbinom(support, mu = q*lambda[i], size = 1/psi)
  }
  colMeans(dens)
}

# full distributions:
rb1 <- rbize(lambda1[1:100000], pars1$psi, 1, support)
rb_q <- rbize(lambda_q[1:100000], psi, q, support)

mean1 <- sum(rb1*support)
var1 <- sum(rb1*support^2) - mean1^2
skew1 <- (sum(support^3*rb1) - 3*var1*mean1 - mean1^3)/
  (var1^1.5)
excurt1 <- (sum(support^4*rb1) - 4*mean1*sum(support^3*rb1) +
              6*sum(support^2*rb1)*mean1^2 - 3*mean1^4)/ # stand. fourth moment
  (var1)^2 - 3


mean_q <- sum(rb_q*support)
var_q <- sum(rb_q*support^2) - mean_q^2
skew_q <- (sum(support^3*rb_q) - 3*var_q*mean_q - mean_q^3)/
  (var_q^1.5)
excurt_q <- (sum(support^4*rb_q) - 4*mean_q*sum(support^3*rb_q) +
               6*sum(support^2*rb_q)*mean_q^2 - 3*mean_q^4)/ # stand. fourth moment
  (var_q)^2 - 3

# how much probability mass would we have to re-distibute?
sum(abs(rb1 - rb_q))/2 # half a percent
# Kullback-Leibler distance
sum(rb_q*log(rb_q/rb1))

# summarize results:
marginal_densities <- data.frame(X = support, density1 = rb1, density_q = rb_q)


# get conditional distributions:
binned_rb <- function(lambda, psi, lambda_tilde, psi1, breaks_lambda_tilde, support){
  lambda_tilde_binned <- cut(lambda_tilde, breaks_lambda_tilde)
  lev <- levels(lambda_tilde_binned)
  cond_means <- cond_var <- cond_skew <- cond_exkurt <-
    cond_kl <- cond_mae_llik <- rep(NA, length(lev))
  names(cond_means) <- names(cond_var) <-
    names(cond_skew) <- names(cond_exkurt) <- names(cond_kl) <- lev
  for(i in seq_along(lev)){
    inds <- which(lambda_tilde_binned == lev[i])
    if(length(inds) <= 100) next()
    if(length(inds) > 5000) inds <- sample(inds, 5000)
    rb <- rbize(lambda[inds], psi, 1, support)
    cond_means[i] <- sum(support*rb)
    cond_var[i] <- sum(support^2*rb) - cond_means[i]^2
    cond_skew[i] <- (sum(support^3*rb) - 3*cond_var[i]*cond_means[i] - cond_means[i]^3)/
      (cond_var[i]^1.5)
    cond_exkurt[i] <- (sum(support^4*rb) - 4*cond_means[i]*sum(support^3*rb) +
                         6*sum(support^2*rb)*cond_means[i]^2 - 3*cond_means[i]^4)/ # stand. fourth moment
      (cond_var[i])^2 - 3
    # get KL distance to approximating distribution:
    approx_distr <- rbize(lambda_tilde[inds], psi1, 1, support)
    cond_kl[i] <- sum(rb*log(rb/approx_distr))
    cond_mae_llik[i] <- sum(rb*abs(log(rb/approx_distr)))
  }
  return(list(mean = cond_means, var = cond_var,
              skew = cond_skew, exkurt = cond_exkurt,
              kl = cond_kl, mae_llik = cond_mae_llik))
}


# helper function for lagging:
lag_vector <- function(vect, lag){
  c(rep(NA, lag), head(vect, length(vect) - lag))
}

# compute lambda_tilde for X_tilde_q:
lambda_tilde <- nu1/(1 - kappa1) +
  phi1*lag_vector(X_tilde_q, 1) +
  phi1*kappa1*lag_vector(X_tilde_q, 2) +
  phi1*kappa1^2*lag_vector(X_tilde_q, 3) +
  phi1*kappa1^3*lag_vector(X_tilde_q, 4) +
  phi1*kappa1^4*lag_vector(X_tilde_q, 5) +
  phi1*kappa1^5*lag_vector(X_tilde_q, 6) +
  phi1*kappa1^6*lag_vector(X_tilde_q, 7) +
  phi1*kappa1^7*lag_vector(X_tilde_q, 8) +
  phi1*kappa1^8*lag_vector(X_tilde_q, 9) +
  phi1*kappa1^9*lag_vector(X_tilde_q, 10)

# define support and bins:
breaks_lambda_tilde <- seq(from = support[1], to = support[length(support)], by = 0.25)
# binned:
lambda_tilde_binned <- cut(lambda_tilde, breaks_lambda_tilde)



# analyse conditional distributions (via moments and KL distanes) for ranges of lambda_tilde
cond_rb_q <- binned_rb(lambda = q*lambda_q, psi = psi,
                       lambda_tilde = lambda_tilde, psi1 = psi1,
                       breaks_lambda_tilde = breaks_lambda_tilde,
                       support = support)
# for q = 1 we can do all this analytically, but just to check things work:
cond_rb1 <- binned_rb(lambda = lambda1, psi = psi1,
                      lambda_tilde = lambda1, psi1 = psi1,
                      breaks_lambda_tilde = breaks_lambda_tilde,
                      support = support)

# summarize results:
conditional_properties <- data.frame(lambda_tilde = breaks_lambda_tilde[-1],
                                  probability =
                                    as.vector(table(lambda_tilde_binned)/length(lambda_tilde_binned)),
                                  cond_mean_q = cond_rb_q$mean,
                                  cond_var_q = cond_rb_q$var,
                                  cond_skew_q = cond_rb_q$skew,
                                  cond_excurt_q = cond_rb_q$exkurt,
                                  cond_kl_q = cond_rb_q$kl,
                                  cond_mae_q = cond_rb_q$mae_llik,
                                  cond_mean1 = cond_rb1$mean,
                                  cond_var1 = cond_rb1$var,
                                  cond_skew1 = cond_rb1$skew,
                                  cond_excurt1 = cond_rb1$exkurt,
                                  cond_kl_1 = cond_rb1$kl,
                                  cond_mae_1 = cond_rb1$mae_llik)

# figure:

layout(matrix(c(1, 2,
                1, 3), ncol = 2, byrow = TRUE))
par(mar = c(5.1, 4.5, 0.2, 1.0), las = 1)

col_q <- "black"
col1 <- "red"
xl <- c(0, 40)

# marginal:
plot(NULL, ylim = c(-0.02, 0.08), xlim = xl,
     xlab = expression(tilde(X)), ylab = expression(f(tilde(X))))
# abline(h = -2:10/100, col = "lightgrey")
# abline(v = quantile(X_tilde1, probs = c(0.01, 0.05, 0.95, 0.99)), col = "lightgrey")
lines(support, marginal_densities$density_q, type = "h", col = col_q)
lines(support + 0.3, marginal_densities$density1, type = "h", col = col1)
# lines(support,
#       (marginal_densities$density1 - marginal_densities$density_q)/
#         marginal_densities$density_q, lty = 2, col = "blue")
legend("topright", col = c("black", "red", "blue"), lty = c(1, 1, 2),
       legend = c("model with q = 0.2", "model with q = 1", "rel. difference"), bty = "b")


# compare conditional Kullback-Leibler distances:
plot(conditional_properties$lambda_tilde,
     conditional_properties$probability, xlim = c(0, 50),
     main = "", xlab = expression(tilde(lambda)[t]), type = "l")
plot(conditional_properties$lambda_tilde, conditional_properties$cond_kl_q,
     xlim = c(0, 50), type = "l",
     xlab = expression(tilde(lambda)[t]), ylab = "KL distance", ylim = c(0, 0.04))


table(lambda_tilde < 7)/length(lambda_tilde)
table(lambda_tilde > 17)/length(lambda_tilde)

# plot of moments:
par(mfrow = c(1, 3))

# compare conditional means:
plot(NULL, xlim = xl, ylim = c(0, 50),
     xlab = expression(lambda^"*"),
     ylab = "cond. mean", las = 1)
lines(0:100, 0:100, col = col1)
lines(breaks_lambda_tilde[-1], conditional_properties$cond_mean_q, type = "l", lty = 2, col = col_q)
lines(breaks_lambda_tilde[-1], conditional_properties$cond_mean1, type = "l", lty = 2, col = "blue")


# compare conditional variances:
plot(NULL, xlim = xl, ylim = xl,
     xlab = expression(lambda^"*"),
     ylab = "cond. sd")
lines(0:100, sqrt(0:100 + (0:100)^2*psi1), col = col1)
lines(breaks_lambda_tilde[-1], sqrt(conditional_properties$cond_var_q), type = "l", lty = 2, col = col_q)
lines(breaks_lambda_tilde[-1], sqrt(conditional_properties$cond_var1), type = "l", lty = 2, col = "blue")

# compare conditional skewness:
p1 <- 1/psi1/(1:100 + 1/psi1)
skew1 <- (2 - p1)/sqrt(1/psi1*(1 - p1))
plot(skew1, col = col1, type = "l", xlim = xl, ylim = c(0.5, 0.8),
     xlab = expression(lambda^"*"),
     ylab = "cond. skewness")
lines(breaks_lambda_tilde[-1], conditional_properties$cond_skew_q, lty = "dashed", col = col_q)
lines(breaks_lambda_tilde[-1], conditional_properties$cond_skew1, lty = "dashed", col = "blue")


# compare conditional excess kurtosis:
curt1 <- 6*psi1 + psi1*p1^2/(1 - p1)
plot(curt1, col = col1, type = "l", xlim = xl, ylim = c(0.5, 1.2),
     xlab = expression(lambda^"*"),
     ylab = "cond. excess kurtosis")
lines(breaks_lambda_tilde[-1], conditional_properties$cond_excurt_q, lty = "dashed", col = col_q)
lines(breaks_lambda_tilde[-1], conditional_properties$cond_excurt1, lty = "dashed", col = "blue")



write.csv(marginal_densities, file = "results_marginal_densities.csv", row.names = FALSE)
write.csv(conditional_properties, file = "results_conditional_properties.csv", row.names = FALSE)
write.csv(rbind(pars_q, pars1), file = "settings_equivalence_distr.csv", row.names = FALSE)
