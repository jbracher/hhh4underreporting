# Comparing properties of an underreported process and the approximating
# fully observed one

# Johannes Bracher, johannes.bracher@uzh.ch

setwd("/home/johannes/Documents/underreporting/Article_hhh4/simulations/equivalance_distr")

library(hhh4underreporting)

# Comparing the stationary and conditional distributions of second-order
# equivalent processes with different reporting probabilities

# specify model parameters
nu <- 20
phi <- 0.5
kappa <- 0.2
psi <- 0.1
q <- 0.2
pars_q <- list(nu = nu, phi = phi, kappa = kappa, psi = psi, q = q)

# parameters of second-order equivalent completely observed process:
pars1 <- hhh4underreporting:::reparam_homog(nu, phi, kappa, psi, q)
nu1 <- pars1$nu
phi1 <- pars1$phi
kappa1 <- pars1$kappa
psi1 <- pars1$psi
q1 <- 1

# length of simulated time series:
lgt <- 100000

# simulate from both versions:
set.seed(8)

# with q = 1
ts1 <- simulate_hhh4u_homog(nu = nu1,
                            phi = phi1,
                            kappa = kappa1,
                            psi = psi1,
                            q = q1,
                            lgt = lgt)
X_tilde1 <- ts1$X_tilde
lambda1 <- ts1$lambda

# with q = 1
ts_q <- simulate_hhh4u_homog(nu = nu,
                             phi = phi,
                             kappa = kappa,
                             psi = psi,
                             q = q,
                             lgt = lgt)
X_tilde_q <- ts_q$X_tilde
lambda_q <- ts_q$lambda

# Rao-Blackwellize using the lambda vectors which are also returned by the simulation functions:
support <- 0:350

rbize <- function(lambda, psi, q, support = 0:100){
  dens <- matrix(ncol = length(support), nrow = length(lambda))
  for(i in 1:nrow(dens)){
    dens[i, ] <- dnbinom(support, mu = q*lambda[i], size = 1/psi)
  }
  colMeans(dens)
}

# full distributions:
rb1 <- rbize(lambda1, pars1$psi, 1, support)
rb_q <- rbize(lambda_q, psi, q, support)

# how much probability mass would we have to re-distibute?
sum(abs(rb1 - rb_q))/2 # half a percent

# summarize results:
marginal_densities <- data.frame(X = support, density1 = rb1, density_q = rb_q)


# get conditional distributions:
binned_rb <- function(lambda, psi, x, breaks_x, support){
  x_binned <- cut(x, breaks_x)
  lev <- levels(x_binned)
  cond_means <- cond_var <- cond_skew <- cond_exkurt <- rep(NA, length(lev))
  names(cond_means) <- names(cond_var) <-
    names(cond_skew) <- names(cond_exkurt) <- lev
  for(i in seq_along(lev)){
    inds <- which(x_binned == lev[i])
    if(length(inds) == 0) next()
    if(length(inds) > 1000) inds <- sample(inds, 1000)
    rb <- rbize(lambda[inds], psi, 1, support)
    cond_means[i] <- sum(support*rb)
    cond_var[i] <- sum(support^2*rb) - cond_means[i]^2
    cond_skew[i] <- (sum(support^3*rb) - 3*cond_var[i]*cond_means[i] - cond_means[i]^3)/
      (cond_var[i]^1.5)
    cond_exkurt[i] <- (sum(support^4*rb) - 4*cond_means[i]*sum(support^3*rb) +
                         6*sum(support^2*rb)*cond_means[i]^2 - 3*cond_means[i]^4)/ # stand. fourth moment
      (cond_var[i])^2 - 3
  }
  return(list(mean = cond_means, var = cond_var, skew = cond_skew, exkurt = cond_exkurt))
}

# helper function for lagging:
lag_vector <- function(vect, lag){
  c(rep(NA, lag), head(vect, length(vect) - lag))
}

# compute lambda_star for X_tilde_q:
lambda_star_from_X_tilde_q <- nu1 +
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
breaks_x <- seq(from = support[1], to = support[length(support)], by = 0.25)

# get conditional distributions for ranges of lambda_star
cond_rb_q <- binned_rb(lambda = lambda_star_from_X_tilde_q, psi = psi1,
                       x = lambda_star_from_X_tilde_q, breaks_x = breaks_x,
                       support = support)
cond_rb1 <- binned_rb(lambda = lambda1, psi = psi1,
                      x = lambda1, breaks_x = breaks_x,
                      support = support)

# summarize results:
conditional_moments <- data.frame(lambda_star = breaks_x[-1],
                                  cond_mean_q = cond_rb_q$mean,
                                  cond_var_q = cond_rb_q$var,
                                  cond_skew_q = cond_rb_q$skew,
                                  cond_excurt_q = cond_rb_q$exkurt,
                                  cond_mean1 = cond_rb1$mean,
                                  cond_var1 = cond_rb1$var,
                                  cond_skeq1 = cond_rb1$skew,
                                  cond_exkurt1 = cond_rb1$exkurt)


# figure:

layout(matrix(c(1, 1, 2,
                3, 4, 5), ncol = 3, byrow = TRUE))
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

# compare conditional means:
plot(NULL, xlim = xl, ylim = c(0, 50),
     xlab = expression(lambda^"*"),
     ylab = "cond. mean", las = 1)
lines(0:100, 0:100, col = col1)
lines(breaks_x[-1], conditional_moments$cond_mean_q, type = "l", lty = 2, col = col_q)

# compare conditional variances:
plot(NULL, xlim = xl, ylim = xl,
     xlab = expression(lambda^"*"),
     ylab = "cond. sd")
lines(0:100, sqrt(0:100 + (0:100)^2*psi1), col = col1)
lines(breaks_x[-1], sqrt(conditional_moments$cond_var_q), type = "l", lty = 2, col = col_q)

# compare conditional skewness:
p1 <- 1/psi1/(1:100 + 1/psi1)
skew1 <- (2 - p1)/sqrt(1/psi1*(1 - p1))
plot(skew1, col = col1, type = "l", xlim = xl, ylim = c(0.5, 0.8),
     xlab = expression(lambda^"*"),
     ylab = "cond. skewness")
lines(breaks_x[-1], conditional_moments$cond_skew_q, lty = "dashed", col = col_q)

# compare conditional excess kurtosis:
curt1 <- 6*psi1 + psi1*p1^2/(1 - p1)
plot(curt1, col = col1, type = "l", xlim = xl, ylim = c(0.5, 1.2),
     xlab = expression(lambda^"*"),
     ylab = "cond. excess kurtosis")
lines(breaks_x[-1], conditional_moments$cond_excurt_q, lty = "dashed", col = col_q)


write.csv(marginal_densities, file = "results_marginal_densities.csv", row.names = FALSE)
write.csv(conditional_moments, file = "results_conditional_moments.csv", row.names = FALSE)
write.csv(rbind(pars_q, pars1), file = "settings_equivalence_distr.csv", row.names = FALSE)
