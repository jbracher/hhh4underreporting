# Comparing the forward pass algorithm and the suggested likelihood approximation

# Johannes Bracher, johannes.bracher@uzh.ch

library(hhh4underreporting)

setwd("/home/johannes/Documents/underreporting/Article_hhh4/simulations/comparison_fp")

#######
### 1. Likelihood evaluation for of hidden INARCH(1) model via forward algorithm

# Transition matrix of an INARCH(1) process
#
# Arguments:
# nu: immigration parameter of the INARCH(1)
# phi: autoregressive parameter of the INARCH(1)
# upper limit: an upper limit at which to truncate the state space of the INAR(1)
# Return: the transition matrix
get_transition_matrix_inarch_nb <- function(nu, phi, psi, upper_limit){
  n_states <- upper_limit + 1
  immigration_probs <- dpois(0:upper_limit, nu)

  matr <- matrix(nrow = n_states, ncol = n_states)

  for(from in 0:upper_limit){
    for(to in 0:upper_limit){
      matr[from + 1, to + 1] = dnbinom(to, mu = nu + phi*from, size = 1/psi)
    }
  }
  return(matr)
}

# Evaluate log-likelihood of hidden INARCH(1) model:
#
# Arguments:
# vect: a vector containing the observed time series
# nu, phi, q: the parameters of the hidden INARCH(1)
# log: should log-likelihood or untransformed likelihood be returned?
# return_fp: should forward probabilities be returned?
# Return: the log-likelihood or (if return_fp) a list containing the
# log-likelihood and a matrix of forward probabilities
llik_fp <- function(vect, lambda1, nu, phi, psi, q, log = TRUE, return_contributions = FALSE){
  lgt <- length(vect)
  upper_limit <- round(1.5*max(vect)/q) # max(qnbinom(0.99, mu = nu/(1 - phi), size = size_stat), round(1.5*max(vect)/q))
  # print(c(qnbinom(0.999, mu = nu/(1 - phi), size = size_stat), round(1.5*max(vect)/q)))

  support <- 0:upper_limit
  transition_matrix <- get_transition_matrix_inarch_nb(nu = nu, phi = phi, psi = psi, upper_limit = upper_limit)

  fp <- matrix(nrow = lgt, ncol = upper_limit + 1)
  normalizing <- numeric(lgt)

  # idea: normalize after ever step and keep normalizing constants
  sop <- hhh4underreporting:::compute_sop_homog(nu = nu, phi = phi, kappa = 0, psi = psi, q = q)
  size_stat <- sop$X$mu^2/(sop$X$sigma2 - sop$X$mu)
  probs_temp <- dnbinom(support, mu = lambda1, size = 1/psi)*dbinom(vect[1], support, q)
  normalizing[1] <- sum(probs_temp)
  fp[1, ] <- probs_temp/normalizing[1]

  for(i in 2:lgt){
    probs_temp <- (fp[i - 1, ] %*% transition_matrix)*dbinom(vect[i], support, q)
    normalizing[i] <- sum(probs_temp)
    fp[i, ] <- probs_temp/normalizing[i]
  }

  return_value <- ifelse(log, sum(log(normalizing)), prod(normalizing))
  if(return_contributions){
    return(list(value = return_value, contributions = normalizing))
  }else{
    return(return_value)
  }
}

# wrapper to simulate many time series from a hidden INAR(1) model:
sim_many_hidden_inarch_nb <- function(lambda1, nu, phi, psi, q, lgt, n_sim){
  # bring into vector form:
  nu <- rep(nu, lgt)
  phi <- rep(phi, lgt)
  kappa <- rep(0, lgt)
  psi <- rep(psi, lgt)
  #initialize matrix for scoring:
  X_samples <- X_tilde_samples <- matrix(ncol = lgt, nrow = n_sim)
  for(i in 1:n_sim){
    sim_temp <- simulate_hhh4u(lambda1 = lambda1, nu = nu, phi = phi, kappa = kappa,
                               psi = psi, q = q, return_sts = FALSE)
    X_samples[i, ] <- sim_temp$X
    X_tilde_samples[i, ] <- sim_temp$X_tilde
  }
  return(list(X = X_samples, X_tilde = X_tilde_samples))
}

# evaluate empirical proportion of a given sequence (used in empirical likelihood
# evaluation of short sequences):
emp_prob <- function(Y_samples, vect){
  mean(colSums(t(Y_samples) == vect) == ncol(Y_samples))
}

# check that likelihood evaluation with forward algorithm works:

# set true parameters of hidden INARCH(1):
lambda1 <- 3
nu <- 3
phi <- 0.65
q <- 0.1
lgt <- 3
psi <- 0.1

# vector for which to evaluate likelihood:
vect <- c(0, 2, 1)

# get many simulations from the model (short chains, otherwise too much computation effort):
sim_hi <- sim_many_hidden_inarch_nb(lambda1, nu, phi, psi, q, lgt = 3, n_sim = 100000)
# evaluate likelihood empirically:
emp_prob(sim_hi$X_tilde, vect)
# evaluate analytically:
llik_fp(vect = vect, lambda1 = lambda1, nu = nu,
        phi = phi, psi = psi, q = q, log = FALSE)
# this works


#######
### 2. Compare results from forward algorithm and the approximation method

n_sim <- 1000
lgt <- 100

# Initialize vectors and lists to store results:
lambda1_true <- nu_true <- phi_true <- q_true <- psi_true <-
  nu_eval <- phi_eval <- q_eval <- psi_eval <- llik_fa <- llik_approx <-
  time_approx <- time_fa <- NA*numeric(n_sim)
sim <- list()

# set seed
# set.seed(1)

# run simulation:
for(i in 1:n_sim){
  set.seed(i)
  # sample true parameter values, discarding caombinations implying non-stationarity:
  while (is.na(phi_true[i]) | (phi_true[i] + phi_true[i]^2*psi_true[i]) > 1) {
    nu_true[i] <- runif(1, 1, 30)
    phi_true[i] <- runif(1, 0.01, 0.99)
    q_true[i] <- runif(1, 0.01, 1)
    psi_true[i] <- runif(1, 0.01, 0.2) # runif(1, 0.1, min(0.2, 0.5*(1 - phi_true[i]^2)/phi_true[i]^2))
    lambda1_true[i] <- runif(1, 0.5*nu_true[i]/(1 - phi_true[i]), 2*nu_true[i]/(1 - phi_true[i]))
  }

  # sample time series:
  sim[[i]] <- simulate_hhh4u(lambda1 = lambda1_true[i],
                             nu = rep(nu_true[i], lgt),
                             phi = rep(phi_true[i], lgt),
                             kappa = rep(0, lgt),
                             psi = rep(psi_true[i], lgt),
                             q = q_true[i], return_sts = FALSE)$X_tilde

  # where to evaluate the log likelihood?
  lambda1_eval <- lambda1_true[i]
  nu_eval[i] <- nu_true[i]
  phi_eval[i] <- phi_true[i]
  q_eval[i] <- q_true[i]
  psi_eval[i] <- psi_true[i]

  # evaluate log-likelihood using approximation:
  t1 <- Sys.time() # to measure computing time
  llik_approx[i] <- -1*hhh4underreporting:::nllik_tv_cpp(
    observed = sim[[i]], lambda1 = lambda1_eval,
    nu = rep(nu_eval[i], lgt),
    phi = rep(phi_eval[i], lgt),
    kappa = rep(0, lgt),
    psi = rep(psi_eval[i], lgt),
    q = rep(q_eval[i], lgt)
  )
  t2 <- Sys.time()
  time_approx[i] <- t2 - t1

  # evaluate using forward algorithm
  t1 <- Sys.time()
  llik_fa[i] <- llik_fp(sim[[i]], lambda1_eval,
                        nu_eval[i], phi_eval[i],
                        psi_eval[i], q_eval[i], log = TRUE)
  t2 <- Sys.time()
  time_fa[i] <- t2 - t1

  # plot every twenty iterations:
  if(i %% 20 == 0){
    print(i)
    boxplot(llik_approx - llik_fa)
  }
}

# assemble results
dat <- cbind(nu = nu_true, phi = phi_true, psi = psi_true, q = q_true,
             llik_approx = llik_approx, llik_fa = llik_fa,
             time_approx = time_approx, time_fa = time_fa)

# store:
write.csv(dat, file = paste0("results_comparison_fp_general_", lgt,  ".csv"))
# dat <- read.csv("results_comparison_fp_general.csv")

save(sim, lambda1_true, nu_true, phi_true, psi_true,
     dat, llik_approx, llik_fa, file = paste0("results_comparison_fp_general_", lgt, ".rda"))

dat <- as.data.frame(dat)
dat$diff <- dat$llik_approx - dat$llik_fa
dat$mu_tilde <- dat$q*dat$nu/(1 - dat$phi)


# compute percentage of cases where diff is below 0.01 and below 1:
percentage_below0.1 <- round(100*mean(abs(dat$diff) < 0.1), 2)
percentage_below1 <- round(100*mean(abs(dat$diff) < 1), 2)

# compute median ratio of computing times:
ratio_computing_times <- round(median(dat$time_fa / dat$time_approx))

# plot as in manuscript:
# structure plot area:
par(mfrow = c(1, 3), las = 1, family = "serif",
    mar = c(4, 4, 0.1, 0.1), cex = 1.1)
yl <- c(-1, 1)*max(abs(dat$diff))

# plot difference in log-likelihood against mu_tilde:
plot(log(dat$mu_tilde), dat$diff, ylab = "Difference in log-likelihood",
     xlab = expression(tilde(mu)), ylim = yl, pch = 15, cex = 0.3, axes = FALSE)
vals_x_axis <- c(0.1, 1, 10, 100)
axis(2); box()
axis(1, at = log(vals_x_axis), labels = vals_x_axis)

# plot difference in log-likelihood against phi*sqrt(1 + psi):
plot(sqrt(dat$phi^2*(1 + dat$psi)), dat$diff, ylab = "", axes = FALSE,
     xlim = 0:1, xlab = expression(phi1*sqrt(1 + psi)), ylim = yl, pch = 16, cex = 0.3)
axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "", "0.5", "", "1"))
axis(2); box()

# plot difference in log-likelihood against pi:
plot(dat$q, dat$diff, ylab = "", axes = FALSE,
     xlab = expression(pi), ylim = yl, pch = 16, cex = 0.3)
axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "", "0.5", "", "1"))
axis(2); box()

