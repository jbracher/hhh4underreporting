% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homogeneous.R
\name{simulate_hhh4u_homog}
\alias{simulate_hhh4u_homog}
\title{Simulate from time-homogeneous underreported model}
\usage{
simulate_hhh4u_homog(
  nu,
  phi,
  kappa,
  psi,
  q = 1,
  lgt = 100,
  start = 10,
  burn_in = 1000
)
}
\arguments{
\item{nu, phi, kappa, psi}{the model parameters (scalars)}

\item{q}{the reporting probability}

\item{lgt}{the length of the generated time series (integer)}

\item{start}{initial value of both $X$ and $lambda$ (at beginning of burn in period)}

\item{burn_in}{number of iterations to discard to reach stationary phase}
}
\value{
A named list with elements \code{"X"}, \code{"X_tilde"} and \code{"lambda"} containing
the unthinned and thinned simulated time series as well as the conditional mean process.
}
\description{
Simulate from a time-homogeneous endemic-epidemic model with underreporting. Includes a burn-in period
to reach stationary phase.
}
\examples{
sim <- simulate_hhh4u(nu = 3, phi = 0.4, kappa = 0.2, psi = 0.1, q = 0.5, lgt = 100)
}
