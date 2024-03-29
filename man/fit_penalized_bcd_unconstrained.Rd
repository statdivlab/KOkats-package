% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_penalized_bcd_unconstrained.R
\name{fit_penalized_bcd_unconstrained}
\alias{fit_penalized_bcd_unconstrained}
\title{Fit Algorithm 2 - DEPRICATED because fit_peanlized_bcd_unconstrained_fast is better
Estimate B and z parameters through block coordinate descent to maximize the Firth penalized log likelihood function, details given in Algorithm 2.}
\usage{
fit_penalized_bcd_unconstrained(
  formula_rhs = NULL,
  Y,
  X = NULL,
  covariate_data = NULL,
  z = NULL,
  tolerance = 1e-05,
  maxit = 100,
  constraint_fn,
  maxit_glm = 100,
  ncores = NULL
)
}
\arguments{
\item{formula_rhs}{The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.}

\item{Y}{An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.}

\item{X}{Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.}

\item{covariate_data}{A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.}

\item{z}{Optional initial parameter estimate for z vector.}

\item{tolerance}{The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.}

\item{maxit}{The maximum number of iterations of the coordinate descent algorithm.}

\item{constraint_fn}{A constraint function to make the B matrix identifiable.}

\item{maxit_glm}{The maximum number of iterations when running the glm to update the block of Bj parameters in the coordinate descent algorithm.}

\item{ncores}{The desired number of cores to optimize block of B parameters in parallel. If not provided, an appropriate number will be chosen for your machine.}
}
\value{
A list including values of the log likelihood, the B matrix, and the z vector at each iteration.
}
\description{
Fit Algorithm 2 - DEPRICATED because fit_peanlized_bcd_unconstrained_fast is better
Estimate B and z parameters through block coordinate descent to maximize the Firth penalized log likelihood function, details given in Algorithm 2.
}
\examples{
n <- 75
J <- 5
X <- cbind(1, rnorm(n))
z <- rnorm(n) + 5
b0 <- rnorm(J)
b1 <- rnorm(J)
b <- rbind(b0, b1)
Y <- matrix(NA, ncol = J, nrow = n)

for (i in 1:n) {
 for (j in 1:J) {
   temp_mean <- exp(X[i, , drop = FALSE] \%*\% b[, j, drop = FALSE] + z[i])
   Y[i,j] <- rpois(1, lambda = temp_mean)
 }
}

start <- proc.time()
res1 <- fit_penalized_bcd_unconstrained(Y = Y, X = X, 
          constraint_fn = function(x) {mean(x)}, ncores = 2)
proc.time() - start

}
