% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reduced_bias_B_alg.R
\name{reduced_bias_B_alg}
\alias{reduced_bias_B_alg}
\title{Get reduced bias B estimates and augmentation values
Implement algorithm from Kosmidis & Firth (2011) to get reduced bias B estimates and augmentation values}
\usage{
reduced_bias_B_alg(X, cc, Y, tolerance = 0.01, maxit = 1000)
}
\arguments{
\item{X}{Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s.}

\item{cc}{Constraint category}

\item{Y}{Data matrix}

\item{tolerance}{The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.}

\item{maxit}{The maximum number of iterations of the coordinate descent algorithm.}

\item{B}{Beta matrix}

\item{z}{z vector}
}
\value{
B estimates and augmentation values
}
\description{
Get reduced bias B estimates and augmentation values
Implement algorithm from Kosmidis & Firth (2011) to get reduced bias B estimates and augmentation values
}
\examples{
n <- 50
J <- 10
X <- cbind(1, rnorm(n))
#X <- cbind(1, rep(0:1, n/2))
b0 <- c(rnorm(J - 1), 0)
b1 <- c(rnorm(J - 1), 0)
B <- rbind(b0, b1)
z_init <- rnorm(n) 
Y <- matrix(NA, ncol = J, nrow = n)

set.seed(2)
for (i in 1:n) {
 for (j in 1:J) {
   temp_mean <- exp(X[i, , drop = FALSE] \%*\% B[, j, drop = FALSE] + z_init[i])
   Y[i,j] <- rpois(1, lambda = temp_mean)
 }
}


res <- reduced_bias_B_alg(X, cc = NULL, Y, tolerance = 1e-6, maxit = 10000)

}
