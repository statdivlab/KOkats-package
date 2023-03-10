% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_X_tilde.R
\name{generate_X_tilde}
\alias{generate_X_tilde}
\title{Generate X tilde
Create expanded design matrix X tilde such that log expected value of Y tilde is equal to X tilde times theta.}
\usage{
generate_X_tilde(X, J, transpose = FALSE)
}
\arguments{
\item{X}{Matrix of X values, n rows by p columns}

\item{J}{The number of KOs.}

\item{transpose}{If true, will generate the transpose of X_tilde. Set to \code{FALSE}.}
}
\value{
X_tilde, a matrix with \code{nJ} rows and \code{pJ + n} columns.
}
\description{
Generate X tilde
Create expanded design matrix X tilde such that log expected value of Y tilde is equal to X tilde times theta.
}
