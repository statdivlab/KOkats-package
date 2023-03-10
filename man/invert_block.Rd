% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invert_block.R
\name{invert_block}
\alias{invert_block}
\title{Invert matrix
Invert a block of the block diagonal information matrix for the beta parameters.}
\usage{
invert_block(j, p, info)
}
\arguments{
\item{j}{Category j.}

\item{p}{Number of covariates in the model.}

\item{info}{Full sparse information matrix for beta.}
}
\value{
The inverse of the block of the info matrix corresponding with category j.
}
\description{
Invert matrix
Invert a block of the block diagonal information matrix for the beta parameters.
}
