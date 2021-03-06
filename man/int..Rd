% Generated by roxygen2 (4.0.1): do not edit by hand
\name{int.}
\alias{int.}
\title{Population-level regression intercept}
\usage{
int.(y, x, p)
}
\arguments{
\item{x}{A numeric vector giving the independent variable of
the population.}

\item{y}{A numeric vector giving the dependent variable of the
population.}

\item{p}{An optional vector of weights.  If missing, equal
weights are used.}
}
\value{
The intercept associated with \code{x} and \code{y}.
}
\description{
Intercept in the regression of one vector on another in a bivariate population.
}

