% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\name{get_edge_constraints}
\alias{get_edge_constraints}
\title{Returns the edge constraints. Not finished yet.}
\usage{
get_edge_constraints(ms, s, directions = NULL)
}
\arguments{
\item{ms}{The order of the Bernstein polynomial.}

\item{s}{The vector of splits.}

\item{directions}{Directions for the derivatives. Can be "equal", "below",
or "above" for each derivative in ps.}

\item{ps}{A list of desired values for the derivatives.}
}
\value{
A matrix of constraints.
}
\description{
Returns the edge constraints. Not finished yet.
}
