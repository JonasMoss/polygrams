% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constraints.R
\name{get_symmetry_constraints}
\alias{get_symmetry_constraints}
\title{Returns the constraints needed to enforce symmetry of a Bernstein polygram.}
\usage{
get_symmetry_constraints(ms, s, symmetric = FALSE)
}
\arguments{
\item{s}{The split points. Must be symmetric for the procedure to work.}

\item{m}{The order of the Bernstein densities.}
}
\value{
A constraint matrix.
}
\description{
Returns the constraints needed to enforce symmetry of a Bernstein polygram.
}
