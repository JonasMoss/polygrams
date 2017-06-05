#ifndef SHARED
#define SHARED

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>

const double MIN_DIFF_EPS = 1e-8;

bool tol_equal(const double& x, const double& y);

arma::vec bernstein_basis_density(const double& x,
                                  const int& m,
                                  const arma::vec& support);

#endif
