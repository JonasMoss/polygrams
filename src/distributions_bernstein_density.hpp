#ifndef DISTRIBUTIONS_BERNSTEIN_DENSITY
#define DISTRIBUTIONS_BERNSTEIN_DENSITY

#include "shared.hpp"
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>

// Helper functions.

arma::mat derivative_coeffient_matrix(const int m,
                                      const int p);

arma::mat integral_coeffient_matrix(const int m,
                                    const int p);


// Exported functions.
//[[Rcpp::export(rng = false)]]
arma::vec dbernstein_cpp(const arma::Col<double>& x,
                         const arma::Col<double>& lambda,
                         const arma::vec& support,
                         const bool& log = false);

//[[Rcpp::export(rng = false)]]
arma::vec pbernstein_cpp(const arma::Col<double>& q,
                         const arma::Col<double>& lambda,
                         const arma::vec& support,
                         const bool& log = false);

//[[Rcpp::export(rng = false)]]
arma::vec dxbernstein_cpp(const arma::Col<double>& x,
                          const arma::Col<double>& lambda,
                          const arma::vec& support,
                          const int& p,
                          const bool& log = false);

#endif
