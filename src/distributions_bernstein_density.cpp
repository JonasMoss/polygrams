#include "shared.hpp"
#include "distributions_bernstein_density.hpp"
// [[Rcpp::export]]
arma::mat derivative_coeffient_matrix(const int m,
                                      const int p) {
  const int ncol = m + 1 - p;
  const int nrow = m + 1;

  arma::mat bounds(nrow, ncol, arma::fill::zeros);
  arma::vec multiplier(p+1);

  for(int index = 0; index <= p; index ++){
    double binom = boost::math::binomial_coefficient<double>(p, index);
    int sign = ((index + p) % 2 == 0) ? 1 : -1;
    multiplier(index) = binom*sign;
  }

  for(int nu = 0; nu <  nrow; nu++) {

    int lower_bound = std::max(0, nu - p);
    int upper_bound = std::min(nu, m - p);

    for(int eta = lower_bound; eta <= upper_bound; eta ++) {
      bounds(arma::span(eta,eta + p), eta) = multiplier;
    }

  }

  return bounds;
}
// [[Rcpp::export]]
arma::mat integral_coeffient_matrix(const int m,
                                    const int p) {
  const int ncol = m + 1 - p;
  const int nrow = m + 1;

  arma::mat bounds(nrow, ncol, arma::fill::zeros);
  arma::vec fillings(m+1, arma::fill::zeros);

  for(int index = 0; index <  (m+1); index++) {
    fillings(index) = boost::math::binomial_coefficient<double>(index - p - 1, -p - 1);
  }

  for(int nu = 0; nu < (m+1); nu++) {
      arma::vec fill = fillings.subvec(0, nu);
      std::reverse(fill.begin(), fill.end());
      bounds(arma::span(0, nu), nu - p) = fill;
  }

  return bounds;
}

arma::vec dbernstein_cpp(const arma::Col<double>& x,
                         const arma::Col<double>& lambda,
                         const arma::vec& support,
                         const bool& log) {

  if(!tol_equal(sum(lambda), 1)) {
    Rcpp::warning("The vector lambda must sum to 1.");
  }

  unsigned int x_size = x.n_elem;
  unsigned int m = lambda.n_elem - 1;

  arma::vec density_values(x_size);
  arma::vec basis_values(m);

  for(unsigned int index = 0; index < x_size; index++) {
    basis_values = bernstein_basis_density((x[index]), m, support);
    density_values(index) = arma::dot(lambda, basis_values);
  }

  return density_values;
}

arma::vec pbernstein_cpp(const arma::Col<double>& q,
                         const arma::Col<double>& lambda,
                         const arma::vec& support,
                         const bool& log) {

  if(!tol_equal(sum(lambda), 1)) {
    Rcpp::warning("The vector lambda must sum to 1.");
  }

  unsigned int q_size = q.n_elem;
  unsigned int m = lambda.n_elem - 1;

  arma::vec density_values(q_size);
  arma::vec basis_values(m+1);

  // In order to calculate the integral of the Bernstein densites, we must use
  // the formula that int b_nu,m = \sum(b_nu,m+1, eta = nu + 1 ... m + 1). For
  // this we will require the vector cumsum_lambda.

  arma::vec cumsum_lambda(m+1);
  cumsum_lambda(0) = lambda(0);

  for(unsigned int index = 1; index <= m; index++) {
    cumsum_lambda(index) = cumsum_lambda(index-1) + lambda(index);
  }

  // This is a vectorized way to calculate the final distribution.

  for(unsigned int index = 0; index < q_size; index++) {
    // The last m+1 Bernstein polynomials are required.
    basis_values = bernstein_basis_density(q[index], m+1, support).subvec(1, m+1);

    density_values(index) = arma::dot(cumsum_lambda, basis_values)/(m+2);
  }

  return density_values;
}

arma::vec dxbernstein_cpp(const arma::Col<double>& x,
                          const arma::Col<double>& lambda,
                          const arma::vec& support,
                          const int& p,
                          const bool& log) {

  if(!tol_equal(sum(lambda), 1)) {
    Rcpp::warning("The vector lambda must sum to 1.");
  }

  const int x_size = x.n_elem;
  const int m = lambda.n_elem - 1;
  double poch;
  arma::mat coefs;

  if(p >= 0) {
    coefs = derivative_coeffient_matrix(m, p);
    poch = boost::math::falling_factorial<double>((m+1), p);
  } else {
    coefs = integral_coeffient_matrix(m, p);
    poch = 1/boost::math::falling_factorial<double>((m - p + 1), -p);
  }

  arma::vec multiplier = coefs.t()*lambda;

  arma::vec derivative_values(x_size);
  arma::vec basis_values(m);

  for(int index = 0; index < x_size; index++) {
    basis_values = bernstein_basis_density((x[index]), m-p, support);
    derivative_values(index) = arma::dot(multiplier, basis_values)*poch;
  }

  derivative_values = log ? arma::log(derivative_values) : derivative_values;

  return derivative_values;
}
