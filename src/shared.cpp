#include "shared.hpp"

bool tol_equal(const double& x, const double& y) {
  return std::abs(x - y) < MIN_DIFF_EPS;
}

// Calculates all of the Bernstein basis densities of order \code{m} at
// \code{x}. This is done recursively.
//
// @param x the value where the density is evaluated.
// @param m degree of the Bernstein density.
// @return the Bernstein basis density at x.

arma::vec bernstein_basis_density(const double& x,
                                  const int& m,
                                  const arma::vec& support) {
  arma::vec bernstein_basis_density_value(m+1);

  // This is a special case where the recursion won't work.
  if (m == 0) {
    bernstein_basis_density_value(0) = 1.0;
    return(bernstein_basis_density_value*1/(support(1) - support(0)));
  }


  double denom = 1/(support(1) - support(0));
  double z = (x - support(0))*denom;
  double y = 1.0 - z;

  bernstein_basis_density_value[0] = y;
  bernstein_basis_density_value[1] = z;

  for (int i = 2; i <= m; i++){
    bernstein_basis_density_value[i] = z*bernstein_basis_density_value[i-1];

    for (int j = i - 1; 1 <= j; j--) {
      bernstein_basis_density_value[j] = z*bernstein_basis_density_value[j-1] + y*bernstein_basis_density_value[j];
    }

    bernstein_basis_density_value[0] = y*bernstein_basis_density_value[0];
  }

  bernstein_basis_density_value = bernstein_basis_density_value*(m+1);

  return bernstein_basis_density_value*denom;
}
