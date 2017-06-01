// The objective vector is created by grouping all elements according to their
// ss. If we are supplied all ss and ms simultaneously, it will be possible to
// reduce the computatiolnal load.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Calculates all of the Bernstein basis densities of order \code{m} at
// \code{x}. This is done recursively.
//
// @param x the value where the density is evaluated.
// @param m degree of the Bernstein density.
// @return the Bernstein basis density at x.

arma::vec bernstein_basis_density(const double& x, const int& m) {
  arma::vec bernstein_basis_density_value(m+1);

  // This is a special case where the recursion won't work.
  if (m == 0) {
    bernstein_basis_density_value(0) = 1.0;
    return(bernstein_basis_density_value);
  }

  double y = 1.0 - x;

  bernstein_basis_density_value[0] = y;
  bernstein_basis_density_value[1] = x;

  for (int i = 2; i <= m; i++){
    bernstein_basis_density_value[i] = x*bernstein_basis_density_value[i-1];

    for (int j = i - 1; 1 <= j; j--) {
      bernstein_basis_density_value[j] = x*bernstein_basis_density_value[j-1] + y*bernstein_basis_density_value[j];
    }

    bernstein_basis_density_value[0] = y*bernstein_basis_density_value[0];
  }

  bernstein_basis_density_value = bernstein_basis_density_value*(m+1);

  return bernstein_basis_density_value;
}

//' A function that calculates the objective vector of a polygram. Needed for
//' the quadratic optimization routine, and can be of independent interest.
//'
//'
//' @param data a Armadillo<double> vector containing the data
//' @param ms an Armadillo<integer> vector containing the Bernstein degrees. Must be
//' of size |s| + 1
//' @param s_aug an Armadillo<double> vector containing the splits, endpoints not included.
//' @return an Armadillo<double> of size \code{sum(ms+1)} that contains the
//' objective vector for the optimization program.
// [[Rcpp::export]]
arma::vec polygram_objective_vector(arma::vec data, const arma::Col<int> ms, const arma::vec s){

  unsigned int bin_number = ms.n_elem;
  unsigned int s_length = s.n_elem;

  if((s_length + 1) != (bin_number)) {
    throw std::range_error("The length of s must be equal to that of ms + 1.");
  }

  int data_length = data.n_elem;

  // The purpose of these lines is to construct (0,s,1) from s, with special attention
  // the case when s = numeric(0).

  arma::vec s_aug = arma::zeros<arma::vec>(s_length + 2);
  s_aug(0) = 0;
  s_aug(s_length + 1) = 1;

  if(s_length > 0) {
    s_aug.subvec(1,s_length) = s;
  }

  data = sort(data); // sorted data is required for the algorithm to work.
  arma::vec bernstein_moments = arma::zeros<arma::vec>(sum(ms+1)); // the vector that should be returned.

  // We must loop through the data, keeping track of which bin we are at. When
  // all data beloning to a bin have been identified, it is appended to
  // bernstein_moments, continuing to the next.

  int bin_index = 0;
  double s_lower = s_aug(bin_index);
  double s_upper = s_aug(bin_index + 1);
  int m = ms(bin_index);
  double local_difference = 1/(s_upper - s_lower);
  int lower_bin_index = 0;
  int upper_bin_index = ms(0);

  //std::cout << data_length << std::endl;

  for (int index = 0; index <  data_length; index++) {
      //std::cout << index << std::endl;
      if (data[index] > s_aug[bin_index+1]) {

        bin_index = bin_index + 1;
        s_lower = s_aug(bin_index);
        s_upper = s_aug(bin_index + 1);
        m = ms(bin_index);
        local_difference = 1/(s_upper - s_lower);
        lower_bin_index = upper_bin_index + 1;
        upper_bin_index += ms(bin_index) + 1;

      } else {

        double current_data_value = (data(index) - s_lower)*local_difference;
        bernstein_moments.subvec(lower_bin_index, upper_bin_index) +=
          bernstein_basis_density(current_data_value, m)*local_difference;

      }
  }

  return bernstein_moments/data_length ;

}
