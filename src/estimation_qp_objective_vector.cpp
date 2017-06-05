// The objective vector is created by grouping all elements according to their
// ss. If we are supplied all ss and ms simultaneously, it will be possible to
// reduce the computatiolnal load.

#include "shared.hpp"
#include "estimation_qp_objective_vector.hpp"

//' A function that calculates the objective vector of a polygram. Needed for
//' the quadratic optimization routine, and can be of independent interest.
//'
//'
//' @param x a Armadillo<double> vector containing the data
//' @param ms an Armadillo<integer> vector containing the Bernstein degrees. Must be
//' of size |s| + 1
//' @param s_aug an Armadillo<double> vector containing the splits, endpoints not included.
//' @return an Armadillo<double> of size \code{sum(ms+1)} that contains the
//' objective vector for the optimization program.
// [[Rcpp::export]]
arma::vec polygram_objective_vector(arma::vec x,
                                    const arma::Col<int>& ms,
                                    const arma::vec& s,
                                    const arma::vec& support){

  const unsigned int bin_number = ms.n_elem;
  const unsigned int s_length = s.n_elem;

  if((s_length + 1) != (bin_number)) {
    throw std::range_error("The length of s must be equal to that of ms + 1.");
  }

  const int x_length = x.n_elem;

  // The purpose of these lines is to construct (0,s,1) from s, with special attention
  // the case when s = numeric(0).

  arma::vec s_aug = arma::zeros<arma::vec>(s_length + 2);
  s_aug(0) = 0;
  s_aug(s_length + 1) = 1;

  if(s_length > 0) {
    s_aug.subvec(1,s_length) = s;
  }

  x = sort(x); // sorted x is required for the algorithm to work.
  arma::vec bernstein_moments = arma::zeros<arma::vec>(sum(ms+1)); // the vector that should be returned.

  // We must loop through the x, keeping track of which bin we are at. When
  // all x beloning to a bin have been identified, it is appended to
  // bernstein_moments, continuing to the next.

  int bin_index = 0;
  double s_lower = s_aug(bin_index);
  double s_upper = s_aug(bin_index + 1);
  int m = ms(bin_index);
  double local_difference = 1/(s_upper - s_lower);
  int lower_bin_index = 0;
  int upper_bin_index = ms(0);

  //std::cout << x_length << std::endl;

  for (int index = 0; index <  x_length; index++) {
      //std::cout << index << std::endl;
      if (x[index] > s_aug[bin_index+1]) {

        bin_index = bin_index + 1;
        s_lower = s_aug(bin_index);
        s_upper = s_aug(bin_index + 1);
        m = ms(bin_index);
        local_difference = 1/(s_upper - s_lower);
        lower_bin_index = upper_bin_index + 1;
        upper_bin_index += ms(bin_index) + 1;

      }

      double current_x_value = (x(index) - s_lower)*local_difference;
      bernstein_moments.subvec(lower_bin_index, upper_bin_index) +=
      bernstein_basis_density(current_x_value, m, support)*local_difference;

  }

  return bernstein_moments/x_length ;

}
