#include "shared.hpp"
#include "distributions_bernstein_density.hpp"
#include "distributions_bernstein_polygrams.hpp"
#include <cmath>

arma::vec dpolygram_cpp(const arma::vec& x,
                        const arma::vec& w,
                        const arma::Col<int>& ms,
                        const arma::vec& s_aug,
                        const bool log) {

  unsigned int x_size = x.n_elem;

  arma::vec density_values(x_size);

  unsigned int bin_index = 0;
  unsigned int m = ms(bin_index);
  unsigned int lower_index = 0;
  unsigned int upper_index = m;

  arma::vec lambda = w.subvec(lower_index, upper_index);

  arma::vec2 support;
  support(0) = s_aug(bin_index);
  support(1) = s_aug(bin_index + 1);

  for (unsigned int index = 0; index <  x_size; index++) {
    if (x[index] > s_aug[bin_index+1]) {

      bin_index = bin_index + 1;
      support(0) = s_aug(bin_index);
      support(1) = s_aug(bin_index + 1);

      lower_index += m + 1;
      m = ms(bin_index);
      upper_index += m + 1;

      lambda = w.subvec(lower_index, upper_index);

    }

    arma::vec basis_values = bernstein_basis_density(x[index], m, support);
    density_values(index) = arma::dot(lambda, basis_values);

  }

  return density_values;
}


arma::vec dxpolygram_cpp(const arma::vec& x,
                        const arma::vec& w,
                        const arma::Col<int>& ms,
                        const arma::vec& s_aug,
                        const int& p,
                        const bool log) {

  unsigned int x_size = x.n_elem;

  arma::vec density_values(x_size);

  unsigned int bin_index = 0;
  unsigned int m = ms(bin_index);
  unsigned int lower_index = 0;
  unsigned int upper_index = m;
  double adder = 0;
  double poch;
  double mult = 1;

  arma::mat coefs;
  arma::vec lambda = w.subvec(lower_index, upper_index);
  arma::vec2 support;

  support(0) = s_aug(bin_index);
  support(1) = s_aug(bin_index + 1);

  if(p >= 0) {
    coefs = derivative_coeffient_matrix(m, p);
    poch = boost::math::falling_factorial<double>((m+1), p);
  } else {
    coefs = integral_coeffient_matrix(m, p);
    poch  = 1/boost::math::falling_factorial<double>((m - p + 1), -p);
    mult  = std::pow((support(1) - support(0)),-p);
  }

  arma::vec multiplier = coefs.t()*lambda;

  for (unsigned int index = 0; index <  x_size; index++) {
    if (x[index] > s_aug[bin_index+1]) {


      bin_index = bin_index + 1;
      support(0) = s_aug(bin_index);
      support(1) = s_aug(bin_index + 1);

      if(p >= 0) {
        coefs = derivative_coeffient_matrix(m, p);
        poch = boost::math::falling_factorial<double>((m+1), p);
      } else {
        coefs = integral_coeffient_matrix(m, p);
        poch  = 1/boost::math::falling_factorial<double>((m - p + 1), -p);
        adder = adder + arma::dot(lambda, coefs.col(m - p))*(m-p+1)*poch;
        //std::cout << adder << std::endl;
        mult  = std::pow((support(1) - support(0)),-p);
      }

      lower_index += m + 1;
      m = ms(bin_index);
      upper_index += m + 1;

      lambda = w.subvec(lower_index, upper_index);

      multiplier = coefs.t()*lambda;

      // for(auto& elem:multiplier) {
      //   std::cout << elem << std::endl;
      // }
      // std::cout << "!" << std::endl;
    }

    arma::vec basis_values = bernstein_basis_density((x[index]), m-p, support);
    density_values(index) = arma::dot(multiplier, basis_values)*mult*poch + adder;

  }

  return density_values;
}




