#ifndef DISTRIBUTIONS_BERNSTEIN_POLYGRAMS
#define DISTRIBUTIONS_BERNSTEIN_POLYGRAMS

// [[Rcpp::export]]
arma::vec dxpolygram_cpp(const arma::vec& x,
                         const arma::vec& w,
                         const arma::Col<int>& ms,
                         const arma::vec& s_aug,
                         const int& p,
                         const bool log = false);

#endif
