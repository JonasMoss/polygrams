#ifndef ESTIMATION_QP_OBJECTIVE_VECTOR
#define ESTIMATION_QP_OBJECTIVE_VECTOR

#include "shared.hpp"

const arma::vec initial_support = {0, 1};

arma::vec polygram_objective_vector(arma::vec x,
                                    const arma::Col<int>& ms,
                                    const arma::vec& s,
                                    const arma::vec& support = {0, 1});


#endif
