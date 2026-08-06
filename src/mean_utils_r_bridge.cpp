// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "mean_utils.h"

// [[Rcpp::export]]
arma::vec stable_col_means_cpp(const arma::mat & x) {
  return stable_col_means_double(x);
}
