#include "mean_utils.h"
#include <stdexcept>

arma::vec stable_col_means_double(const arma::mat & x) {
  if (x.n_rows == 0) {
    throw std::invalid_argument("x must contain at least one row.");
  }

  const arma::rowvec origin = x.row(0);
  arma::mat centered = x;
  centered.each_row() -= origin;
  return origin.t() + arma::mean(centered, 0).t();
}
