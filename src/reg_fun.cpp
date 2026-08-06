// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "reg_fun.h"
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
arma::vec nmcd_individual_loss(const arma::vec & x,
                               int l,
                               int r,
                               int l_end,
                               int r_end,
                               const arma::vec & sorted_reference,
                               int tail_truncation) {
  const int n = static_cast<int>(x.n_elem);
  const int reference_n = static_cast<int>(sorted_reference.n_elem);
  if (n == 0 || !x.is_finite()) {
    Rcpp::stop("x must contain finite values.");
  }
  if (l < 1 || l > r || r > n ||
      l_end < 1 || l_end > r_end || r_end > n) {
    Rcpp::stop("NMCD interval endpoints are out of range.");
  }
  if (reference_n < 2 || !sorted_reference.is_finite() ||
      !std::is_sorted(sorted_reference.begin(), sorted_reference.end())) {
    Rcpp::stop(
      "sorted_reference must contain at least two finite sorted values."
    );
  }
  if (tail_truncation < 0 ||
      tail_truncation > (reference_n - 2) / 2) {
    Rcpp::stop(
      "tail_truncation leaves too few empirical-CDF cutpoints."
    );
  }

  std::vector<double> training(
    x.begin() + l - 1,
    x.begin() + r
  );
  std::sort(training.begin(), training.end());
  const double training_n = static_cast<double>(training.size());
  const int n_cutpoints = reference_n - 1 - 2 * tail_truncation;

  arma::vec loss_if_greater(n_cutpoints);
  arma::vec loss_if_less_equal(n_cutpoints);
  for (int cutpoint = 0; cutpoint < n_cutpoints; ++cutpoint) {
    const int reference_id = tail_truncation + cutpoint;
    const double threshold = sorted_reference(reference_id);
    double probability = static_cast<double>(
      std::upper_bound(training.begin(), training.end(), threshold) -
      training.begin()
    ) / training_n;
    if (std::abs(probability) <= 1e-15) {
      probability = 1.0 / (2.0 * training_n);
    }
    if (std::abs(probability - 1.0) <= 1e-15) {
      probability = 1.0 - 1.0 / (2.0 * training_n);
    }

    const double weight = static_cast<double>(reference_n) /
      ((reference_id + 1.0) * (reference_n - reference_id - 1.0));
    loss_if_greater(cutpoint) = -std::log1p(-probability) * weight;
    loss_if_less_equal(cutpoint) = -std::log(probability) * weight;
  }

  // Cutpoints below x use -log(1-p); the rest use -log(p). Prefix and
  // suffix sums reduce every evaluation row to one binary search.
  arma::vec greater_prefix(n_cutpoints + 1, arma::fill::zeros);
  arma::vec less_equal_suffix(n_cutpoints + 1, arma::fill::zeros);
  for (int i = 0; i < n_cutpoints; ++i) {
    greater_prefix(i + 1) = greater_prefix(i) + loss_if_greater(i);
  }
  for (int i = n_cutpoints - 1; i >= 0; --i) {
    less_equal_suffix(i) =
      less_equal_suffix(i + 1) + loss_if_less_equal(i);
  }

  const double * cutpoint_begin =
    sorted_reference.begin() + tail_truncation;
  const double * cutpoint_end = cutpoint_begin + n_cutpoints;
  arma::vec loss(r_end - l_end + 1);
  for (int row = l_end - 1; row < r_end; ++row) {
    const int first_less_equal = static_cast<int>(
      std::lower_bound(cutpoint_begin, cutpoint_end, x(row)) - cutpoint_begin
    );
    loss(row - l_end + 1) =
      greater_prefix(first_less_equal) +
      less_equal_suffix(first_less_equal);
  }
  return loss;
}
