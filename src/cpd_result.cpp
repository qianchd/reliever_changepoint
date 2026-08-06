// Keep the implementation compiled with the same Armadillo configuration used
// by the RcppArmadillo translation units that pass arma objects across here.
#include <RcppArmadillo.h>
#include "cpd_result.h"
#include <cmath>

void set_cpd_candidates_from_matrix(SingleCpdResult & out,
                                    const arma::mat & cpd_mat,
                                    const bool include_baseline) {
  const int n_candidate =
    static_cast<int>(cpd_mat.n_rows) + static_cast<int>(include_baseline);
  std::vector<int> flat;
  std::vector<int> start;
  std::vector<int> len;
  start.reserve(n_candidate);
  len.reserve(n_candidate);

  if (include_baseline) {
    start.push_back(0);
    len.push_back(0);
  }
  for (arma::uword row = 0; row < cpd_mat.n_rows; ++row) {
    start.push_back(static_cast<int>(flat.size()));
    int row_len = 0;
    for (arma::uword col = 0; col < cpd_mat.n_cols; ++col) {
      const double value = cpd_mat(row, col);
      if (std::isfinite(value)) {
        flat.push_back(static_cast<int>(value));
        ++row_len;
      }
    }
    len.push_back(row_len);
  }

  out.cpd.flat.set_size(flat.size());
  for (arma::uword i = 0; i < flat.size(); ++i) {
    out.cpd.flat(i) = flat[i];
  }

  out.cpd.start.set_size(start.size());
  out.cpd.len.set_size(len.size());
  for (arma::uword i = 0; i < start.size(); ++i) {
    out.cpd.start(i) = start[i];
    out.cpd.len(i) = len[i];
  }
}
