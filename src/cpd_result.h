#ifndef RELIEVER_CPD_RESULT_H
#define RELIEVER_CPD_RESULT_H

#include <armadillo>
#include <cstddef>
#include <vector>

struct CpdCandidates {
  arma::ivec flat;
  arma::ivec start;
  arma::ivec len;
};

struct SingleCpdResult {
  CpdCandidates cpd;
  arma::vec loss;
  arma::vec objective;
  arma::ivec cps_num;
  arma::vec path_score;
  arma::mat num_pruned;
  int n_model_fit = 0;
  double model_fit_time = 0.0;
  double total_time = 0.0;
};

struct AllCpdResults {
  explicit AllCpdResults(const std::size_t n_run = 0) {
    runs.reserve(n_run);
  }

  std::vector<SingleCpdResult> runs;
};

void set_cpd_candidates_from_matrix(SingleCpdResult & out,
                                    const arma::mat & cpd_mat,
                                    const bool include_baseline);

#endif
