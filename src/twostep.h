#ifndef RELIEVER_TWOSTEP_H
#define RELIEVER_TWOSTEP_H

#include "cost_engine.h"
#include "cpd_algorithms.h"

class TwoStepSplitEvaluator;

// Pure C++ two-step search state. R callbacks and native losses both implement
// RegLossFunction before reaching this class.
class TwoStepSearch {
public:
  TwoStepSearch(RegLossFunction & reg_loss,
                const int & n_loss_outputs,
                const int & n_search_intervals,
                const arma::vec & initial_split_fractions);

  SingleCpdResult wbs_one_loss_output(const int & n,
                                      const int & L,
                                      const int & dm,
                                      const arma::imat & search_intervals,
                                      const int & loss_output,
                                      const bool & recursive = false);

  const arma::mat & gain_matrix() const;
  const arma::mat & split_matrix() const;

private:
  friend class TwoStepSplitEvaluator;

  RegLossFunction & reg_loss;
  int n_search_intervals;
  arma::vec initial_split_fractions;
  arma::mat gains;
  arma::mat splits;
  int model_fit_calls = 0;

  arma::ivec unique_initial_splits(const int & left,
                                   const int & right,
                                   const int & dm) const;
  arma::mat fit_individual_loss(const int & l,
                                const int & r,
                                const int & l_end,
                                const int & r_end);
  BestSplitResult best_split(const int & loss_output,
                             const int & left,
                             const int & right,
                             const int & dm,
                             const int & interval_slot);
  void compute_interval(const int & left,
                        const int & right,
                        const int & dm,
                        const int & interval_slot);
};

#endif
