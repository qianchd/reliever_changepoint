#ifndef RELIEVER_CPD_ALGORITHMS_H
#define RELIEVER_CPD_ALGORITHMS_H

#include <RcppArmadillo.h>
#include "cost_engine.h"
#include "cpd_result.h"

struct BestSplitResult {
  double gain = -arma::datum::inf;
  int tau = -1;
};

class SplitEvaluator {
public:
  virtual ~SplitEvaluator() {}
  virtual BestSplitResult best_split(const int & left,
                                     const int & right,
                                     const int & dm,
                                     const int & interval_slot) = 0;
};

struct WbsSearchPath {
  arma::vec tau_hat;
  arma::vec gain;
  double total_time = 0.0;
};

enum class CpdSearchMethod {
  SegmentNeighbourhood,
  Wbs,
  WbsRecursive,
  Pelt
};

arma::mat sorted_cpd_path(const arma::vec & tau_hat);

SingleCpdResult sn_one_loss_output(const int & n,
                                   const int & L,
                                   const int & dm,
                                   CostEngine & cost_engine);

WbsSearchPath wbs_search_path(const int n,
                              const int L,
                              const int dm,
                              const arma::imat & lr_m,
                              SplitEvaluator & split_evaluator,
                              const bool & prefer_last_tie = true);

WbsSearchPath wbs_recursive_search_path(const int n,
                                        const int L,
                                        const int dm,
                                        const arma::imat & lr_m,
                                        SplitEvaluator & split_evaluator,
                                        const bool & prefer_last_tie = true);

SingleCpdResult wbs_one_loss_output(const int n,
                                    const int L,
                                    const int dm,
                                    const arma::imat & lr_m,
                                    CostEngine & cost_engine,
                                    const bool & recursive = false);

SingleCpdResult pelt_one_loss_output(const int & n,
                                     const arma::vec & pen_val,
                                     const int & dm,
                                     const double & prune_value,
                                     CostEngine & cost_engine);

SingleCpdResult cpd_one_loss_output(
  const CpdSearchMethod method,
  const int & n,
  const int & L,
  const int & dm,
  const arma::imat & search_intervals,
  const arma::vec & pen_val,
  const double & prune_value,
  CostEngine & cost_engine
);

#endif
