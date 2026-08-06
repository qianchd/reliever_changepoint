// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "cost_engine.h"
#include "cost_engine_r_bridge.h"
#include "cpd_algorithms.h"
#include "cpd_result_r_bridge.h"
#include <string>
#include <utility>

using namespace Rcpp;

namespace {

void require_loss_block_interval_set(
  const arma::ivec & miss_cover_len,
  const arma::ivec & int_len,
  const arma::ivec & layer_point,
  const arma::imat & int_eps,
  const char * caller
) {
  if (miss_cover_len.n_elem == 0 || int_len.n_elem == 0 ||
      layer_point.n_elem == 0 || int_eps.n_elem == 0) {
    Rcpp::stop(
      "%s requires a non-full Reliever interval set.",
      caller
    );
  }
}

Rcpp::List loss_block_results_to_list(
  const AllCpdResults & results,
  LossBlockCache & loss_block_cache,
  const bool & return_cache_profile
) {
  warn_loss_block_exact_fallbacks(loss_block_cache);

  Rcpp::List out = all_cpd_results_to_list(results);
  Rcpp::RObject cache_profile = R_NilValue;
  if (return_cache_profile) {
    cache_profile = Rcpp::List::create(
      Rcpp::Named("cache_state") = loss_block_cache_state(loss_block_cache)
    );
  }
  out["cache_profile"] = cache_profile;
  return out;
}

void restore_loss_block_cache(
  LossBlockCache & loss_block_cache,
  Rcpp::Nullable<Rcpp::List> cache_state
) {
  if (cache_state.isNull()) {
    return;
  }
  Rcpp::List cache_state_list(cache_state);
  loss_block_cache.restore(loss_block_cache_state_from_r(cache_state_list));
}

CpdSearchMethod cpd_search_method_from_r(const std::string & method) {
  if (method == "SN") {
    return CpdSearchMethod::SegmentNeighbourhood;
  }
  if (method == "WBS_recursive") {
    return CpdSearchMethod::WbsRecursive;
  }
  if (method == "WBS" || method == "SeedBS" || method == "BS") {
    return CpdSearchMethod::Wbs;
  }
  if (method == "PELT" || method == "OP") {
    return CpdSearchMethod::Pelt;
  }
  Rcpp::stop("Unsupported changepoint search method: %s.", method);
  return CpdSearchMethod::SegmentNeighbourhood;
}

AllCpdResults search_loss_outputs_by_loss_block(
  const CpdSearchMethod method,
  const int & n,
  const int & L,
  const int & dm,
  const arma::imat & search_intervals,
  const arma::vec & pen_val,
  const double & prune_value,
  const arma::ivec & run_loss_outputs,
  LossBlockCache & loss_block_cache
) {
  const int n_run = run_loss_outputs.n_elem;
  AllCpdResults results(n_run);
  double relief_updates_before = loss_block_cache.relief_update_calls;
  double full_updates_before = loss_block_cache.full_update_calls;
  double model_fit_time_before = loss_block_cache.model_fit_time;

  for (int run_id = 0; run_id < n_run; ++run_id) {
    Rcpp::checkUserInterrupt();
    CostEngineByLossBlock cost_engine(
      run_loss_outputs(run_id), loss_block_cache
    );
    SingleCpdResult result = cpd_one_loss_output(
      method, n, L, dm, search_intervals, pen_val, prune_value, cost_engine
    );
    result.n_model_fit = static_cast<int>(
      loss_block_cache.relief_update_calls - relief_updates_before +
      loss_block_cache.full_update_calls - full_updates_before
    );
    result.model_fit_time =
      loss_block_cache.model_fit_time - model_fit_time_before;
    results.runs.push_back(std::move(result));

    relief_updates_before = loss_block_cache.relief_update_calls;
    full_updates_before = loss_block_cache.full_update_calls;
    model_fit_time_before = loss_block_cache.model_fit_time;
  }
  return results;
}

AllCpdResults search_loss_outputs_by_cost_mat(
  const CpdSearchMethod method,
  const int & n,
  const int & L,
  const int & dm,
  const arma::imat & search_intervals,
  const arma::vec & pen_val,
  const double & prune_value,
  const arma::ivec & run_loss_outputs,
  CostMatCache & cost_mat_cache,
  RegLossFunction & reg_loss,
  const bool & is_full
) {
  const int n_run = run_loss_outputs.n_elem;
  AllCpdResults results(n_run);

  for (int run_id = 0; run_id < n_run; ++run_id) {
    Rcpp::checkUserInterrupt();
    CostEngineByCostMat cost_engine(
      cost_mat_cache, reg_loss, run_loss_outputs(run_id), is_full
    );
    SingleCpdResult result = cpd_one_loss_output(
      method, n, L, dm, search_intervals, pen_val, prune_value, cost_engine
    );
    results.runs.push_back(std::move(result));
  }
  return results;
}

arma::ivec one_loss_output() {
  arma::ivec out(1);
  out(0) = 1;
  return out;
}

} // namespace

//' R loss-function backend using the loss-block cache.
//' @noRd
// [[Rcpp::export]]
List cpd_r_by_loss_block(const std::string & method,
                         const int & n,
                         const int & L,
                         const int & dm,
                         const arma::imat & search_intervals,
                         const arma::vec & pen_val,
                         const double & prune_value,
                         const arma::ivec & run_loss_outputs,
                         Rcpp::Function & reg_fun_wrap,
                         const arma::ivec & miss_cover_len,
                         const arma::ivec & int_len,
                         const arma::ivec & layer_point,
                         const arma::imat & int_eps,
                         Rcpp::Nullable<Rcpp::List> cache_state = R_NilValue,
                         bool return_cache_profile = true,
                         bool use_owner_key = true) {
  RRegLossFunction reg_loss(reg_fun_wrap);
  LossBlockCache loss_block_cache(
    reg_loss, miss_cover_len, int_len, layer_point, int_eps, n, use_owner_key
  );
  restore_loss_block_cache(loss_block_cache, cache_state);
  AllCpdResults results = search_loss_outputs_by_loss_block(
    cpd_search_method_from_r(method), n, L, dm, search_intervals,
    pen_val, prune_value, run_loss_outputs, loss_block_cache
  );
  return loss_block_results_to_list(
    results, loss_block_cache, return_cache_profile
  );
}

//' Native mean-loss backend using the loss-block cache.
//' @noRd
// [[Rcpp::export]]
List cpd_mean_by_loss_block_cpp(const std::string & method,
                                const int & n,
                                const int & L,
                                const int & dm,
                                const arma::imat & search_intervals,
                                const arma::vec & pen_val,
                                const double & prune_value,
                                const arma::mat & data,
                                const double & ratio,
                                const arma::ivec & miss_cover_len,
                                const arma::ivec & int_len,
                                const arma::ivec & layer_point,
                                const arma::imat & int_eps,
                                bool return_cache_profile = true,
                                bool use_owner_key = true) {
  require_loss_block_interval_set(
    miss_cover_len, int_len, layer_point, int_eps,
    "cpd_mean_by_loss_block_cpp"
  );
  MeanSquareRegLossFunction reg_loss(ratio, data);
  LossBlockCache loss_block_cache(
    reg_loss, miss_cover_len, int_len, layer_point, int_eps, n, use_owner_key
  );
  AllCpdResults results = search_loss_outputs_by_loss_block(
    cpd_search_method_from_r(method), n, L, dm, search_intervals,
    pen_val, prune_value, one_loss_output(), loss_block_cache
  );
  return loss_block_results_to_list(
    results, loss_block_cache, return_cache_profile
  );
}

//' R loss-function backend using the cost-matrix cache.
//' @noRd
// [[Rcpp::export]]
List cpd_r_by_cost_mat(const std::string & method,
                       const int & n,
                       const int & L,
                       const int & dm,
                       const arma::imat & search_intervals,
                       const arma::vec & pen_val,
                       const double & prune_value,
                       const arma::ivec & run_loss_outputs,
                       arma::mat & cost_mat,
                       Rcpp::Function & reg_fun_wrap,
                       const arma::ivec & miss_cover_len,
                       const arma::ivec & int_len,
                       const arma::ivec & layer_point,
                       const arma::imat & int_eps,
                       bool is_full = false) {
  RRegLossFunction reg_loss(reg_fun_wrap);
  CostMatCache cost_mat_cache(
    cost_mat, miss_cover_len, int_len, layer_point, int_eps, n
  );
  AllCpdResults results = search_loss_outputs_by_cost_mat(
    cpd_search_method_from_r(method), n, L, dm, search_intervals,
    pen_val, prune_value, run_loss_outputs, cost_mat_cache, reg_loss, is_full
  );
  return all_cpd_results_to_list(results);
}

//' Native mean-loss backend using the cost-matrix cache.
//' @noRd
// [[Rcpp::export]]
List cpd_mean_by_cost_mat_cpp(const std::string & method,
                              const int & n,
                              const int & L,
                              const int & dm,
                              const arma::imat & search_intervals,
                              const arma::vec & pen_val,
                              const double & prune_value,
                              const arma::mat & data,
                              arma::mat & cost_mat,
                              const double & ratio,
                              const arma::ivec & miss_cover_len,
                              const arma::ivec & int_len,
                              const arma::ivec & layer_point,
                              const arma::imat & int_eps,
                              bool is_full = false) {
  MeanSquareRegLossFunction reg_loss(ratio, data);
  CostMatCache cost_mat_cache(
    cost_mat, miss_cover_len, int_len, layer_point, int_eps, n
  );
  AllCpdResults results = search_loss_outputs_by_cost_mat(
    cpd_search_method_from_r(method), n, L, dm, search_intervals,
    pen_val, prune_value, one_loss_output(), cost_mat_cache, reg_loss, is_full
  );
  return all_cpd_results_to_list(results);
}
