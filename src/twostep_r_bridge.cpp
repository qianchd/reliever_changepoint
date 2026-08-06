// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "cost_engine_r_bridge.h"
#include "cpd_result_r_bridge.h"
#include "twostep.h"
#include <utility>

//' Run WBS-family twostep search for a set of loss outputs.
//' @noRd
// [[Rcpp::export]]
Rcpp::List wbs_r_twostep_loss_outputs(const int & n,
                                      const int & L,
                                      const int & dm,
                                      const arma::imat & lr_m,
                                      const arma::ivec & run_loss_outputs,
                                      Rcpp::Function & individual_loss_fun,
                                      const int & n_loss_outputs,
                                      const arma::vec & init_cand,
                                      bool recursive = false) {
  RRegLossFunction reg_loss(individual_loss_fun);
  TwoStepSearch search(
    reg_loss, n_loss_outputs, lr_m.n_rows, init_cand
  );

  const int n_run = run_loss_outputs.n_elem;
  AllCpdResults results(n_run);
  for (int run_id = 0; run_id < n_run; ++run_id) {
    Rcpp::checkUserInterrupt();
    SingleCpdResult cpd_res = search.wbs_one_loss_output(
      n, L, dm, lr_m, run_loss_outputs(run_id), recursive
    );
    results.runs.push_back(std::move(cpd_res));
  }

  Rcpp::List out = all_cpd_results_to_list(results);
  out["gain_mat"] = search.gain_matrix();
  out["split_mat"] = search.split_matrix();
  return out;
}
