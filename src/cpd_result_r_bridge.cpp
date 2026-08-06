#include "cpd_result_r_bridge.h"

namespace {

struct FlatCpdBridge {
  explicit FlatCpdBridge(const int n_run) : num_pruned(n_run) {}

  std::vector<int> cpd_flat;
  std::vector<int> cpd_start;
  std::vector<int> cpd_len;
  std::vector<int> run_start;
  std::vector<int> run_len;
  std::vector<double> loss;
  std::vector<double> objective;
  std::vector<int> cps_num;
  std::vector<double> path_score_flat;
  std::vector<int> path_score_start;
  std::vector<int> path_score_len;
  std::vector<double> n_model_fit;
  std::vector<double> model_fit_time;
  std::vector<double> total_time;
  Rcpp::List num_pruned;
};

Rcpp::IntegerVector to_integer_vector(const std::vector<int> & x) {
  Rcpp::IntegerVector out(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    out[i] = x[i];
  }
  return out;
}

Rcpp::NumericVector to_numeric_vector(const std::vector<double> & x) {
  Rcpp::NumericVector out(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    out[i] = x[i];
  }
  return out;
}

void append_run(FlatCpdBridge & flat,
                const SingleCpdResult & run,
                const int run_id) {
  flat.run_start.push_back(static_cast<int>(flat.cpd_len.size()));
  flat.run_len.push_back(static_cast<int>(run.cpd.len.n_elem));

  for (arma::uword j = 0; j < run.cpd.len.n_elem; ++j) {
    const int start = run.cpd.start(j);
    const int len = run.cpd.len(j);
    flat.cpd_start.push_back(static_cast<int>(flat.cpd_flat.size()));
    flat.cpd_len.push_back(len);
    for (int k = 0; k < len; ++k) {
      flat.cpd_flat.push_back(run.cpd.flat(start + k));
    }
  }

  for (arma::uword j = 0; j < run.cpd.len.n_elem; ++j) {
    flat.loss.push_back(j < run.loss.n_elem ? run.loss(j) : NA_REAL);
    flat.objective.push_back(
      j < run.objective.n_elem ? run.objective(j) : NA_REAL
    );
    flat.cps_num.push_back(
      j < run.cps_num.n_elem ? run.cps_num(j) : NA_INTEGER
    );
  }

  flat.path_score_start.push_back(
    static_cast<int>(flat.path_score_flat.size())
  );
  flat.path_score_len.push_back(static_cast<int>(run.path_score.n_elem));
  for (arma::uword j = 0; j < run.path_score.n_elem; ++j) {
    flat.path_score_flat.push_back(run.path_score(j));
  }

  flat.n_model_fit.push_back(run.n_model_fit);
  flat.model_fit_time.push_back(run.model_fit_time);
  flat.total_time.push_back(run.total_time);

  flat.num_pruned[run_id] = run.num_pruned.n_elem > 0 ?
    Rcpp::wrap(run.num_pruned) : R_NilValue;
}

Rcpp::List flat_bridge_to_list(const FlatCpdBridge & flat) {
  return Rcpp::List::create(
    Rcpp::Named("cpd_flat") = to_integer_vector(flat.cpd_flat),
    Rcpp::Named("cpd_start") = to_integer_vector(flat.cpd_start),
    Rcpp::Named("cpd_len") = to_integer_vector(flat.cpd_len),
    Rcpp::Named("run_start") = to_integer_vector(flat.run_start),
    Rcpp::Named("run_len") = to_integer_vector(flat.run_len),
    Rcpp::Named("loss") = to_numeric_vector(flat.loss),
    Rcpp::Named("objective") = to_numeric_vector(flat.objective),
    Rcpp::Named("cps_num") = to_integer_vector(flat.cps_num),
    Rcpp::Named("path_score_flat") = to_numeric_vector(flat.path_score_flat),
    Rcpp::Named("path_score_start") =
      to_integer_vector(flat.path_score_start),
    Rcpp::Named("path_score_len") = to_integer_vector(flat.path_score_len),
    Rcpp::Named("n_model_fit") = to_numeric_vector(flat.n_model_fit),
    Rcpp::Named("model_fit_time") = to_numeric_vector(flat.model_fit_time),
    Rcpp::Named("total_time") = to_numeric_vector(flat.total_time),
    Rcpp::Named("num_pruned") = flat.num_pruned
  );
}

} // namespace

Rcpp::List all_cpd_results_to_list(const AllCpdResults & results) {
  const int n_run = static_cast<int>(results.runs.size());
  FlatCpdBridge flat(n_run);
  for (int i = 0; i < n_run; i++) {
    append_run(flat, results.runs[i], i);
  }
  return flat_bridge_to_list(flat);
}
