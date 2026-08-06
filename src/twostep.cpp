#include "twostep.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

arma::mat suffix_cumsum(const arma::mat & x) {
  return arma::flipud(arma::cumsum(arma::flipud(x), 0));
}

} // namespace

class TwoStepSplitEvaluator : public SplitEvaluator {
public:
  TwoStepSplitEvaluator(TwoStepSearch & search,
                        const int & loss_output)
    : search(search),
      loss_output(loss_output) {}

  BestSplitResult best_split(const int & left,
                             const int & right,
                             const int & dm,
                             const int & interval_slot) override {
    return search.best_split(
      loss_output, left, right, dm, interval_slot
    );
  }

private:
  TwoStepSearch & search;
  int loss_output;
};

TwoStepSearch::TwoStepSearch(
  RegLossFunction & reg_loss,
  const int & n_loss_outputs,
  const int & n_search_intervals,
  const arma::vec & initial_split_fractions
) : reg_loss(reg_loss),
    n_search_intervals(n_search_intervals),
    initial_split_fractions(initial_split_fractions) {
  const int n_slots = n_search_intervals + 1;
  this->gains.set_size(n_loss_outputs, n_slots);
  this->gains.fill(arma::datum::inf);
  this->splits.set_size(n_loss_outputs, n_slots);
  this->splits.fill(-1);
}

arma::ivec TwoStepSearch::unique_initial_splits(
  const int & left,
  const int & right,
  const int & dm
) const {
  const int n_obs = right - left + 1;
  const int first_split = left + dm - 1;
  const int last_split = right - dm;
  std::vector<int> out;
  out.reserve(this->initial_split_fractions.n_elem);
  for (arma::uword i = 0; i < this->initial_split_fractions.n_elem; ++i) {
    const int raw_split = static_cast<int>(
      this->initial_split_fractions(i) * n_obs
    ) + left - 1;
    const int split = std::max(
      first_split, std::min(raw_split, last_split)
    );
    if (std::find(out.begin(), out.end(), split) == out.end()) {
      out.push_back(split);
    }
  }
  return arma::ivec(out);
}

arma::mat TwoStepSearch::fit_individual_loss(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  ++this->model_fit_calls;
  return this->reg_loss.individual_loss(l, r, l_end, r_end);
}

BestSplitResult TwoStepSearch::best_split(
  const int & loss_output,
  const int & left,
  const int & right,
  const int & dm,
  const int & interval_slot
) {
  if (loss_output < 1 || loss_output > static_cast<int>(this->gains.n_rows)) {
    throw std::runtime_error("twostep loss output index is out of range.");
  }
  if (interval_slot < 0 ||
      interval_slot >= static_cast<int>(this->gains.n_cols)) {
    throw std::runtime_error("twostep interval slot is out of range.");
  }

  const bool dynamic_bs_slot = interval_slot >= this->n_search_intervals;
  if (dynamic_bs_slot ||
      !std::isfinite(this->gains(loss_output - 1, interval_slot))) {
    this->compute_interval(left + 1, right, dm, interval_slot);
  }

  BestSplitResult out;
  out.gain = this->gains(loss_output - 1, interval_slot);
  out.tau = this->splits(loss_output - 1, interval_slot);
  return out;
}

void TwoStepSearch::compute_interval(
  const int & left,
  const int & right,
  const int & dm,
  const int & interval_slot
) {
  const int n_obs = right - left + 1;
  if (n_obs < 2 * dm) {
    throw std::runtime_error("the twostep interval is too narrow.");
  }

  const int n_candidates = n_obs - 2 * dm + 1;
  arma::mat best_loss(n_candidates, this->gains.n_rows);
  best_loss.fill(arma::datum::inf);

  arma::ivec initial_splits = this->unique_initial_splits(left, right, dm);
  for (arma::uword i = 0; i < initial_splits.n_elem; ++i) {
    const int split = initial_splits(i);
    arma::mat loss_left = arma::cumsum(
      this->fit_individual_loss(left, split, left, right - dm), 0
    );
    arma::mat loss_right = suffix_cumsum(
      this->fit_individual_loss(split + 1, right, left + dm, right)
    );
    arma::mat candidate_loss =
      loss_left.rows(dm - 1, dm + n_candidates - 2) +
      loss_right.rows(0, n_candidates - 1);
    best_loss = arma::min(best_loss, candidate_loss);
  }

  arma::rowvec full_cost = arma::sum(
    this->fit_individual_loss(left, right, left, right), 0
  );
  for (arma::uword row = 0; row < best_loss.n_cols; ++row) {
    const arma::uword best_idx = best_loss.col(row).index_min();
    this->gains(row, interval_slot) =
      full_cost(row) - best_loss(best_idx, row);
    this->splits(row, interval_slot) = left + dm - 1 + best_idx;
  }
}

SingleCpdResult TwoStepSearch::wbs_one_loss_output(
  const int & n,
  const int & L,
  const int & dm,
  const arma::imat & search_intervals,
  const int & loss_output,
  const bool & recursive
) {
  const int fit_before = this->model_fit_calls;
  const double time_before = this->reg_loss.model_fit_time();
  TwoStepSplitEvaluator split_evaluator(*this, loss_output);
  WbsSearchPath path = recursive ?
    wbs_recursive_search_path(
      n, L, dm, search_intervals, split_evaluator, false
    ) :
    wbs_search_path(n, L, dm, search_intervals, split_evaluator, false);

  SingleCpdResult out;
  set_cpd_candidates_from_matrix(
    out, sorted_cpd_path(path.tau_hat), true
  );
  out.loss.set_size(out.cpd.len.n_elem);
  out.loss.fill(arma::datum::nan);
  out.objective = out.loss;
  out.cps_num.set_size(out.cpd.len.n_elem);
  for (arma::uword i = 0; i < out.cps_num.n_elem; ++i) {
    out.cps_num(i) = static_cast<int>(i);
  }
  out.path_score = path.gain;
  out.n_model_fit = this->model_fit_calls - fit_before;
  out.model_fit_time = this->reg_loss.model_fit_time() - time_before;
  out.total_time = path.total_time;
  return out;
}

const arma::mat & TwoStepSearch::gain_matrix() const {
  return this->gains;
}

const arma::mat & TwoStepSearch::split_matrix() const {
  return this->splits;
}
