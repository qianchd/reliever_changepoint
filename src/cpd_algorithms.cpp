#include "cpd_algorithms.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

class ExactSplitEvaluator : public SplitEvaluator {
public:
  explicit ExactSplitEvaluator(CostEngine & cost_engine)
    : cost_engine(cost_engine) {}

  BestSplitResult best_split(const int & left,
                             const int & right,
                             const int & dm,
                             const int & interval_slot) override {
    (void) interval_slot;
    const int n_split = right - left - 2 * dm + 1;
    arma::vec split_points(n_split);
    arma::vec gain(n_split);
    split_points.fill(arma::datum::nan);
    gain.fill(arma::datum::nan);

    for (int t = left + dm; t <= right - dm; t++) {
      split_points(t - left - dm) = t;
      gain(t - left - dm) = -this->cost_engine.get_cost(left + 1, t);
    }
    const double full_cost = this->cost_engine.get_cost(left + 1, right);
    for (int t = right - dm; t >= left + dm; t--) {
      gain(t - left - dm) +=
        full_cost - this->cost_engine.get_cost(t + 1, right);
    }

    BestSplitResult out;
    arma::uvec finite_idx = arma::find_finite(gain);
    if (finite_idx.n_elem == 0) {
      return out;
    }

    arma::vec finite_gain = gain.elem(finite_idx);
    arma::uword best_pos = finite_gain.index_max();
    double best_gain = finite_gain(best_pos);
    arma::uvec ties = arma::find(
      arma::abs(finite_gain - best_gain) < 1e-8, 1, "last"
    );
    arma::uword best_idx = finite_idx(ties(0));

    out.gain = gain(best_idx);
    out.tau = split_points(best_idx);
    return out;
  }

private:
  CostEngine & cost_engine;
};

arma::vec finite_inner_changepoints(const arma::vec & tau_hat_all,
                                    const int & n) {
  arma::vec tau_hat_finite = tau_hat_all(arma::find_finite(tau_hat_all));
  return tau_hat_finite(arma::find(
    (tau_hat_finite > 0) && (tau_hat_finite < n)
  ));
}

void set_cost_engine_metrics(SingleCpdResult & result,
                             CostEngine & cost_engine,
                             const double & total_time) {
  result.n_model_fit = cost_engine.get_n_model_fit();
  result.model_fit_time = cost_engine.get_model_fit_time();
  result.total_time = total_time;
}

struct SegmentSplit {
  int left = 0;
  int right = 0;
  int interval_col = -1;
  double gain = arma::datum::nan;
  int tau = -1;

  bool finite() const {
    return std::isfinite(gain) && tau > left && tau < right;
  }
};

SegmentSplit best_split_for_segment(const int & left,
                                    const int & right,
                                    const int & row_id,
                                    const int & M,
                                    const int & dm,
                                    const arma::imat & lr_m,
                                    SplitEvaluator & split_evaluator,
                                    const bool & prefer_last_tie) {
  SegmentSplit best;
  best.left = left;
  best.right = right;

  for (int i = 0; i <= M; i++) {
    int li, ri;
    int interval_slot;
    if (i == M) {
      li = left;
      ri = right;
      interval_slot = M + row_id;
    } else {
      li = lr_m(i, 0);
      ri = lr_m(i, 1);
      interval_slot = i;
    }

    if (left <= li && ri <= right && ri - li >= 2 * dm) {
      const BestSplitResult split =
        split_evaluator.best_split(li, ri, dm, interval_slot);
      if (!std::isfinite(split.gain)) {
        continue;
      }

      const bool better_gain = !best.finite() || split.gain > best.gain;
      const bool tied_gain =
        best.finite() && std::abs(split.gain - best.gain) < 1e-8;
      const bool better_tie = prefer_last_tie && tied_gain &&
        i > best.interval_col;

      if (better_gain || better_tie) {
        best.interval_col = i;
        best.gain = split.gain;
        best.tau = split.tau;
      }
    }
  }

  return best;
}

bool segment_order(const SegmentSplit & a, const SegmentSplit & b) {
  if (a.left != b.left) {
    return a.left < b.left;
  }
  return a.right < b.right;
}

arma::vec selected_tau_vector(const std::vector<int> & tau_hat) {
  arma::vec out(tau_hat.size());
  for (arma::uword i = 0; i < out.n_elem; ++i) {
    out(i) = tau_hat[i];
  }
  return out;
}

arma::vec selected_gain_vector(const std::vector<double> & gain) {
  arma::vec out(gain.size());
  for (arma::uword i = 0; i < out.n_elem; ++i) {
    out(i) = gain[i];
  }
  return out;
}

void append_recursive_wbs_path(const int & left,
                               const int & right,
                               const int & L,
                               const int & M,
                               const int & dm,
                               const arma::imat & lr_m,
                               SplitEvaluator & split_evaluator,
                               const bool & prefer_last_tie,
                               std::vector<int> & tau_hat_all,
                               std::vector<double> & gain_all) {
  if (static_cast<int>(tau_hat_all.size()) >= L) {
    return;
  }

  SegmentSplit selected = best_split_for_segment(
    left, right, 0, M, dm, lr_m, split_evaluator, prefer_last_tie
  );
  if (!selected.finite()) {
    return;
  }

  tau_hat_all.push_back(selected.tau);
  gain_all.push_back(selected.gain);

  append_recursive_wbs_path(
    selected.left, selected.tau, L, M, dm, lr_m, split_evaluator,
    prefer_last_tie, tau_hat_all, gain_all
  );
  append_recursive_wbs_path(
    selected.tau, selected.right, L, M, dm, lr_m, split_evaluator,
    prefer_last_tie, tau_hat_all, gain_all
  );
}

}

arma::mat sorted_cpd_path(const arma::vec & tau_hat) {
  const int n_tau_hat = tau_hat.n_rows;
  arma::mat tau_cand(n_tau_hat, n_tau_hat);
  tau_cand.fill(arma::datum::nan);
  for (int i = 0; i < n_tau_hat; i++) {
    arma::vec tau_hat_i = tau_hat(arma::span(0, i));
    tau_cand(i, arma::span(0, i)) = arma::sort(tau_hat_i.t());
  }
  return tau_cand;
}

SingleCpdResult sn_one_loss_output(const int & n,
                                   const int & L,
                                   const int & dm,
                                   CostEngine & cost_engine) {
  SingleCpdResult out;
  auto tstart = std::chrono::steady_clock::now();
  if (L == 0) {
    arma::vec loss_vec(1);
    arma::ivec cps_num(1);
    loss_vec(0) = cost_engine.get_cost(1, n);
    cps_num(0) = 0;
    auto tend = std::chrono::steady_clock::now();
    set_cpd_candidates_from_matrix(out, arma::mat(0, 0), true);
    set_cost_engine_metrics(
      out, cost_engine, std::chrono::duration<double>(tend - tstart).count()
    );
    out.loss = loss_vec;
    out.objective = loss_vec;
    out.cps_num = cps_num;
    return out;
  }

  int Q = L + 1;
  double total_time = 0;

  arma::mat loss(Q + 1, n + 1); loss.fill(arma::datum::nan);
  arma::vec loss_vec(Q); loss_vec.fill(arma::datum::nan);
  arma::ivec cps_num(Q);
  for (int i = 0; i < Q; ++i) {
    cps_num(i) = i;
  }
  arma::uvec v_min_vec;
  arma::uword v_min;

  for(int j = dm; j <= n; j++) {
    loss(1, j) = cost_engine.get_cost(1, j);
  }
  loss_vec(0) = loss(1, n);
  arma::mat V(Q + 1, n + 1); V.fill(arma::datum::nan);
  for (int q = 2; q <= Q - 1; q++) {
    for (int j = q * dm; j <= n; j++) {
      arma::vec loss_temp(j - q * dm + 1); loss_temp.fill(arma::datum::nan);
      for (int v = (j - dm); v >= (q - 1) * dm; v--) {
        loss_temp(v - (q - 1) * dm) =
          loss(q - 1, v) + cost_engine.get_cost(v + 1, j);
      }
      v_min = loss_temp.index_min();
      v_min_vec = arma::find(arma::abs(loss_temp - loss_temp(v_min)) < 1e-8, 1, "last");
      v_min = v_min_vec(0);
      loss(q, j) = loss_temp(v_min);
      V(q, j) = v_min + (q - 1) * dm;
    }
    loss_vec(q - 1) = loss(q, n);
  }

  arma::vec loss_temp(n - Q * dm + 1); loss_temp.fill(arma::datum::nan);
  for (int v = (n - dm); v >= (Q - 1) * dm; v--) {
    loss_temp(v - (Q - 1) * dm) =
      loss(Q - 1, v) + cost_engine.get_cost(v + 1, n);
  }
  v_min = loss_temp.index_min();
  v_min_vec = arma::find(arma::abs(loss_temp - loss_temp(v_min)) < 1e-8, 1, "last");
  v_min = v_min_vec(0);
  loss(Q, n) = loss_temp(v_min);
  V(Q, n) = v_min + (Q - 1) * dm;
  loss_vec(Q - 1) = loss(Q, n);

  arma::mat cps(Q + 1, L + 1); cps.fill(arma::datum::nan);
  cps.col(1) = V.col(n);
  if (Q >= 3) {
    for (int q = 3; q <= Q; q++) {
      for (int i = 2; i <= q - 1; i++) {
        cps(q, i) = V(q - (i - 1), cps(q, i - 1));
      }
    }
  }
  arma::mat cpd_cand = cps(arma::span(2, Q), arma::span(1, L));
  arma::mat temp = cpd_cand;
  int n_cpd_cand = cpd_cand.n_rows;
  for (int i = 0; i < n_cpd_cand; i++) {
    for (int j = 0; j <= i; j++) {
      cpd_cand(i, j) = temp(i, i - j);
    }
  }
  auto tend = std::chrono::steady_clock::now();
  total_time = std::chrono::duration<double>(tend - tstart).count();
  set_cpd_candidates_from_matrix(out, cpd_cand, true);
  set_cost_engine_metrics(out, cost_engine, total_time);
  out.loss = loss_vec;
  out.objective = loss_vec;
  out.cps_num = cps_num;
  return out;
}

WbsSearchPath wbs_search_path(const int n,
                              const int L,
                              const int dm,
                              const arma::imat & lr_m,
                              SplitEvaluator & split_evaluator,
                              const bool & prefer_last_tie) {
  WbsSearchPath out;
  int M = lr_m.n_rows;
  auto tstart = std::chrono::steady_clock::now();

  std::vector<SegmentSplit> active_segments;
  SegmentSplit initial = best_split_for_segment(
    0, n, 0, M, dm, lr_m, split_evaluator, prefer_last_tie
  );
  if (initial.finite()) {
    active_segments.push_back(initial);
  }

  std::vector<int> tau_hat_all;
  std::vector<double> gain_all;
  tau_hat_all.reserve(L);
  gain_all.reserve(L);

  for (int k = 1; k <= L; k++) {
    if (active_segments.empty()) {
      break;
    }

    std::sort(active_segments.begin(), active_segments.end(), segment_order);

    int best_segment = -1;
    for (int j = 0; j < static_cast<int>(active_segments.size()); ++j) {
      const SegmentSplit & candidate = active_segments[j];
      if (!candidate.finite()) {
        continue;
      }
      if (best_segment < 0) {
        best_segment = j;
        continue;
      }

      const SegmentSplit & current = active_segments[best_segment];
      const bool better_gain = candidate.gain > current.gain;
      const bool tied_gain = std::abs(candidate.gain - current.gain) < 1e-8;
      const bool better_tie = prefer_last_tie && tied_gain &&
        (
          candidate.interval_col > current.interval_col ||
          (
            candidate.interval_col == current.interval_col &&
            segment_order(current, candidate)
          )
        );

      if (better_gain || better_tie) {
        best_segment = j;
      }
    }

    if (best_segment < 0) {
      break;
    }

    const SegmentSplit selected = active_segments[best_segment];
    tau_hat_all.push_back(selected.tau);
    gain_all.push_back(selected.gain);
    active_segments.erase(active_segments.begin() + best_segment);

    if (k == L) {
      continue;
    }

    SegmentSplit left_segment = best_split_for_segment(
      selected.left, selected.tau, 0, M, dm, lr_m, split_evaluator,
      prefer_last_tie
    );
    if (left_segment.finite()) {
      active_segments.push_back(left_segment);
    }

    SegmentSplit right_segment = best_split_for_segment(
      selected.tau, selected.right, 0, M, dm, lr_m, split_evaluator,
      prefer_last_tie
    );
    if (right_segment.finite()) {
      active_segments.push_back(right_segment);
    }
  }

  arma::vec tau_hat = finite_inner_changepoints(selected_tau_vector(tau_hat_all), n);
  auto tend = std::chrono::steady_clock::now();
  out.total_time = std::chrono::duration<double>(tend - tstart).count();
  out.tau_hat = tau_hat;
  out.gain = selected_gain_vector(gain_all);
  return out;
}

WbsSearchPath wbs_recursive_search_path(const int n,
                                        const int L,
                                        const int dm,
                                        const arma::imat & lr_m,
                                        SplitEvaluator & split_evaluator,
                                        const bool & prefer_last_tie) {
  WbsSearchPath out;
  int M = lr_m.n_rows;
  auto tstart = std::chrono::steady_clock::now();

  std::vector<int> tau_hat_all;
  std::vector<double> gain_all;
  tau_hat_all.reserve(L);
  gain_all.reserve(L);

  append_recursive_wbs_path(
    0, n, L, M, dm, lr_m, split_evaluator, prefer_last_tie,
    tau_hat_all, gain_all
  );

  arma::vec tau_hat = finite_inner_changepoints(
    selected_tau_vector(tau_hat_all), n
  );
  auto tend = std::chrono::steady_clock::now();
  out.total_time = std::chrono::duration<double>(tend - tstart).count();
  out.tau_hat = tau_hat;
  out.gain = selected_gain_vector(gain_all);
  return out;
}

SingleCpdResult wbs_one_loss_output(const int n,
                                    const int L,
                                    const int dm,
                                    const arma::imat & lr_m,
                                    CostEngine & cost_engine,
                                    const bool & recursive) {
  SingleCpdResult out;
  ExactSplitEvaluator split_evaluator(cost_engine);
  WbsSearchPath path = recursive ?
    wbs_recursive_search_path(n, L, dm, lr_m, split_evaluator) :
    wbs_search_path(n, L, dm, lr_m, split_evaluator);
  arma::mat tau_cand = sorted_cpd_path(path.tau_hat);
  const int n_tau_hat = path.tau_hat.n_rows;

  arma::vec loss_vec(n_tau_hat + 1); loss_vec.fill(0);
  arma::ivec cps_num(n_tau_hat + 1);
  for (int i = 0; i <= n_tau_hat; ++i) {
    cps_num(i) = i;
  }

  loss_vec(0) = cost_engine.get_cost(1, n);
  for(int i = 1; i <= n_tau_hat; i++) {
    for(int j = 0; j <= i; j++) {
      if(j == 0) {
        loss_vec(i) += cost_engine.get_cost(1, tau_cand(i - 1, j));
      } else if(j < i) {
        loss_vec(i) += cost_engine.get_cost(
          tau_cand(i - 1, j - 1) + 1, tau_cand(i - 1, j)
        );
      } else {
        loss_vec(i) += cost_engine.get_cost(
          tau_cand(i - 1, j - 1) + 1, n
        );
      }
    }
  }

  set_cpd_candidates_from_matrix(out, tau_cand, true);
  set_cost_engine_metrics(out, cost_engine, path.total_time);
  out.loss = loss_vec;
  out.objective = loss_vec;
  out.cps_num = cps_num;
  out.path_score = path.gain;
  return out;
}

SingleCpdResult pelt_one_loss_output(const int & n,
                                     const arma::vec & pen_val,
                                     const int & dm,
                                     const double & prune_value,
                                     CostEngine & cost_engine) {
  SingleCpdResult out;
  double total_time = 0;

  int n_pen_val = pen_val.n_elem;
  auto tstart = std::chrono::steady_clock::now();

  arma::vec loss_vec(n_pen_val + 1); loss_vec.fill(arma::datum::nan);
  arma::vec objective_vec(n_pen_val + 1); objective_vec.fill(arma::datum::nan);
  arma::ivec cps_num(n_pen_val + 1); cps_num.fill(0);
  arma::mat num_pruned(n_pen_val, n + 1, arma::fill::zeros);
  const double full_loss = cost_engine.get_cost(1, n);
  loss_vec(0) = full_loss;
  objective_vec(0) = full_loss;

  if (n < 2 * dm) {
    loss_vec.fill(full_loss);
    objective_vec.fill(full_loss);
    arma::mat empty_cpd(n_pen_val, 0);
    auto tend = std::chrono::steady_clock::now();
    set_cpd_candidates_from_matrix(out, empty_cpd, true);
    set_cost_engine_metrics(
      out, cost_engine, std::chrono::duration<double>(tend - tstart).count()
    );
    out.loss = loss_vec;
    out.objective = objective_vec;
    out.cps_num = cps_num;
    out.num_pruned = num_pruned;
    return out;
  }

  arma::mat F_cost(n + 1, n_pen_val); F_cost.fill(arma::datum::nan);
  F_cost.row(0) = -pen_val.t();
  for (int ts = dm; ts < 2 * dm; ts++) {
    for (int b = 0; b < n_pen_val; b++) {
      F_cost(ts, b) = cost_engine.get_cost(1, ts);
    }
  }
  arma::mat cp(n + 1, n_pen_val); cp.fill(0);
  int ncp_max = n / dm;
  int count, cp0;
  arma::mat cps_final(n_pen_val, ncp_max); cps_final.fill(arma::datum::nan);
  arma::vec R_ts(n + 1);
  int n_R;
  arma::vec temp;
  for (int b = 0; b < n_pen_val; b++) {
    arma::ivec prune_after(n + 1);
    prune_after.fill(n + dm + 1);
    R_ts(0) = 0;
    R_ts(1) = dm;
    n_R = 2;
    for (int ts = 2 * dm; ts <= n; ts++) {
      int n_retained = 0;
      int n_removed = 0;
      for (int i = 0; i < n_R; i++) {
        const int candidate = static_cast<int>(R_ts(i));
        if (prune_after(candidate) > ts) {
          R_ts(n_retained++) = candidate;
        } else {
          n_removed += 1;
        }
      }
      n_R = n_retained;
      num_pruned(b, ts) = n_removed;

      temp.set_size(n_R);
      temp.fill(arma::datum::nan);
      for (int i = 0; i < n_R; i++) {
        temp(i) = F_cost(R_ts(i), b) + cost_engine.get_cost(R_ts(i) + 1, ts) + pen_val(b);
      }
      arma::uword idx = temp.index_min();
      F_cost(ts, b) = temp(idx);
      cp(ts, b) = R_ts(idx);
      if (!(std::isinf(prune_value) && prune_value < 0)) {
        for (int i = 0; i < n_R; i++) {
          if (!(temp(i) + prune_value <
                F_cost(ts, b) + pen_val(b))) {
            const int candidate = static_cast<int>(R_ts(i));
            prune_after(candidate) = std::min(
              prune_after(candidate), ts + dm
            );
          }
        }
      }
      if (ts < n) {
        R_ts(n_R) = ts - (dm - 1);
        n_R += 1;
      }
    }
    objective_vec(b + 1) = F_cost(n, b);
    cp0 = cp(n, b);
    count = 0;
    while (cp0 > 0) {
      cps_final(b, count) = cp0;
      cp0 = cp(cp0, b);
      count += 1;
    }
    if(count > 0) {
      arma::rowvec tau_hat_i = cps_final(b, arma::span(0, count - 1));
      cps_final(b, arma::span(0, count - 1)) = arma::sort(tau_hat_i);
    }
    loss_vec(b + 1) = objective_vec(b + 1) - count * pen_val(b);
    cps_num(b + 1) = count;
  }
  auto tend = std::chrono::steady_clock::now();
  total_time = std::chrono::duration<double>(tend - tstart).count();

  set_cpd_candidates_from_matrix(out, cps_final, true);
  set_cost_engine_metrics(out, cost_engine, total_time);
  out.loss = loss_vec;
  out.objective = objective_vec;
  out.cps_num = cps_num;
  out.num_pruned = num_pruned;
  return out;
}

SingleCpdResult cpd_one_loss_output(
  const CpdSearchMethod method,
  const int & n,
  const int & L,
  const int & dm,
  const arma::imat & search_intervals,
  const arma::vec & pen_val,
  const double & prune_value,
  CostEngine & cost_engine
) {
  switch (method) {
  case CpdSearchMethod::SegmentNeighbourhood:
    return sn_one_loss_output(n, L, dm, cost_engine);
  case CpdSearchMethod::Wbs:
    return wbs_one_loss_output(
      n, L, dm, search_intervals, cost_engine, false
    );
  case CpdSearchMethod::WbsRecursive:
    return wbs_one_loss_output(
      n, L, dm, search_intervals, cost_engine, true
    );
  case CpdSearchMethod::Pelt:
    return pelt_one_loss_output(
      n, pen_val, dm, prune_value, cost_engine
    );
  }
  throw std::invalid_argument("unsupported changepoint search method.");
}
