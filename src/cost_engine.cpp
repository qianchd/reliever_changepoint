#include "cost_engine.h"
#include "mean_utils.h"
#include "relief_interval.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

long long cost_interval_index(const int & i, const int & j) {
  return static_cast<long long>(i) +
    static_cast<long long>(j) * (j - 1) / 2;
}

int exact2relief_itv_routine_id_only(
  const int & l,
  const int & r,
  const arma::ivec & miss_cover_len,
  const arma::ivec & int_len,
  const arma::ivec & layer_point,
  const arma::imat & int_eps
) {
  const int len = r - l;
  const int k_begin = static_cast<int>(
    std::upper_bound(
      miss_cover_len.begin(), miss_cover_len.end(), len - 1
    ) - miss_cover_len.begin()
  ) - 1;
  if (k_begin < 0) {
    return 0;
  }

  const int k_end = static_cast<int>(
    std::upper_bound(int_len.begin(), int_len.end(), len) - int_len.begin()
  ) - 1;
  if (k_end < k_begin) {
    return 0;
  }

  for (int k = k_end; k >= k_begin; --k) {
    const int id_first = k == 0 ? 0 : layer_point(k - 1);
    const int id_last = layer_point(k) - 1;
    const int max_left = r - int_len(k);

    int lo = id_first;
    int hi = id_last;
    int best = -1;
    while (lo <= hi) {
      const int mid = lo + (hi - lo) / 2;
      if (int_eps(mid, 0) <= max_left) {
        best = mid;
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }

    if (best >= id_first && int_eps(best, 0) >= l &&
        int_eps(best, 1) <= r) {
      return best + 1;
    }
  }
  return 0;
}

template <typename Callback>
void for_each_sawtooth_interval(
  const arma::imat & sawtooth_points,
  const int & l,
  const int & r,
  Callback callback
) {
  for (unsigned int s = 0; s < sawtooth_points.n_rows; ++s) {
    const int right_begin = s == 0 ? r : sawtooth_points(s - 1, 1) + 1;
    for (int ll = l; ll >= sawtooth_points(s, 0); --ll) {
      for (int rr = right_begin; rr <= sawtooth_points(s, 1); ++rr) {
        callback(ll, rr);
      }
    }
  }
}

}

arma::mat RegLossFunction::individual_loss(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  (void) l;
  (void) r;
  (void) l_end;
  (void) r_end;
  throw std::runtime_error(
    "This RegLossFunction does not provide individual_loss."
  );
}

arma::mat RegLossFunction::block_loss(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  const int left_n = l - l_end;
  const int core_n = r - l + 1;
  const int right_n = r_end - r;
  if (left_n < 0 || core_n < 1 || right_n < 0) {
    throw std::runtime_error("invalid loss-block interval bounds.");
  }

  const int n_block_cols = 1 + left_n + right_n;
  arma::mat individual = individual_loss(l, r, l_end, r_end);
  if (static_cast<int>(individual.n_rows) != left_n + core_n + right_n) {
    throw std::runtime_error(
      "individual_loss returned an incompatible evaluation length."
    );
  }

  arma::mat loss(individual.n_cols, n_block_cols, arma::fill::zeros);
  loss.col(0) = arma::sum(
    individual.rows(left_n, left_n + core_n - 1),
    0
  ).t();
  if (left_n > 0) {
    loss.cols(1, left_n) =
      arma::cumsum(arma::flipud(individual.rows(0, left_n - 1)), 0).t();
  }
  if (right_n > 0) {
    loss.cols(1 + left_n, left_n + right_n) =
      arma::cumsum(
        individual.rows(left_n + core_n, left_n + core_n + right_n - 1),
        0
      ).t();
  }
  return loss;
}

double RegLossFunction::model_fit_time() const {
  return 0.0;
}

MeanSquareRegLossFunction::MeanSquareRegLossFunction(
  const double & ratio,
  const arma::mat & data
) : data(data),
    sq_smry(arma::mean(arma::square(data), 1)),
    ratio(ratio) {
  mu_hat = arma::vec(this->data.n_cols, arma::fill::zeros);
}

void MeanSquareRegLossFunction::fit_core_interval(
  const unsigned int & l,
  const unsigned int & r
) {
  if (this->len_current > 0 &&
      this->l_current == l &&
      this->r_current == r) {
    return;
  }
  this->live_update(l, r);
  this->stabilize_core_mean();
}

void MeanSquareRegLossFunction::live_update(
  const unsigned int & l,
  const unsigned int & r
) {
  unsigned int len_new = r - l + 1;
  if (this->l_current == 0) {
    this->l_current = l;
    this->r_current = l - 1;
  }

  arma::uvec set_add;
  arma::uvec set_minus;
  if (l < l_current) {
    set_add = arma::linspace<arma::uvec>(l - 1, l_current - 2, l_current - l);
  } else if (l > l_current) {
    set_minus = arma::linspace<arma::uvec>(l_current - 1, l - 2, l - l_current);
  }
  if (r > r_current) {
    set_add = arma::join_vert(
      set_add,
      arma::linspace<arma::uvec>(r_current, r - 1, r - r_current)
    );
  } else if (r < r_current) {
    set_minus = arma::join_vert(
      set_minus,
      arma::linspace<arma::uvec>(r, r_current - 1, r_current - r)
    );
  }

  if (static_cast<double>(set_add.n_elem + set_minus.n_elem) >
      static_cast<double>(len_new) * ratio) {
    this->mu_hat = arma::mean(data.rows(l - 1, r - 1), 0).t();
    this->squ = arma::accu(sq_smry.subvec(l - 1, r - 1));
  } else {
    const double sum_add = set_add.n_elem > 0 ? arma::accu(sq_smry(set_add)) : 0.0;
    const double sum_minus = set_minus.n_elem > 0 ? arma::accu(sq_smry(set_minus)) : 0.0;
    this->squ += sum_add - sum_minus;

    arma::vec col_sum_add;
    arma::vec col_sum_minus;
    if (set_add.n_elem > 0) {
      col_sum_add = arma::sum(data.rows(set_add), 0).t();
    } else {
      col_sum_add = arma::zeros<arma::vec>(data.n_cols);
    }
    if (set_minus.n_elem > 0) {
      col_sum_minus = arma::sum(data.rows(set_minus), 0).t();
    } else {
      col_sum_minus = arma::zeros<arma::vec>(data.n_cols);
    }
    this->mu_hat =
      (this->mu_hat * len_current + col_sum_add - col_sum_minus) / len_new;
  }

  this->l_current = l;
  this->r_current = r;
  this->len_current = len_new;
}

void MeanSquareRegLossFunction::refresh_core_loss_assessment() {
  const double fitted_square =
    arma::dot(this->mu_hat, this->mu_hat) /
    static_cast<double>(this->mu_hat.n_elem) *
    this->len_current;
  this->core_fast_loss = this->squ - fitted_square;
  const double cancellation_scale =
    std::abs(this->squ) + std::abs(fitted_square);
  const double error_bound =
    16.0 * std::numeric_limits<double>::epsilon() *
    std::max(1.0, cancellation_scale) *
    std::max(1.0, static_cast<double>(this->len_current));
  this->core_loss_risky =
    !std::isfinite(this->core_fast_loss) ||
    this->core_fast_loss <= error_bound;
}

void MeanSquareRegLossFunction::stabilize_core_mean() {
  this->refresh_core_loss_assessment();
  if (this->core_loss_risky) {
    this->mu_hat = this->stable_core_mean();
    this->refresh_core_loss_assessment();
  }
}

arma::vec MeanSquareRegLossFunction::stable_core_mean() const {
  return stable_col_means_double(
    data.rows(this->l_current - 1, this->r_current - 1)
  );
}

arma::vec MeanSquareRegLossFunction::squared_residual_loss(
  const arma::uvec & row_idx
) const {
  arma::mat residual = data.rows(row_idx);
  residual.each_row() -= this->mu_hat.t();
  return arma::mean(arma::square(residual), 1);
}

double MeanSquareRegLossFunction::stable_core_loss() const {
  if (!this->core_loss_risky) {
    return this->core_fast_loss;
  }

  arma::uvec core_idx = arma::regspace<arma::uvec>(
    this->l_current - 1, this->r_current - 1
  );
  return arma::accu(this->squared_residual_loss(core_idx));
}

arma::mat MeanSquareRegLossFunction::individual_loss(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  if (l_end > r_end) {
    return arma::mat(0, 1);
  }

  const auto tstart = std::chrono::steady_clock::now();
  this->fit_core_interval(l, r);
  arma::uvec eval_idx =
    arma::linspace<arma::uvec>(l_end - 1, r_end - 1, r_end - l_end + 1);
  arma::vec loss_obs = this->squared_residual_loss(eval_idx);
  const auto tend = std::chrono::steady_clock::now();
  this->model_fit_time_ +=
    std::chrono::duration<double>(tend - tstart).count();
  return loss_obs;
}

arma::mat MeanSquareRegLossFunction::block_loss(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  const auto tstart = std::chrono::steady_clock::now();
  this->fit_core_interval(l, r);

  arma::mat loss(1, l - l_end + r_end - r + 1, arma::fill::zeros);
  loss(0) = this->stable_core_loss();
  if (l == l_end && r == r_end) {
    const auto tend = std::chrono::steady_clock::now();
    this->model_fit_time_ +=
      std::chrono::duration<double>(tend - tstart).count();
    return loss;
  }

  if (l_end < l && r_end > r) {
    arma::uvec l_idx = arma::linspace<arma::uvec>(l - 2, l_end - 1, l - l_end);
    arma::uvec r_idx = arma::linspace<arma::uvec>(r, r_end - 1, r_end - r);
    arma::uvec tem_idx = arma::join_vert(l_idx, r_idx);
    arma::vec loss_obs = this->squared_residual_loss(tem_idx);
    loss.cols(1, l - l_end) =
      arma::cumsum(loss_obs.head(l - l_end)).t();
    loss.tail_cols(r_end - r) =
      arma::cumsum(loss_obs.tail(r_end - r)).t();
  } else if (l_end < l && r_end == r) {
    arma::uvec tem_idx =
      arma::linspace<arma::uvec>(l - 2, l_end - 1, l - l_end);
    arma::vec loss_obs = this->squared_residual_loss(tem_idx);
    loss.tail_cols(l - l_end) = arma::cumsum(loss_obs).t();
  } else {
    arma::uvec tem_idx = arma::linspace<arma::uvec>(r, r_end - 1, r_end - r);
    arma::vec loss_obs = this->squared_residual_loss(tem_idx);
    loss.tail_cols(r_end - r) = arma::cumsum(loss_obs).t();
  }
  const auto tend = std::chrono::steady_clock::now();
  this->model_fit_time_ +=
    std::chrono::duration<double>(tend - tstart).count();
  return loss;
}

double MeanSquareRegLossFunction::model_fit_time() const {
  return this->model_fit_time_;
}

CostEngine::CostEngine(const int & cost_row_) : cost_row(cost_row_) {}

int CostEngine::get_n_model_fit() {
  return this->n_model_fit;
}

double CostEngine::get_model_fit_time() {
  return this->model_fit_time;
}

CostMatCache::CostMatCache(
  arma::mat & cost_mat,
  const arma::ivec & miss_cover_len,
  const arma::ivec & int_len,
  const arma::ivec & layer_point,
  const arma::imat & int_eps,
  const int & n
) : cost_mat(cost_mat),
    miss_cover_len(miss_cover_len),
    int_len(int_len),
    layer_point(layer_point),
    int_eps(int_eps),
    n(n) {}

double CostMatCache::cached_cost(
  const int & row,
  const int & col
) {
  return this->cost_mat(row, col);
}

void CostMatCache::write_cached_cost(
  const int & col,
  const arma::vec & value
) {
  this->cost_mat.col(col) = value;
}

arma::ivec CostMatCache::find_relief_interval(
  const int & l,
  const int & r
) {
  return exact2relief_itv_routine_c(
    l, r, this->miss_cover_len, this->int_len, this->layer_point, this->int_eps
  );
}

arma::imat CostMatCache::covered_intervals(const int & id) {
  return relief2exact_itv_routine_c(
    id, this->int_len, this->layer_point, this->int_eps, this->n
  );
}

CostEngineByCostMat::CostEngineByCostMat(
  CostMatCache & cost_mat_cache,
  RegLossFunction & reg_loss,
  const int & cost_row,
  const bool & is_full_search
) : CostEngine(cost_row),
    cost_mat_cache(cost_mat_cache),
    reg_loss(reg_loss),
    is_full_search(is_full_search) {}

double CostEngineByCostMat::get_cost(const int & i, const int & j) {
  const long long interval_index = cost_interval_index(i, j);
  if (std::isinf(this->cost_mat_cache.cached_cost(cost_row - 1, interval_index))) {
    auto tstart = std::chrono::steady_clock::now();
    update_cache(i, j);
    n_model_fit += 1;
    auto tend = std::chrono::steady_clock::now();
    model_fit_time += std::chrono::duration<double>(tend - tstart).count();
  }
  return this->cost_mat_cache.cached_cost(cost_row - 1, interval_index);
}

arma::mat CostEngineByCostMat::fit_loss_block(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  return reg_loss.block_loss(l, r, l_end, r_end);
}

void CostEngineByCostMat::update_exact_cache(
  const int & i,
  const int & j
) {
  arma::mat loss = fit_loss_block(i, j, i, j);
  this->cost_mat_cache.write_cached_cost(cost_interval_index(i, j), loss.col(0));
}

void CostEngineByCostMat::update_cache(const int & i, const int & j) {
  if (is_full_search) {
    update_exact_cache(i, j);
    return;
  }

  arma::ivec relief_int = this->cost_mat_cache.find_relief_interval(i - 1, j);

  const int id = relief_int(0);
  const int l = relief_int(1) + 1;
  const int r = relief_int(2);

  arma::imat sawtooth_points = this->cost_mat_cache.covered_intervals(id);

  const int l_end = arma::min(sawtooth_points.col(0));
  const int r_end = arma::max(sawtooth_points.col(1));
  arma::mat loss = fit_loss_block(l, r, l_end, r_end);

  for_each_sawtooth_interval(
    sawtooth_points,
    l,
    r,
    [&](const int & ll, const int & rr) {
      arma::vec value = loss.col(0);
      if (ll < l) {
        value += loss.col(l - ll);
      }
      if (rr > r) {
        value += loss.col(rr - r + l - l_end);
      }
      this->cost_mat_cache.write_cached_cost(
        cost_interval_index(ll, rr), value
      );
    }
  );
}

LossBlockCache::LossBlockCache(
  RegLossFunction & reg_loss,
  const arma::ivec & miss_cover_len,
  const arma::ivec & int_len,
  const arma::ivec & layer_point,
  const arma::imat & int_eps,
  const int & n,
  const bool & use_owner_key
) : reg_loss(reg_loss),
    miss_cover_len(miss_cover_len),
    int_len(int_len),
    layer_point(layer_point),
    int_eps(int_eps),
    n(n),
    use_owner_key(use_owner_key) {
  relief_blocks.resize(static_cast<std::size_t>(int_eps.n_rows) + 1U);
  relief_block_ready.assign(static_cast<std::size_t>(int_eps.n_rows) + 1U, 0);
  if (use_owner_key) {
    owner_key_ready.assign(static_cast<std::size_t>(int_eps.n_rows) + 1U, 0);
    const std::size_t n_owner_key =
      static_cast<std::size_t>(n) * (static_cast<std::size_t>(n) + 1U) / 2U + 1U;
    owner_key.assign(n_owner_key, 0);
  }
}

double LossBlockCache::get_cost(
  const int & cost_row,
  const int & i,
  const int & j
) {
  get_cost_calls += 1.0;

  if (i < 1 || j < i || j > n) {
    throw std::runtime_error(
      "interval bounds are out of range for the by_loss_block cache."
    );
  }

  const long long interval_index = cost_interval_index(i, j);

  int key = 0;
  if (use_owner_key) {
    key = owner_key[static_cast<std::size_t>(interval_index)];
  }
  if (!use_owner_key || key == 0) {
    key = locate_block_key(i, j, interval_index);
  }
  const LossBlockRecord & block = cached_block(key);
  const int row = cost_row - 1;
  if (row < 0 || row >= static_cast<int>(block.loss.n_rows)) {
    throw std::runtime_error(
      "loss output index is out of range for the cached relief block."
    );
  }

  double out = block.loss(row, 0);
  if (i < block.l) {
    out += block.loss(row, block.l - i);
  }
  if (j > block.r) {
    out += block.loss(row, j - block.r + block.l - block.l_end);
  }
  return out;
}

void LossBlockCache::restore(const LossBlockCacheState & cache_state) {
  if (cache_state.n != n) {
    throw std::runtime_error(
      "loss_block_cache was built for a different n."
    );
  }

  relief_blocks.assign(static_cast<std::size_t>(int_eps.n_rows) + 1U,
                       LossBlockRecord());
  relief_block_ready.assign(static_cast<std::size_t>(int_eps.n_rows) + 1U, 0);
  if (use_owner_key) {
    owner_key_ready.assign(static_cast<std::size_t>(int_eps.n_rows) + 1U, 0);
  }
  exact_block_cache.clear();
  if (use_owner_key) {
    std::fill(owner_key.begin(), owner_key.end(), 0);
  }

  for (const LossBlockStateRecord & record : cache_state.blocks) {
    const int key = record.key;
    const LossBlockRecord & block = record.block;
    if (key > 0) {
      if (key >= static_cast<int>(relief_blocks.size())) {
        throw std::runtime_error(
          "loss_block_cache contains a relief block id outside int_eps."
        );
      }
      relief_blocks[key] = block;
      relief_block_ready[key] = 1;
    } else {
      exact_block_cache[key] = block;
    }
  }

  if (use_owner_key && cache_state.has_owner_key) {
    if (cache_state.owner_key.size() != owner_key.size()) {
      throw std::runtime_error(
        "loss_block_cache owner_key has incompatible length."
      );
    }
    owner_key = cache_state.owner_key;
    for (std::size_t id = 1; id < relief_block_ready.size(); ++id) {
      if (relief_block_ready[id]) {
        owner_key_ready[id] = 1;
      }
    }
  }

  raw_loss_block_cells = cache_state.raw_loss_block_cells;
  get_cost_calls = cache_state.get_cost_calls;
  expanded_cost_cells = cache_state.expanded_cost_cells;
  expanded_interval_writes = cache_state.expanded_interval_writes;
  full_update_calls = cache_state.full_update_calls;
  relief_update_calls = cache_state.relief_update_calls;
  model_fit_time = cache_state.model_fit_time;
  exact_fallback_calls = 0.0;
  exact_fallback_i.clear();
  exact_fallback_j.clear();
}

LossBlockCacheState LossBlockCache::state() const {
  std::size_t block_count = exact_block_cache.size();
  for (std::size_t id = 1; id < relief_block_ready.size(); ++id) {
    if (relief_block_ready[id]) {
      block_count += 1U;
    }
  }

  LossBlockCacheState out;
  out.n = n;
  out.blocks.reserve(block_count);
  for (std::size_t id = 1; id < relief_blocks.size(); ++id) {
    if (!relief_block_ready[id]) {
      continue;
    }
    LossBlockStateRecord record;
    record.key = static_cast<int>(id);
    record.block = relief_blocks[id];
    out.blocks.push_back(record);
  }
  for (const auto & entry : exact_block_cache) {
    LossBlockStateRecord record;
    record.key = entry.first;
    record.block = entry.second;
    out.blocks.push_back(record);
  }

  out.raw_loss_block_cells = raw_loss_block_cells;
  out.get_cost_calls = get_cost_calls;
  out.expanded_cost_cells = expanded_cost_cells;
  out.expanded_interval_writes = expanded_interval_writes;
  out.full_update_calls = full_update_calls;
  out.relief_update_calls = relief_update_calls;
  out.model_fit_time = model_fit_time;
  if (use_owner_key) {
    out.has_owner_key = true;
    out.owner_key = owner_key;
  }
  return out;
}

std::string LossBlockCache::exact_fallback_warning() const {
  if (exact_fallback_calls <= 0.0) {
    return "";
  }

  std::ostringstream message;
  message
    << "Reliever by_loss_block computed exact losses for "
    << static_cast<long long>(exact_fallback_calls)
    << " interval(s) not covered by any relief block. "
    << "This suggests create_relief_itv/exact2relief_itv_routine_c/"
    << "relief2exact_itv_routine_c coverage is inconsistent. "
    << "Results remain exact because these intervals were fitted directly; "
    << "please report this warning with a reproducible example.";
  if (!exact_fallback_i.empty()) {
    message << " First intervals:";
    for (std::size_t k = 0; k < exact_fallback_i.size(); ++k) {
      message << (k == 0 ? " " : ", ")
              << "(" << exact_fallback_i[k] << "," << exact_fallback_j[k] << ")";
    }
    if (exact_fallback_calls > static_cast<double>(exact_fallback_i.size())) {
      message << ", ...";
    }
  }
  return message.str();
}

int LossBlockCache::locate_block_key(
  const int & i,
  const int & j,
  const long long & interval_index
) {
  const int id = exact2relief_itv_routine_id_only(
    i - 1,
    j,
    miss_cover_len,
    int_len,
    layer_point,
    int_eps
  );
  if (id > 0) {
    ensure_relief_block(id);
    if (use_owner_key) {
      owner_key[static_cast<std::size_t>(interval_index)] = id;
    }
    return id;
  }

  const int key = -static_cast<int>(interval_index);
  const bool existed = exact_block_cache.find(key) != exact_block_cache.end();
  ensure_exact_block(key, i, j);
  if (use_owner_key) {
    owner_key[static_cast<std::size_t>(interval_index)] = key;
  }
  if (!existed) {
    expanded_interval_writes += 1.0;
    expanded_cost_cells += exact_block_cache.at(key).loss.n_rows;
    exact_fallback_calls += 1.0;
    if (exact_fallback_i.size() < 20U) {
      exact_fallback_i.push_back(i);
      exact_fallback_j.push_back(j);
    }
  }
  return key;
}

const LossBlockRecord & LossBlockCache::cached_block(const int & key) const {
  if (key > 0) {
    if (key >= static_cast<int>(relief_blocks.size()) ||
        !relief_block_ready[key]) {
      throw std::runtime_error(
        "relief block is missing from the by_loss_block cache."
      );
    }
    return relief_blocks[key];
  }
  return exact_block_cache.at(key);
}

void LossBlockCache::ensure_relief_block(const int & id) {
  if (id <= 0 || id >= static_cast<int>(relief_blocks.size())) {
    throw std::runtime_error("relief block id is outside int_eps.");
  }
  if (relief_block_ready[id] &&
      (!use_owner_key || owner_key_ready[id])) {
    return;
  }

  const int l = int_eps(id - 1, 0) + 1;
  const int r = int_eps(id - 1, 1);
  arma::imat sawtooth_points = relief2exact_itv_routine_c(
    id,
    int_len,
    layer_point,
    int_eps,
    n
  );
  const int l_end = arma::min(sawtooth_points.col(0));
  const int r_end = arma::max(sawtooth_points.col(1));

  if (relief_block_ready[id]) {
    fill_owner_keys(id, sawtooth_points, l, r);
    owner_key_ready[id] = 1;
    return;
  }

  LossBlockRecord & block = relief_blocks[id];
  block.l = l;
  block.r = r;
  block.l_end = l_end;
  block.r_end = r_end;
  block.loss = fit_loss_block(l, r, l_end, r_end);
  const arma::uword loss_n_rows = block.loss.n_rows;
  const double interval_writes = count_covered_intervals(sawtooth_points, l, r);
  if (use_owner_key) {
    fill_owner_keys(id, sawtooth_points, l, r);
  }
  relief_block_ready[id] = 1;
  if (use_owner_key) {
    owner_key_ready[id] = 1;
  }
  expanded_interval_writes += interval_writes;
  expanded_cost_cells += interval_writes * loss_n_rows;
  relief_update_calls += 1.0;
}

void LossBlockCache::fill_owner_keys(
  const int & id,
  const arma::imat & sawtooth_points,
  const int & l,
  const int & r
) {
  for_each_sawtooth_interval(
    sawtooth_points,
    l,
    r,
    [&](const int & ll, const int & rr) {
      owner_key[static_cast<std::size_t>(cost_interval_index(ll, rr))] = id;
    }
  );
}

double LossBlockCache::count_covered_intervals(
  const arma::imat & sawtooth_points,
  const int & l,
  const int & r
) const {
  double interval_writes = 0.0;
  for (unsigned int s = 0; s < sawtooth_points.n_rows; ++s) {
    const double left_count = l - sawtooth_points(s, 0) + 1.0;
    double right_count = 0.0;
    if (s == 0) {
      right_count = sawtooth_points(s, 1) - r + 1.0;
    } else {
      right_count = sawtooth_points(s, 1) - sawtooth_points(s - 1, 1);
    }
    interval_writes += left_count * right_count;
  }
  return interval_writes;
}

void LossBlockCache::ensure_exact_block(
  const int & key,
  const int & i,
  const int & j
) {
  if (exact_block_cache.find(key) != exact_block_cache.end()) {
    return;
  }

  LossBlockRecord block;
  block.l = i;
  block.r = j;
  block.l_end = i;
  block.r_end = j;
  block.loss = fit_loss_block(i, j, i, j);
  exact_block_cache.emplace(key, std::move(block));
  full_update_calls += 1.0;
}

arma::mat LossBlockCache::fit_loss_block(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  const double time_before = reg_loss.model_fit_time();
  arma::mat loss = reg_loss.block_loss(l, r, l_end, r_end);
  model_fit_time += reg_loss.model_fit_time() - time_before;
  raw_loss_block_cells += static_cast<double>(loss.n_rows) * loss.n_cols;
  return loss;
}

CostEngineByLossBlock::CostEngineByLossBlock(
  const int & cost_row,
  LossBlockCache & loss_block_cache
) : CostEngine(cost_row), loss_block_cache(loss_block_cache) {}

double CostEngineByLossBlock::get_cost(const int & i, const int & j) {
  return loss_block_cache.get_cost(cost_row, i, j);
}
