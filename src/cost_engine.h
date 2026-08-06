#ifndef RELIEVER_COST_ENGINE_H
#define RELIEVER_COST_ENGINE_H

#include <RcppArmadillo.h>
#include <string>
#include <unordered_map>
#include <vector>

class RegLossFunction {
public:
  virtual ~RegLossFunction() {}
  virtual arma::mat individual_loss(const int & l, const int & r,
                                    const int & l_end, const int & r_end);
  virtual arma::mat block_loss(const int & l, const int & r,
                               const int & l_end, const int & r_end);
  virtual double model_fit_time() const;
};

class MeanSquareRegLossFunction : public RegLossFunction {
public:
  MeanSquareRegLossFunction(const double & ratio,
                            const arma::mat & data);
  arma::mat individual_loss(const int & l, const int & r,
                            const int & l_end, const int & r_end) override;
  arma::mat block_loss(const int & l, const int & r,
                       const int & l_end, const int & r_end) override;
  double model_fit_time() const override;
private:
  void fit_core_interval(const unsigned int & l, const unsigned int & r);
  void live_update(const unsigned int & l, const unsigned int & r);
  void refresh_core_loss_assessment();
  void stabilize_core_mean();
  arma::vec stable_core_mean() const;
  arma::vec squared_residual_loss(const arma::uvec & row_idx) const;
  double stable_core_loss() const;

  const arma::mat & data;
  const arma::vec sq_smry;
  arma::vec mu_hat;
  double model_fit_time_ = 0.0;
  double squ = 0.0;
  double core_fast_loss = 0.0;
  double ratio = 0.9;
  bool core_loss_risky = true;
  unsigned int l_current = 0;
  unsigned int r_current = 0;
  unsigned int len_current = 0;
};

class CostEngine {
public:
  explicit CostEngine(const int & cost_row_);
  virtual ~CostEngine() {}
  virtual double get_cost(const int & i, const int & j) = 0;
  double get_model_fit_time();
  int get_n_model_fit();
protected:
  const int cost_row;
  int n_model_fit = 0;
  double model_fit_time = 0.0;
};

class CostMatCache {
public:
  CostMatCache(arma::mat & cost_mat,
               const arma::ivec & miss_cover_len,
               const arma::ivec & int_len,
               const arma::ivec & layer_point,
               const arma::imat & int_eps,
               const int & n);

  double cached_cost(const int & row, const int & col);
  void write_cached_cost(const int & col, const arma::vec & value);
  arma::ivec find_relief_interval(const int & l, const int & r);
  arma::imat covered_intervals(const int & id);
private:
  arma::mat & cost_mat;
  const arma::ivec & miss_cover_len;
  const arma::ivec & int_len;
  const arma::ivec & layer_point;
  const arma::imat & int_eps;
  const int & n;
};

class CostEngineByCostMat : public CostEngine {
public:
  CostEngineByCostMat(CostMatCache & cost_mat_cache,
                      RegLossFunction & reg_loss,
                      const int & cost_row,
                      const bool & is_full_search);
  double get_cost(const int & i, const int & j) override;
protected:
  CostMatCache & cost_mat_cache;
  RegLossFunction & reg_loss;
  bool is_full_search;

  virtual void update_cache(const int & i, const int & j);
  void update_exact_cache(const int & i, const int & j);
  arma::mat fit_loss_block(const int & l, const int & r,
                           const int & l_end, const int & r_end);
};

struct LossBlockRecord {
  arma::mat loss;
  int l;
  int r;
  int l_end;
  int r_end;
};

struct LossBlockStateRecord {
  int key = 0;
  LossBlockRecord block;
};

struct LossBlockCacheState {
  int n = 0;
  std::vector<LossBlockStateRecord> blocks;
  bool has_owner_key = false;
  std::vector<int> owner_key;
  double raw_loss_block_cells = 0.0;
  double get_cost_calls = 0.0;
  double expanded_cost_cells = 0.0;
  double expanded_interval_writes = 0.0;
  double full_update_calls = 0.0;
  double relief_update_calls = 0.0;
  double model_fit_time = 0.0;
};

class LossBlockCache {
public:
  LossBlockCache(RegLossFunction & reg_loss,
                 const arma::ivec & miss_cover_len,
                 const arma::ivec & int_len,
                 const arma::ivec & layer_point,
                 const arma::imat & int_eps,
                 const int & n,
                 const bool & use_owner_key);

  double get_cost(const int & cost_row, const int & i, const int & j);
  void restore(const LossBlockCacheState & cache_state);
  LossBlockCacheState state() const;
  std::string exact_fallback_warning() const;

  double raw_loss_block_cells = 0.0;
  double get_cost_calls = 0.0;
  double expanded_cost_cells = 0.0;
  double expanded_interval_writes = 0.0;
  double exact_fallback_calls = 0.0;
  double full_update_calls = 0.0;
  double relief_update_calls = 0.0;
  double model_fit_time = 0.0;

private:
  RegLossFunction & reg_loss;
  const arma::ivec & miss_cover_len;
  const arma::ivec & int_len;
  const arma::ivec & layer_point;
  const arma::imat & int_eps;
  const int & n;
  bool use_owner_key;
  std::vector<LossBlockRecord> relief_blocks;
  std::vector<unsigned char> relief_block_ready;
  std::vector<unsigned char> owner_key_ready;
  std::vector<int> owner_key;
  std::unordered_map<int, LossBlockRecord> exact_block_cache;
  std::vector<int> exact_fallback_i;
  std::vector<int> exact_fallback_j;

  int locate_block_key(const int & i, const int & j,
                       const long long & interval_index);
  const LossBlockRecord & cached_block(const int & key) const;
  void ensure_relief_block(const int & id);
  void fill_owner_keys(const int & id,
                       const arma::imat & sawtooth_points,
                       const int & l,
                       const int & r);
  double count_covered_intervals(const arma::imat & sawtooth_points,
                                 const int & l,
                                 const int & r) const;
  void ensure_exact_block(const int & key, const int & i, const int & j);
  arma::mat fit_loss_block(const int & l, const int & r,
                           const int & l_end, const int & r_end);
};

class CostEngineByLossBlock : public CostEngine {
public:
  CostEngineByLossBlock(const int & cost_row, LossBlockCache & loss_block_cache);
  double get_cost(const int & i, const int & j) override;
private:
  LossBlockCache & loss_block_cache;
};

#endif
