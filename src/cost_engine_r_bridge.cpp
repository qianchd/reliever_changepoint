#include "cost_engine_r_bridge.h"
#include <chrono>

namespace {

double list_double_or_zero(const Rcpp::List & x, const char * name) {
  return x.containsElementNamed(name) ? Rcpp::as<double>(x[name]) : 0.0;
}

} // namespace

RRegLossFunction::RRegLossFunction(Rcpp::Function & individual_loss_fun)
  : individual_loss_fun(individual_loss_fun) {}

arma::mat RRegLossFunction::individual_loss(
  const int & l,
  const int & r,
  const int & l_end,
  const int & r_end
) {
  arma::mat loss;
  {
    auto tstart = std::chrono::steady_clock::now();
    Rcpp::NumericMatrix loss_r = Rcpp::as<Rcpp::NumericMatrix>(
      individual_loss_fun(l, r, l_end, r_end)
    );
    auto tend = std::chrono::steady_clock::now();
    model_fit_time_ += std::chrono::duration<double>(tend - tstart).count();
    loss = arma::mat(loss_r.begin(), loss_r.nrow(), loss_r.ncol(), true);
  }
  return loss;
}

double RRegLossFunction::model_fit_time() const {
  return model_fit_time_;
}

Rcpp::List loss_block_cache_state(
  const LossBlockCache & loss_block_cache
) {
  const LossBlockCacheState state = loss_block_cache.state();
  Rcpp::List blocks(state.blocks.size());
  for (std::size_t i = 0; i < state.blocks.size(); ++i) {
    const LossBlockStateRecord & record = state.blocks[i];
    const LossBlockRecord & block = record.block;
    blocks[i] = Rcpp::List::create(
      Rcpp::Named("key") = record.key,
      Rcpp::Named("loss") = block.loss,
      Rcpp::Named("l") = block.l,
      Rcpp::Named("r") = block.r,
      Rcpp::Named("l_end") = block.l_end,
      Rcpp::Named("r_end") = block.r_end
    );
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("n") = state.n,
    Rcpp::Named("blocks") = blocks,
    Rcpp::Named("raw_loss_block_cells") = state.raw_loss_block_cells,
    Rcpp::Named("get_cost_calls") = state.get_cost_calls,
    Rcpp::Named("expanded_cost_cells") = state.expanded_cost_cells,
    Rcpp::Named("expanded_interval_writes") =
      state.expanded_interval_writes,
    Rcpp::Named("full_update_calls") = state.full_update_calls,
    Rcpp::Named("relief_update_calls") = state.relief_update_calls,
    Rcpp::Named("model_fit_time") = state.model_fit_time
  );
  if (state.has_owner_key) {
    out["owner_key"] = Rcpp::IntegerVector(
      state.owner_key.begin(),
      state.owner_key.end()
    );
  }
  return out;
}

LossBlockCacheState loss_block_cache_state_from_r(
  const Rcpp::List & cache_state
) {
  LossBlockCacheState out;
  if (cache_state.containsElementNamed("n")) {
    out.n = Rcpp::as<int>(cache_state["n"]);
  }
  if (cache_state.containsElementNamed("blocks")) {
    Rcpp::List blocks = cache_state["blocks"];
    out.blocks.reserve(blocks.size());
    for (int i = 0; i < blocks.size(); ++i) {
      Rcpp::List block_r = blocks[i];
      LossBlockStateRecord record;
      record.key = Rcpp::as<int>(block_r["key"]);
      record.block.loss = Rcpp::as<arma::mat>(block_r["loss"]);
      record.block.l = Rcpp::as<int>(block_r["l"]);
      record.block.r = Rcpp::as<int>(block_r["r"]);
      record.block.l_end = Rcpp::as<int>(block_r["l_end"]);
      record.block.r_end = Rcpp::as<int>(block_r["r_end"]);
      out.blocks.push_back(record);
    }
  }
  if (cache_state.containsElementNamed("owner_key")) {
    Rcpp::IntegerVector owner_key_r = cache_state["owner_key"];
    out.has_owner_key = true;
    out.owner_key.assign(owner_key_r.begin(), owner_key_r.end());
  }
  out.raw_loss_block_cells =
    list_double_or_zero(cache_state, "raw_loss_block_cells");
  out.get_cost_calls = list_double_or_zero(cache_state, "get_cost_calls");
  out.expanded_cost_cells =
    list_double_or_zero(cache_state, "expanded_cost_cells");
  out.expanded_interval_writes =
    list_double_or_zero(cache_state, "expanded_interval_writes");
  out.full_update_calls =
    list_double_or_zero(cache_state, "full_update_calls");
  out.relief_update_calls =
    list_double_or_zero(cache_state, "relief_update_calls");
  out.model_fit_time = list_double_or_zero(cache_state, "model_fit_time");
  return out;
}

void warn_loss_block_exact_fallbacks(
  const LossBlockCache & loss_block_cache
) {
  const std::string warning = loss_block_cache.exact_fallback_warning();
  if (!warning.empty()) {
    Rcpp::warning("%s", warning.c_str());
  }
}
