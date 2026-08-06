#ifndef RELIEVER_COST_ENGINE_R_BRIDGE_H
#define RELIEVER_COST_ENGINE_R_BRIDGE_H

#include <RcppArmadillo.h>
#include "cost_engine.h"

class RRegLossFunction : public RegLossFunction {
public:
  explicit RRegLossFunction(Rcpp::Function & individual_loss_fun);
  arma::mat individual_loss(const int & l, const int & r,
                            const int & l_end, const int & r_end) override;
  double model_fit_time() const override;
private:
  Rcpp::Function & individual_loss_fun;
  double model_fit_time_ = 0.0;
};

Rcpp::List loss_block_cache_state(const LossBlockCache & loss_block_cache);
LossBlockCacheState loss_block_cache_state_from_r(
  const Rcpp::List & cache_state
);
void warn_loss_block_exact_fallbacks(
  const LossBlockCache & loss_block_cache
);

#endif
