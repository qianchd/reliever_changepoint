#ifndef RELIEVER_RELIEF_INTERVAL_H
#define RELIEVER_RELIEF_INTERVAL_H

#include <armadillo>

arma::ivec exact2relief_itv_routine_c(
  const int & l,
  const int & r,
  const arma::ivec & miss_cover_len,
  const arma::ivec & int_len,
  const arma::ivec & layer_point,
  const arma::imat & int_eps
);

arma::imat relief2exact_itv_routine_c(
  const int & id,
  const arma::ivec & int_len,
  const arma::ivec & layer_point,
  const arma::imat & int_eps,
  const int & n
);

#endif
