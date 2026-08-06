#ifndef RELIEVER_REG_FUN_H
#define RELIEVER_REG_FUN_H

#include <armadillo>

arma::vec nmcd_individual_loss(const arma::vec & x,
                               int l,
                               int r,
                               int l_end,
                               int r_end,
                               const arma::vec & sorted_reference,
                               int tail_truncation = 0);

#endif
