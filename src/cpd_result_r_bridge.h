#ifndef RELIEVER_CPD_RESULT_R_BRIDGE_H
#define RELIEVER_CPD_RESULT_R_BRIDGE_H

#include <RcppArmadillo.h>
#include <vector>
#include "cpd_result.h"

Rcpp::List all_cpd_results_to_list(const AllCpdResults & results);

#endif
