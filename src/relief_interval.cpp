// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "relief_interval.h"
#include <algorithm>
#include <stdexcept>

// [[Rcpp::export]]
arma::ivec exact2relief_itv_routine_c(const int & l, const int & r, const arma::ivec & miss_cover_len,
                       const arma::ivec & int_len,
                       const arma::ivec & layer_point,
                       const arma::imat & int_eps) {
  int len = r - l;

  // Locate the lowest interval layer whose missing length can cover this
  // exact interval length.
  unsigned int k = std::upper_bound(miss_cover_len.begin(), miss_cover_len.end(), len - 1) - miss_cover_len.begin();
  if (k == 0) {
    throw std::runtime_error("Reliever interval not found for (" + std::to_string(l) + "," + std::to_string(r) + "]");
  }
  k -= 1;

  // Restrict the search to layers whose relief interval length is no larger
  // than the queried exact interval length.
  const unsigned int kmax = std::upper_bound(int_len.begin(), int_len.end(), len) - int_len.begin() - 1;

  int id;
  int left, right, id_lend, id_rend, it_loc;

  // layer_point stores cumulative counts from R. Use it to form an inclusive
  // 0-based C++ index range in int_eps.
  id_lend = (k == 0) ? 0 : layer_point(k - 1);
  id_rend = layer_point(kmax) - 1;

  // Search from later intervals to match relief2exact_itv_routine_c's ownership rule.
  it_loc = -1;
  for (int i = id_rend; i >= id_lend; --i) {
    if (int_eps(i, 1) <= r && int_eps(i, 0) >= l) {
      it_loc = i;
      break;
    }
  }

  if (it_loc == -1) {
    throw std::runtime_error("No valid interval found within the specified range");
  }

  id = it_loc;
  left = int_eps(id, 0);
  right = int_eps(id, 1);
  return arma::ivec({id + 1, left, right}); // Return id as R's 1-based index.
}

// [[Rcpp::export]]
arma::imat relief2exact_itv_routine_c(const int & id,
  const arma::ivec & int_len,
  const arma::ivec & layer_point,
  const arma::imat & int_eps, const int & n) {

    // Find the first layer point that is >= id
    unsigned int k = std::lower_bound(layer_point.begin(), layer_point.end(), id) - layer_point.begin();

    // Get the interval for the given id. The input id is R's 1-based index.
    const int it_left = int_eps(id - 1, 0);
    const int it_right = int_eps(id - 1, 1);

    if (k + 1 == layer_point.n_elem) {
        if (id < layer_point(k)) {
            // Later intervals in the same last layer win ties in exact2relief_itv_routine_c.
            // Stop before the next interval can also contain the expanded segment.
            const int next_right = int_eps(id, 1) - 1;
            return arma::imat({{1, next_right}});
        }
        return arma::imat({{1, n}});
    }

    if (k > 0 && layer_point(k) - layer_point(k - 1) == 1) {
        // Special case where the layer point difference is 1
        const int next_left = int_eps(id, 0) + 2; // +1 because int_eps saves (l,r], converting to the matrix index will be [l+1, r]; and the other +1 is for the reason that the next interval covers l+1.
        const int next_right = int_eps(id, 1) - 1;
        return arma::imat({{1, next_right}, {next_left, n}});
    }

    // Get the upper intervals
    const int start_idx = layer_point(k);
    const int end_idx = layer_point(k + 1);
    const arma::imat it_up = int_eps.rows(start_idx, end_idx - 1);

    // Calculate wriggle sizes
    int wriggle_size = int_eps(layer_point(k) - 1, 0) - int_eps(layer_point(k) - 2, 0);
    int wriggle_size_up = -1;
    if (k + 2 < layer_point.n_elem && layer_point(k + 1) - layer_point(k) > 1) {
        wriggle_size_up = int_eps(layer_point(k + 1) - 1, 0) - int_eps(layer_point(k + 1) - 2, 0);
    }

    if (wriggle_size_up == 1) {
        // Special case for wriggle size up
        int new_left = std::max(1, it_left + int_len(k) + 2 - int_len(k + 1));
        return arma::imat({{new_left, it_right}});
    }
    // Find the last interval where it_right >= it_up[, 2]
    int loc = -1;
    for (unsigned int i = 0; i < it_up.n_rows; ++i) {
        if (it_right >= it_up(i, 1)) {
            loc = i;
        }
    }

    int i0 = (loc == -1) ? 1 : it_up(loc, 0) + 2;
    // Find intervals in the middle
    arma::uvec loc_mid = arma::find((it_up.col(1) > it_right) && (it_up.col(1) < (it_right + wriggle_size)));
    arma::ivec ilist(loc_mid.n_elem + 1);
    arma::uvec id_tem = {0};
    if(loc_mid.n_elem > 0) ilist.subvec(1, loc_mid.n_elem) = it_up(loc_mid, id_tem) + 2;
    ilist(0) = i0;
    arma::ivec jlist(loc_mid.n_elem + 1);
    id_tem(0) = 1;
    if(loc_mid.n_elem > 0) jlist.subvec(0, loc_mid.n_elem - 1) = it_up(loc_mid, id_tem) - 1;
    jlist(loc_mid.n_elem) = std::min(n + 1, it_right + wriggle_size) - 1;
    // Remove intervals where it_left + 1 < ilist
    int s = 0;
    for(unsigned int ss = 0; ss < ilist.n_elem; ss++) {
      if(ilist(ss) > it_left + 1) {
        s = ss + 1;
        break;
      }
    }
    if (s > 0) {
        ilist = ilist.subvec(0, s - 2);
        jlist = jlist.subvec(0, s - 2);
    }
    // Combine results
    return arma::join_horiz(ilist, jlist);
}
