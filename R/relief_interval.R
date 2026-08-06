# Core relief interval constructors and C++ wrappers used by reliever.
#' Create deterministic relief intervals
#'
#' Construct the representative intervals on which Reliever fits models when
#' \code{cov_rate < 1}. The intervals are arranged in layers of increasing
#' length so that an exact interval requested by a changepoint algorithm can be
#' mapped to a containing representative interval. Endpoints use the half-open
#' convention \code{(l, r]}, meaning observations \code{l + 1} through
#' \code{r}; consequently, \code{l} may be 0 and \code{r} may be \code{n}.
#'
#' @param n Number of observations.
#' @param cov_rate Requested Reliever coverage rate in `(0, 1]`.
#' @param dm Smallest exact-interval length that must be covered. This normally
#'   equals the minimum segment length supplied to \code{reliever()}.
#'
#' @return A list describing the layered interval set:
#'   \itemize{
#'     \item \code{int_eps}: two-column matrix of all \code{(l, r]} endpoints,
#'     ordered by layer.
#'     \item \code{layer_point}: cumulative row number at which each layer ends
#'     in \code{int_eps}.
#'     \item \code{int_len}: representative-interval length in each layer.
#'     \item \code{target_cover_len}: exact-interval length from which each
#'     layer is designed to provide the requested coverage.
#'     \item \code{miss_cover_len}: largest exact-interval length that may still
#'     fail to contain a representative interval from that layer; the mapping
#'     routines use this cutoff to choose the first eligible layer.
#'     \item \code{n}: number of observations.
#'   }
#' @seealso \code{\link{reliever_generic}()},
#'   \code{\link{reliever_mean}()}, \code{\link{create_wbs_itv}()},
#'   \code{\link{create_seed_itv}()}
#' @export
#'
#' @examples
#' int_set <- create_relief_itv(n = 900, cov_rate = 0.8, dm = 30)
#' head(int_set$int_eps)
#' int_set$target_cover_len
create_relief_itv <- function(n, cov_rate, dm) {
  n <- .reliever_validate_positive_integer(n, "n")
  dm <- .reliever_validate_positive_integer(dm, "dm")
  if (dm > n) {
    stop("dm must be no larger than n.", call. = FALSE)
  }
  if (!is.numeric(cov_rate) || length(cov_rate) != 1 || is.na(cov_rate) || cov_rate <= 0 || cov_rate > 1) {
    stop("cov_rate must be a single number in (0, 1].", call. = FALSE)
  }
  if (cov_rate >= (n - 1) / n && cov_rate <= 1) {
    warning(
      "cov_rate >= (n - 1) / n triggers a full search; use cov_rate < ",
      format((n - 1) / n),
      " when a partial Reliever search is intended.",
      call. = FALSE
    )
  }
  int_len <- max(floor(dm * sqrt(cov_rate)), ceiling(dm * cov_rate))
  target_cover_len <- dm
  int_old <- 0

  i <- 1
  int_endpoints <- matrix(0, 0, 2)
  int_len_vec <- layer_point <- miss_cover_len_vec <- target_cover_len_vec <- 0
  while(target_cover_len <= n) { # note that (0, n] is always not missed.
    miss_cover_len <- target_cover_len - 1
    wriggling_size <- target_cover_len - int_len + 1

    l_list <- seq(0, n - int_len - wriggling_size + 1, wriggling_size)
    r_list <- l_list + int_len
    end_size <- floor((n - r_list[length(r_list)]) / 2)
    if (end_size > 0) {
      l_list <- l_list + end_size
      r_list <- r_list + end_size
    }
    int_endpoints <- rbind(int_endpoints, cbind(l_list, r_list))

    if(length(l_list) == 1) {
      miss_cover_len <- int_len + max(l_list, n - r_list) - 1
    }
    int_len_vec[i] <- int_len
    layer_point[i] <- nrow(int_endpoints)
    miss_cover_len_vec[i] <- miss_cover_len
    target_cover_len_vec[i] <- target_cover_len
    i <- i + 1

    target_cover_next <- floor(int_len / cov_rate) + 1
    int_next <- max(
      round(int_len / sqrt(cov_rate)),
      int_len + 1,
      floor(int_old / cov_rate) + 1
    )
    int_old <- int_len
    int_len <- int_next
    target_cover_len <- target_cover_next
  }
  mode(int_endpoints) <- mode(layer_point) <- mode(int_len_vec) <- mode(miss_cover_len_vec) <- mode(target_cover_len_vec) <- 'integer'
  colnames(int_endpoints) <- c("l", "r")
  int_set <- list(
    int_eps = int_endpoints,
    layer_point = layer_point,
    int_len = (int_len_vec),
    miss_cover_len = (miss_cover_len_vec),
    target_cover_len = target_cover_len_vec,
    n = as.integer(n)
  )
  return(int_set)
}

#' Create random WBS intervals
#'
#' Create random intervals for wild binary segmentation style searches. The
#' returned matrix has two columns \code{l} and \code{r}; each row represents
#' observations \code{l + 1} through \code{r}. Intervals shorter than
#' \code{dm} are discarded, so the function can return fewer than \code{M}
#' rows.
#'
#' @param n Number of observations.
#' @param dm Minimum retained interval length.
#' @param M Maximum number of retained intervals.
#' @param wbs_seed Optional seed used for this interval sample.
#'
#' @return An integer matrix with columns `l` and `r`.
#' @seealso \code{\link{create_seed_itv}()},
#'   \code{\link{twostep}()}, \code{\link{reliever_generic}()}
#' @export
#'
#' @examples
#' create_wbs_itv(n = 900, dm = 30, M = 20, wbs_seed = 2026)
create_wbs_itv <- function(n, dm, M, wbs_seed = NULL) {
  n <- .reliever_validate_positive_integer(n, "n")
  dm <- .reliever_validate_positive_integer(dm, "dm")
  M <- .reliever_validate_positive_integer(M, "M", allow_zero = TRUE)
  if (M == 0) {
    return(matrix(
      integer(), nrow = 0L, ncol = 2L,
      dimnames = list(NULL, c("l", "r"))
    ))
  }
  lr_m <- .reliever_with_seed(wbs_seed, {
    lr_m <- matrix(NA, 2 * M, 2)
    for (i in seq_len(nrow(lr_m))) {
      lr_m[i, ] <- sort(sample(0:n, 2, replace = FALSE))
    }
    lr_m
  })
  lr_m <- lr_m[which((lr_m[, 2] - lr_m[, 1]) >= dm), , drop = FALSE]
  if (nrow(lr_m) > M) {
    lr_m <- lr_m[seq_len(M), , drop = FALSE]
  }
  mode(lr_m) <- 'integer'
  colnames(lr_m) <- c("l", "r")
  lr_m
}

#' Create deterministic SeedBS intervals
#'
#' Create the multiscale deterministic intervals used by
#' \code{method = "SeedBS"}. The returned matrix has columns \code{l} and
#' \code{r}; each row represents observations \code{l + 1} through \code{r}.
#'
#' @param n Number of observations.
#' @param dm Target width of the base partition used to build the multiscale
#'   seed intervals. The shortest returned intervals span about two base
#'   blocks, rather than exactly \code{dm} observations.
#'
#' @return An integer matrix with columns `l` and `r`.
#' @seealso \code{\link{create_wbs_itv}()},
#'   \code{\link{twostep}()}, \code{\link{reliever_generic}()}
#' @export
#'
#' @examples
#' head(create_seed_itv(n = 900, dm = 30))
create_seed_itv <- function(n, dm) {
  n <- .reliever_validate_positive_integer(n, "n")
  dm <- .reliever_validate_positive_integer(dm, "dm")
  if (dm > n) {
    stop("dm must be no larger than n.", call. = FALSE)
  }
  n_grid <- floor(n / dm)
  dm <- floor(n / n_grid)
  int_len <- rep(dm, n_grid)

  remainder <- n - sum(int_len)
  if (remainder > 0L) {
    ind_expand <- floor(
      (seq_len(remainder) - 0.5) * n_grid / remainder
    ) + 1L
    int_len[ind_expand] <- int_len[ind_expand] + 1L
  }
  cum_int_len <- c(0, cumsum(int_len))
  seed_int <- matrix(0, 0, 2)
  if (n_grid < 2L) {
    seed_int <- matrix(c(0L, n), nrow = 1L)
    mode(seed_int) <- 'integer'
    colnames(seed_int) <- c("l", "r")
    return(seed_int)
  }
  s <- 1
  while (n_grid >= 2 * s) {
    if (n_grid %% s == 0) {
      seed_int <- rbind(
        seed_int,
        cbind(
          cum_int_len[seq(1, n_grid + 1 - 2 * s, s)],
          cum_int_len[seq(2 * s + 1, n_grid + 1, s)]
        )
      )
    } else {
      seed_int <- rbind(
        seed_int,
        cbind(
          cum_int_len[seq(1, n_grid + 1 - 2 * s, s)],
          cum_int_len[seq(2 * s + 1, n_grid + 1, s)]
        ),
        c(cum_int_len[n_grid + 1 - 2 * s], cum_int_len[n_grid + 1])
      )
    }
    if (n_grid >= 3 * s) {
      if (n_grid %% s == 0) {
        seed_int <- rbind(
          seed_int,
          cbind(
            cum_int_len[seq(1, n_grid + 1 - 3 * s, s)],
            cum_int_len[seq(3 * s + 1, n_grid + 1, s)]
          )
        )
      } else {
        seed_int <- rbind(
          seed_int,
          cbind(
            cum_int_len[seq(1, n_grid + 1 - 3 * s, s)],
            cum_int_len[seq(3 * s + 1, n_grid + 1, s)]
          ),
          c(
            cum_int_len[n_grid + 1 - 3 * s],
            cum_int_len[n_grid + 1]
          )
        )
      }
    }
    s <- 2 * s
  }
  if (!(seed_int[nrow(seed_int), 1] == 0 && seed_int[nrow(seed_int), 2] == n)) {
    seed_int <- rbind(seed_int, c(0, n))
  }
  mode(seed_int) <- 'integer'
  colnames(seed_int) <- c("l", "r")
  seed_int
}
