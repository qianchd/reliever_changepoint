.reliever_run_native_mean <- function(data, cpn_max, dm, cov_rate,
                                      method, pen_val, prune_value,
                                      M, wbs_seed, cpn_crit,
                                      wbs_stop_crit,
                                      cache_backend, owner_key, ratio, detail,
                                      echo) {
  detail_mode <- .reliever_detail_mode(detail)
  return_cache_profile <- detail_mode != "none"

  setup <- .reliever_prepare_search_setup(
    data = data,
    cpn_max = cpn_max,
    dm = dm,
    pen_val = pen_val,
    cov_rate = cov_rate,
    method = method,
    M = M,
    prune_value = prune_value,
    dc_grid = NULL,
    cache_backend = cache_backend,
    wbs_stop_crit = wbs_stop_crit
  )
  cpn_penalty_info <- .cpn_penalty(cpn_crit, setup$n)

  search_intervals <- .reliever_search_intervals(
    setup$method, setup$n, setup$dm, setup$M, wbs_seed
  )
  cost_mat <- NULL
  int_args <- if (setup$is_full) {
    list(
      miss_cover_len = integer(),
      int_len = integer(),
      layer_point = integer(),
      int_eps = matrix(integer(), 0L, 2L)
    )
  } else {
    setup$int_set[c("miss_cover_len", "int_len", "layer_point", "int_eps")]
  }
  if (setup$cache_backend == "by_cost_mat") {
    cost_mat <- matrix(Inf, 1L, 1L + setup$n * (setup$n + 1L) / 2L)
  }

  if (echo) timer_start <- proc.time()[["elapsed"]]
  native_res <- if (setup$cache_backend == "by_loss_block") {
    cpd_mean_by_loss_block_cpp(
      setup$method, setup$n, setup$cpn_max, setup$dm, search_intervals,
      pen_val, setup$prune_value, data, ratio,
      int_args$miss_cover_len, int_args$int_len,
      int_args$layer_point, int_args$int_eps, return_cache_profile,
      owner_key
    )
  } else {
    cpd_mean_by_cost_mat_cpp(
      setup$method, setup$n, setup$cpn_max, setup$dm, search_intervals,
      pen_val, setup$prune_value, data, cost_mat, ratio,
      int_args$miss_cover_len, int_args$int_len,
      int_args$layer_point, int_args$int_eps, setup$is_full
    )
  }
  if (echo) {
    message(sprintf(
      "reliever_mean-%s: %.3f sec elapsed",
      setup$cache_backend,
      proc.time()[["elapsed"]] - timer_start
    ))
  }

  cache_profile <- NULL
  if (return_cache_profile) {
    cache_profile <- if (setup$cache_backend == "by_loss_block") {
      .reliever_by_loss_block_cache_profile(
        native_res$cache_profile, setup$int_set
      )
    } else {
      .reliever_by_cost_mat_cache_profile(
        cost_mat, setup$int_set
      )
    }
  }

  cpd_result <- .reliever_format_backend_cpd_result(
    backend_res = native_res,
    setup = setup,
    pen_val = pen_val
  )
  model_fit_time <- as.numeric(native_res$model_fit_time)
  total_time <- as.numeric(native_res$total_time)
  num_pruned <- if (setup$method %in% c("PELT", "OP")) {
    native_res$num_pruned
  } else {
    list()
  }

  out <- list(
    cpd_path = cpd_result$cpd_path,
    run_meta = .reliever_run_meta(
      data.frame(
        loss_output_id = 1L, row_type = "mean", loss_kind = "rss"
      ),
      1L
    ),
    n_model_fit = as.numeric(native_res$n_model_fit),
    model_fit_time = model_fit_time,
    cpd_time = total_time - model_fit_time,
    total_time = total_time,
    diagnostics = if (detail_mode == "none") {
      NULL
    } else {
      list(path_score = cpd_result$path_score, num_pruned = num_pruned)
    },
    cache_profile = cache_profile,
    settings = .reliever_result_settings(
      setup = setup,
      cpn_crit = cpn_crit,
      cpn_penalty_info = cpn_penalty_info,
      pen_val = pen_val,
      wbs_seed = wbs_seed,
      owner_key = owner_key,
      ratio = ratio
    ),
    loss_spec = list(reg_fun = reg_fun_mean, args = list())
  )
  .reliever_finalize_result(out)
}

#' Mean-square interval-loss function
#'
#' Estimate the segment mean from rows \code{l:r} and return one mean-square
#' loss per row in \code{l_end:r_end}. For multivariate observations, each loss
#' is the average squared error across columns. This R interval-loss function is intended for
#' \code{reliever_generic()} and \code{cv.reliever_generic()}. For ordinary
#' mean-change detection without outer cross-validation,
#' \code{reliever(X, cpd_family = "mean")} dispatches to the faster native C++
#' implementation.
#'
#' @param data Numeric vector or matrix with observations in rows.
#' @param l,r Rows used to estimate the segment mean.
#' @param l_end,r_end Rows whose mean-square losses are returned.
#' @param save_model Return the fitted segment mean as
#'   \code{list(center = ...)} when \code{TRUE}.
#' @param is_virtual_run Query flag used by \code{reliever_generic()}. When
#'   \code{TRUE}, return metadata for the single loss output without estimating
#'   a mean.
#'
#' @return With \code{is_virtual_run = TRUE}, metadata for one loss output.
#'   Otherwise, a list containing a one-column \code{loss} matrix and
#'   \code{model}. When \code{save_model = TRUE}, \code{model} is
#'   \code{list(center = ...)}; otherwise it is \code{NULL}. See the
#'   \emph{Writing a custom reg_fun} section of
#'   \code{\link{reliever_generic}()} for the common return convention.
#' @seealso \code{\link{reliever_mean}()},
#'   \code{\link{cv.reliever_generic}()}, \code{\link{reliever_generic}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' out <- reg_fun_mean(x, l = 241, r = 660, l_end = 211, r_end = 690)
#' dim(out$loss)
reg_fun_mean <- function(data, l, r, l_end = l, r_end = r,
                         save_model = FALSE, is_virtual_run = FALSE) {
  if (is_virtual_run) {
    return(list(
      n_loss_outputs = 1L,
      loss_output_meta = data.frame(
        loss_output_id = 1L, row_type = "mean", loss_kind = "rss"
      )
    ))
  }
  .reg_fun_mean_loss(
    data = data, l = l, r = r, l_end = l_end, r_end = r_end,
    save_model = save_model
  )
}

#' Reliever for mean-change loss
#'
#' A focused wrapper for mean-change detection. Ungridded searches use the
#' package's native C++ mean-loss engine. For a gridded search, prefer
#' \code{dc_grid_size}; the same mean loss then runs through the generic
#' gridded backend and returns changepoints on the original time scale. The
#' main built-in entry is \code{reliever(X)}; this wrapper exposes all
#' mean-specific arguments.
#'
#' @section Selection:
#' Use \code{cv.reliever()} for CPSS-style outer sample-splitting selection of
#' K. As an information-criterion alternative,
#' \code{select_by_run(fit, cpn_crit = "rss_sic")} minimizes
#' \eqn{n\log(RSS/n)/2 + \log(n)K} on an existing fit.
#'
#' @param data A numeric matrix or vector with observations in rows.
#' @param cpn_max,dm,cov_rate,method,cpn_crit Common path controls; see
#'   \code{\link{reliever}()}.
#' @param pen_val,prune_value,M,wbs_seed,wbs_stop_crit Common penalty and
#'   search controls; see \code{\link{reliever}()}.
#' @param detail,cache_backend,owner_key,echo Common output and cache controls;
#'   see \code{\link{reliever}()}.
#' @param dc_grid_size,dc_grid Common candidate-grid controls; see
#'   \code{\link{reliever}()}.
#' @param ratio Computational tuning parameter in \code{(0, 1]}. When the next
#'   fitted interval differs from the previous one in no more than
#'   \code{ratio} times the new interval length, the C++ engine updates its
#'   mean and sum of squares incrementally; otherwise it recomputes them from
#'   the new interval. This affects speed, not the fitted criterion, and is
#'   ignored for a gridded search.
#'
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @seealso \code{\link{reliever}()}, \code{\link{cv.reliever}()},
#'   \code{\link{reliever_mean_crossfit}()},
#'   \code{\link{reg_fun_mean}()},
#'   \code{\link{select_by_run}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' res <- reliever_mean(data = x, cpn_max = 7, dm = 30, cov_rate = 0.8,
#'                      method = "SN", echo = FALSE)
#' summary(res) # empty: fitting did not choose K
#'
#' # Optional information-criterion selection from the stored path:
#' selected_sic <- select_by_run(result = res, cpn_crit = "rss_sic")
#' selected_sic
#' stopifnot(identical(selected_sic$K_hat, 2L))
#' res$cpd_path$candidates
reliever_mean <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                          method = "SN", cpn_crit = "none",
                          pen_val = 1, prune_value = 0,
                          M = 100, wbs_seed = NULL,
                          wbs_stop_crit = NULL,
                          detail = FALSE,
                          cache_backend = "by_loss_block",
                          owner_key = TRUE,
                          echo = FALSE,
                          dc_grid_size = NULL,
                          dc_grid = NULL,
                          ratio = 0.9) {
  if (!is.null(dc_grid_size) || !is.null(dc_grid)) {
    .reliever_resolve_dc_grid(
      data = data,
      dc_grid_size = dc_grid_size,
      dc_grid = dc_grid
    )
    if (!missing(ratio)) {
      warning("ratio is ignored for a gridded search.", call. = FALSE)
    }
    return(reliever_generic(
      data = data,
      reg_fun = reg_fun_mean,
      cpn_max = cpn_max,
      dm = dm,
      cov_rate = cov_rate,
      method = method,
      cpn_crit = cpn_crit,
      pen_val = pen_val,
      prune_value = prune_value,
      M = M,
      wbs_seed = wbs_seed,
      wbs_stop_crit = wbs_stop_crit,
      detail = detail,
      cache_backend = cache_backend,
      owner_key = owner_key,
      echo = echo,
      dc_grid_size = dc_grid_size,
      dc_grid = dc_grid
    ))
  }
  method <- match.arg(
    method,
    c("SN", "WBS", "WBS_recursive", "SeedBS", "BS", "PELT", "OP")
  )
  cache_backend <- match.arg(
    cache_backend, c("by_loss_block", "by_cost_mat")
  )
  if (!is.logical(owner_key) || length(owner_key) != 1L ||
      is.na(owner_key)) {
    stop("owner_key must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(ratio) || length(ratio) != 1 || is.na(ratio) ||
      ratio <= 0 || ratio > 1) {
    stop("ratio must be a single number in (0, 1].", call. = FALSE)
  }
  if (is.vector(data)) {
    data <- matrix(data, ncol = 1L)
  } else {
    data <- as.matrix(data)
  }
  if (!is.numeric(data) || nrow(data) < 1L || ncol(data) < 1L || anyNA(data) ||
      any(!is.finite(data))) {
    stop("data must contain finite numeric observations.", call. = FALSE)
  }

  .reliever_run_native_mean(
    data = data,
    cpn_max = cpn_max,
    dm = dm,
    cov_rate = cov_rate,
    method = method,
    pen_val = pen_val,
    prune_value = prune_value,
    M = M,
    wbs_seed = wbs_seed,
    wbs_stop_crit = wbs_stop_crit,
    cpn_crit = cpn_crit,
    cache_backend = cache_backend,
    owner_key = owner_key,
    ratio = ratio,
    detail = detail,
    echo = echo
  )
}
