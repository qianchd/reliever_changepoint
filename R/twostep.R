#' Two-step changepoint search
#'
#' Alternative to \code{reliever()} for WBS-family searches. Within each search
#' interval, fit left and right models at only a small set of initial split
#' fractions. Their observation-level losses are then reused to score every
#' admissible split location. This can reduce model fitting when a custom model
#' is expensive, at the cost of using the two-step approximation rather than
#' exact interval fits.
#'
#' @param data Data object accepted by \code{reg_fun}.
#' @param reg_fun Custom interval-loss function following the generic fitting
#'   contract. It must return observation-level losses;
#'   a fast block-loss-only C++ implementation is not sufficient.
#' @param cpn_max Maximum number of changepoints.
#' @param dm Minimum segment length.
#' @param num_init Number of equally spaced initial split fractions used when
#'   \code{init_cand = NULL}.
#' @param init_cand Optional numeric vector of initial split fractions in
#'   \code{(0, 1)}. Within each search interval, a fraction is mapped to the
#'   nearest split that leaves at least \code{dm} observations on both sides;
#'   duplicate mapped splits are fitted once.
#' @param method Search algorithm. Supported methods are `"WBS"`,
#'   `"WBS_recursive"`, `"SeedBS"`, and `"BS"`.
#'   \code{"WBS"} uses Reliever's global greedy WBS path: after each split, it
#'   compares all current child segments and splits the segment with the
#'   largest gain. \code{"WBS_recursive"} follows the original recursive WBS
#'   style used by the WBS package.
#' @param M Maximum retained random intervals for \code{"WBS"} and
#'   \code{"WBS_recursive"}. Ignored by \code{"SeedBS"} and \code{"BS"}.
#' @param wbs_seed Optional seed for random WBS interval generation.
#' @param wbs_stop_crit Optional stopping thresholds. For each threshold, keep
#'   the splits before the first proposed split whose gain is no larger than the
#'   threshold. The search remains limited to \code{cpn_max} splits.
#' @param echo Print timing messages.
#' @param detail When \code{TRUE}, retain \code{gain_mat}, the best gain found
#'   for each loss output and search interval, and \code{split_mat}, the
#'   corresponding split locations.
#' @param ... Additional arguments passed to \code{reg_fun}.
#'
#' @return A \code{reliever_result} object with the common components
#'   documented in \code{\link{reliever}()}. With \code{wbs_stop_crit},
#'   \code{summary()} reports the segmentation for each threshold and loss
#'   output. Otherwise the summary is empty because this approximation has no
#'   whole-segmentation loss for selecting K; the candidates remain in
#'   \code{cpd_path$candidates}. Use \code{select_holdout()} when independent
#'   evaluation observations are available.
#' @references Kaul, A., Jandhyala, V. K., and Fotopoulos, S. B. (2019). An
#'   efficient two step algorithm for high dimensional change point regression
#'   models without grid search. \emph{Journal of Machine Learning Research},
#'   20(111), 1--40.
#' @seealso \code{\link{reliever_generic}()},
#'   \code{\link{select_holdout}()}, and \code{\link{local_refine}()}.
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
#' fit <- twostep(x, cpn_max = 7, dm = 30, reg_fun = reg_fun_mean,
#'                method = "SeedBS", echo = FALSE)
#' fit$cpd_path$candidates
twostep <- function(data, reg_fun, cpn_max = 3, dm = 50,
                    num_init = 3, init_cand = NULL,
                    method = "WBS", M = 100, wbs_seed = NULL,
                    wbs_stop_crit = NULL, echo = FALSE, detail = FALSE, ...) {
  if (missing(reg_fun)) {
    stop("reg_fun must be supplied to twostep().", call. = FALSE)
  }
  para_list <- list(...)
  .reliever_reject_renamed_cpn_max(para_list, "twostep")
  method <- match.arg(method, c("WBS", "WBS_recursive", "SeedBS", "BS"))
  wbs_stop_crit <- .reliever_validate_wbs_stop_crit(wbs_stop_crit, method)
  cpn_max <- .reliever_validate_positive_integer(
    cpn_max, "cpn_max", allow_zero = TRUE
  )
  dm <- .reliever_validate_positive_integer(dm, "dm")
  if (is.null(init_cand)) {
    num_init <- .reliever_validate_positive_integer(num_init, "num_init")
    init_cand <- seq(0, 1, length.out = num_init + 2)
    init_cand <- init_cand[2:(num_init + 1)]
  } else {
    if (!is.numeric(init_cand) || length(init_cand) == 0L ||
        any(is.na(init_cand)) || any(init_cand <= 0 | init_cand >= 1)) {
      stop("init_cand must contain numbers in (0, 1).", call. = FALSE)
    }
    init_cand <- sort(unique(as.numeric(init_cand)))
  }
  if (method == "BS") M <- 0
  n <- .reliever_nobs(data)
  virtual_info <- .reliever_reg_fun_virtual_info(reg_fun, data, para_list)
  n_loss_outputs <- virtual_info$n_loss_outputs
  lr_m <- .reliever_search_intervals(method, n, dm, M, wbs_seed)

  individual_loss_fun <- .reliever_make_individual_loss_fun(
    data = data,
    reg_fun = reg_fun,
    dc_grid = NULL,
    para_list = para_list,
    n_loss_outputs = n_loss_outputs
  )
  if (echo) timer_start <- proc.time()[["elapsed"]]
  res_tem <- wbs_r_twostep_loss_outputs(
    n, cpn_max, dm, lr_m, as.integer(seq_len(n_loss_outputs)),
    individual_loss_fun, n_loss_outputs,
    as.numeric(init_cand), method == "WBS_recursive"
  )
  if (echo) {
    message(sprintf(
      "twostep-%s: %.3f sec elapsed",
      method, proc.time()[["elapsed"]] - timer_start
    ))
  }

  setup <- list(
    method = method,
    cpn_max = cpn_max,
    n = n,
    dm = dm,
    dc_grid = NULL,
    wbs_stop_crit = wbs_stop_crit
  )
  cpd_result <- .reliever_format_backend_cpd_result(
    backend_res = res_tem,
    setup = setup,
    pen_val = 1
  )
  cpd_result$cpd_path$candidates$loss[] <- NA_real_
  n_model_fit <- as.numeric(res_tem$n_model_fit)
  model_fit_time <- as.numeric(res_tem$model_fit_time)
  total_time <- as.numeric(res_tem$total_time)
  cpd_time <- total_time - model_fit_time
  out <- list(
    cpd_path = cpd_result$cpd_path,
    run_meta = .reliever_run_meta(
      virtual_info$loss_output_meta, seq_len(n_loss_outputs)
    ),
    n_model_fit = n_model_fit,
    model_fit_time = model_fit_time,
    cpd_time = cpd_time,
    total_time = total_time,
    diagnostics = if (detail) {
      list(
        path_score = cpd_result$path_score,
        gain_mat = res_tem$gain_mat,
        split_mat = res_tem$split_mat
      )
    },
    settings = list(
      method = method,
      n = n,
      dm = dm,
      cpn_max = cpn_max,
      M = if (method %in% c("WBS", "WBS_recursive")) M,
      wbs_seed = if (method %in% c("WBS", "WBS_recursive")) wbs_seed,
      wbs_stop_crit = wbs_stop_crit,
      init_cand = init_cand,
      cpn_crit = "none"
    ),
    loss_spec = list(reg_fun = reg_fun, args = para_list)
  )
  .reliever_finalize_result(out)
}


#' Locally refine a changepoint set
#'
#' Refine each supplied changepoint by minimizing left/right split loss in a
#' local window bounded by its neighboring initial changepoints. Refinement is
#' performed separately for every loss output returned by \code{reg_fun}.
#'
#' @param data Data object accepted by \code{reg_fun}.
#' @param tau_init Integer initial changepoint locations. The values are sorted
#'   and duplicates are removed before refinement; each value represents one
#'   changepoint to refine.
#' @param reg_fun Custom interval-loss function following the generic fitting
#'   contract.
#' @param dm Minimum segment length used in local split evaluation.
#' @param r Neighborhood shrinkage parameter in `[0, 1]`. Larger values search
#'   closer to each initial changepoint. Zero uses the full neighboring
#'   interval; one returns the initial locations unchanged.
#' @param echo Print timing messages.
#' @param ... Additional arguments passed to \code{reg_fun}.
#'
#' @return A list with \code{tau_est}, \code{run_meta}, and
#'   \code{total_time}. \code{tau_est} has one row per loss output returned by
#'   \code{reg_fun} and one column per initial changepoint. \code{run_meta}
#'   describes those rows using the same metadata contract as
#'   \code{reliever_generic()}. For a single-loss model, read the refined
#'   changepoints from \code{tau_est[1, ]}.
#' @seealso \code{\link{twostep}()}, \code{\link{reliever_generic}()},
#'   \code{\link{reg_fun_lasso_solpath}()},
#'   \code{\link{reg_fun_mean}()}
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
#' # Refine two initial changepoints.
#' refined <- local_refine(x, tau_init = c(295, 605), dm = 30,
#'                         reg_fun = reg_fun_mean)
#' refined$tau_est
local_refine <- function(data, tau_init, reg_fun, dm = 50, r = 0.1,
                         echo = FALSE, ...) {
  if (missing(reg_fun)) {
    stop("reg_fun must be supplied to local_refine().", call. = FALSE)
  }
  para_list <- list(...)
  n <- .reliever_nobs(data)
  dm <- .reliever_validate_positive_integer(dm, "dm")
  if (!is.numeric(r) || length(r) != 1L || is.na(r) || r < 0 || r > 1) {
    stop("r must be a single number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(tau_init) || length(tau_init) == 0L ||
      any(is.na(tau_init)) || any(tau_init != floor(tau_init)) ||
      any(tau_init <= 0 | tau_init >= n)) {
    stop("tau_init must contain changepoint locations in (0, n).", call. = FALSE)
  }
  tau_init <- sort(unique(as.integer(tau_init)))
  virtual_info <- .reliever_reg_fun_virtual_info(reg_fun, data, para_list)
  n_loss_outputs <- virtual_info$n_loss_outputs
  run_meta <- .reliever_run_meta(
    virtual_info$loss_output_meta, seq_len(n_loss_outputs)
  )
  if (r == 1) {
    return(list(
      tau_est = matrix(
        rep(tau_init, each = n_loss_outputs), nrow = n_loss_outputs
      ),
      run_meta = run_meta,
      total_time = 0
    ))
  }
  individual_loss_fun <- .reliever_make_individual_loss_fun(
    data = data,
    reg_fun = reg_fun,
    dc_grid = NULL,
    para_list = para_list,
    n_loss_outputs = n_loss_outputs
  )
  tau_init <- c(0, tau_init, n)

  refine_one <- function(i, j, tau) {
    if (j - i + 1 < 2 * dm) {
      stop("the interval [i, j] is too narrow", call. = FALSE)
    }

    loss_left <- matrix(
      individual_loss_fun(i, tau, i, j - dm),
      nrow = j - dm - i + 1
    )
    loss_left <- apply(loss_left, 2, cumsum)
    loss_right <- matrix(
      individual_loss_fun(tau + 1, j, i + dm, j),
      nrow = j - dm - i + 1
    )
    loss_right <- apply(
      matrix(loss_right[nrow(loss_right):1, ], nrow = nrow(loss_right)),
      2,
      cumsum
    )
    split_loss <- loss_left[dm:(j - i - dm + 1), , drop = FALSE] +
      loss_right[(j - i - dm + 1):dm, , drop = FALSE]
    as.vector(apply(split_loss, 2L, which.min)) + i + dm - 2L
  }

  start_time <- proc.time()[["elapsed"]]
  tau_est <- matrix(NA, n_loss_outputs, length(tau_init) - 2)
  for (i in 2:(length(tau_init) - 1L)) {
    tau_m <- tau_init[i]
    tau_l <- max(1L, floor(tau_init[i - 1] * (1 - r) + tau_m * r))
    tau_r <- ceiling(tau_init[i + 1] * (1 - r) + tau_m * r)
    tau_est[, i - 1L] <- refine_one(tau_l, tau_r, tau_m)
  }
  total_time <- proc.time()[["elapsed"]] - start_time
  if (echo) message(sprintf("refining: %.3f sec elapsed", total_time))
  list(
    tau_est = tau_est,
    run_meta = run_meta,
    total_time = total_time
  )
}
