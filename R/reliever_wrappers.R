# Shared wrapper helper --------------------------------------------------------

.reliever_with_reg_fun <- function(data, reg_fun, reg_args,
                                   cpn_max, dm, cov_rate, method, cpn_crit,
                                   pen_val, prune_value,
                                   M, wbs_seed, wbs_stop_crit,
                                   detail, cache_backend, owner_key, echo,
                                   dc_grid_size,
                                   cache_profile, run_cpd_ids,
                                   dc_grid) {
  reliever_args <- list(
    data = data,
    reg_fun = reg_fun,
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
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids
  )
  do.call(
    reliever_generic,
    c(reliever_args, reg_args, list(dc_grid = dc_grid))
  )
}

# Solution-path wrappers -------------------------------------------------------

#' Reliever with lasso solution-path losses
#'
#' Run Reliever with glmnet lasso losses over a common normalized lambda path.
#' A stored value \code{lam_set} is sample-size independent: an interval with
#' \eqn{m} rows passes \eqn{lam_set/\sqrt m} to \code{glmnet}. Use this
#' when the response is in the first column of \code{data} and the remaining
#' columns are predictors. For each lambda, the function computes candidate
#' segmentations with different numbers of changepoints.
#' The main built-in entry is
#' \code{reliever(X, y, cpd_family = "lasso")}. This wrapper exposes all lasso
#' arguments and accepts the equivalent response-first \code{data} matrix.
#' If \code{lam_set = NULL}, it asks a full-data glmnet fit for a path controlled
#' by \code{nlambda}, rescales it to Reliever's convention, and extends both ends
#' by two values. This constructs a candidate range; it does not select lambda.
#' Supply \code{lam_set} to use a fixed user-defined path instead.
#' Built-in lasso fits use \code{intercept = FALSE} and
#' \code{standardize = FALSE}; put predictors on their intended scale before
#' analysis.
#'
#' @section Selection:
#' Use \code{cv.reliever(..., cpd_family = "lasso")} to select lambda and K by
#' CPSS-style outer cross-validation. For Gaussian loss,
#' \code{select_by_run(fit, cpn_crit = "rss_sic")} selects K separately for
#' each fitted lambda path; binomial and Poisson losses use \code{"sic"} instead.
#'
#' @param cpn_max,dm,cov_rate,method,cpn_crit Common path controls; see
#'   \code{\link{reliever}()}.
#' @param pen_val,prune_value,M,wbs_seed,wbs_stop_crit Common penalty and
#'   search controls; see \code{\link{reliever}()}.
#' @param detail,cache_backend,owner_key,echo Common output and cache controls;
#'   see \code{\link{reliever}()}.
#' @param dc_grid_size,dc_grid Common candidate-grid controls; see
#'   \code{\link{reliever}()}.
#' @param cache_profile,run_cpd_ids Advanced reusable-cache and loss-output
#'   controls; see \code{\link{reliever_generic}()}.
#' @inheritParams reg_fun_lasso_solpath
#' @param data Numeric matrix with the response in the first column and
#'   at least two predictors in the remaining columns, as required by
#'   \code{glmnet}.
#' @param lam_set Optional finite positive normalized lambda values. The default
#'   \code{NULL} constructs one data-adaptive path and reuses it throughout the
#'   changepoint search. The low-level lasso-path help page defines the scale.
#' @param nlambda Integer scalar of at least 2 passed to the full-data
#'   \code{glmnet} fit that constructs the automatic path when
#'   \code{lam_set = NULL}. The default is 30. It is ignored when an explicit
#'   \code{lam_set} is supplied.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @seealso \code{\link{reliever}()}, \code{\link{cv.reliever}()}, and
#'   \code{\link{reg_fun_lasso_solpath}()}.
#' @export
#'
#' @examples
#' \donttest{
#' # This high-dimensional lasso example takes more than five seconds.
#' set.seed(2026)
#' n <- 900
#' p <- 100
#' tau <- c(300, 600)
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' data <- dgp_linear_regression(n, p, tau, b0, delta, sig = 1)$data
#'
#' # Fit the in-sample lambda paths; no K or model setting is selected.
#' fit_in <- reliever_lasso(
#'   data = data, cpn_max = 7, dm = 30, cov_rate = 0.8, method = "SN"
#' )
#' summary(fit_in)
#'
#' # For Gaussian loss, RSS-SIC can select K separately at every lambda.
#' selected <- select_by_run(result = fit_in, cpn_crit = "rss_sic")
#' selected
#' # At least one candidate lambda should recover the true K.
#' stopifnot(any(selected$K_hat == 2L))
#' }
reliever_lasso <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                           method = "SN", cpn_crit = "none",
                           pen_val = 1, prune_value = 0,
                           M = 100, wbs_seed = NULL,
                           wbs_stop_crit = NULL,
                           detail = FALSE,
                           cache_backend = "by_loss_block",
                           owner_key = TRUE,
                           echo = FALSE,
                           dc_grid_size = NULL,
                           cache_profile = NULL, run_cpd_ids = NULL,
                           dc_grid = NULL,
                           lam_set = NULL, nlambda = 30L,
                           family = "gaussian", thresh = 1e-7) {
  lam_set <- .reliever_resolve_lam_set(
    data, lam_set, family, thresh, nlambda = nlambda
  )
  out <- .reliever_with_reg_fun(
    data = data,
    reg_fun = reg_fun_lasso_solpath,
    reg_args = list(lam_set = lam_set, family = family, thresh = thresh),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
  .reliever_set_input_spec(
    out, .reliever_response_input_spec(data, y = NULL)
  )
}

# Cross-fitting reg_fun wrappers ----------------------------------------------

#' Reliever with cross-fitted lasso losses
#'
#' Run Reliever with lasso losses built by interval-level cross-fitting. The
#' result contains a cross-fitted path for every normalized lambda setting
#' held fixed across all intervals, plus an interval-adaptive recycled-CV
#' (\code{recv}) path in which that setting is chosen separately within every
#' fitted interval. A fold with \eqn{m} training rows passes
#' \eqn{lam_set/\sqrt m} to \code{glmnet}. The default
#' \code{cpn_crit = "none"} fits and retains these paths without choosing K or
#' one global normalized lambda setting. If \code{lam_set = NULL}, a
#' data-adaptive full-data glmnet path is converted to this normalized scale
#' once and reused for every interval.
#' Built-in lasso fits use \code{intercept = FALSE} and
#' \code{standardize = FALSE}; put predictors on their intended scale before
#' analysis.
#' The main built-in entry is
#' \code{reliever(X, y, cpd_family = "lasso_crossfit")}. This wrapper exposes
#' all crossfit-lasso arguments and accepts response-first \code{data}.
#'
#' @section Selection:
#' \itemize{
#'   \item Interval-adaptive ReCV uses \code{select_by_run()} with run type
#'   \code{"recv"} and criterion \code{"loss"}.
#'   \item Joint selection of one normalized lambda and K applies
#'   \code{select_across_runs()} to run type \code{"crossfit_homo_hyper"}.
#'   Request both crossfit output types when fitting.
#' }
#' Replace \code{"loss"} with \code{"sic"} in either selector to minimize
#' \eqn{ReCV\ loss + \log(n)K}. Inspect fixed-setting losses with a
#' hyperparameter plot at the chosen K.
#'
#' @inheritParams reliever_lasso
#' @inheritParams reg_fun_lasso_crossfit
#' @param data Numeric matrix with the response in the first column and
#'   at least two predictors in the remaining columns, as required by
#'   \code{glmnet}.
#' @param lam_set Optional finite positive normalized lambda values. The default
#'   \code{NULL} constructs one data-adaptive path and reuses it throughout the
#'   changepoint search. The low-level crossfit-lasso help page defines the scale.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @seealso \code{\link{reliever}()},
#'   \code{\link{reg_fun_lasso_crossfit}()}, and
#'   \code{\link{select_by_run}()}.
#' @export
#'
#' @examples
#' \donttest{
#' # This cross-fitted lasso example takes more than five seconds.
#' set.seed(2026)
#' n <- 450
#' p <- 20
#' tau <- c(150, 300)
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' data <- dgp_linear_regression(n, p, tau, b0, delta, sig = 1)$data
#' fit <- reliever_lasso_crossfit(
#'   data = data, cpn_max = 5, dm = 15, cov_rate = 0.8,
#'   nfolds = 2,
#'   loss_output_types = c("recv", "crossfit_homo_hyper"),
#'   method = "SN"
#' )
#' selected <- select_by_run(
#'   result = fit, run_type = "recv", cpn_crit = "loss"
#' )
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#'
#' # Homogeneous-hyperparameter ReCV: select one normalized setting and K.
#' selected_fixed <- select_across_runs(
#'   result = fit, run_type = "crossfit_homo_hyper", cpn_crit = "loss"
#' )
#' selected_fixed
#' stopifnot(identical(selected_fixed$K_hat, 2L))
#' plot(x = fit, run_type = "recv")
#' plot(
#'   x = fit, x_axis = "hyperparameter", K = 2,
#'   run_type = "crossfit_homo_hyper", cpn_crit = "loss", log = "x"
#' )
#' }
reliever_lasso_crossfit <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                                   method = "SN", cpn_crit = "none",
                                   pen_val = 1, prune_value = 0,
                                   M = 100, wbs_seed = NULL,
                                   wbs_stop_crit = NULL,
                                   detail = FALSE,
                                   cache_backend = "by_loss_block",
                                   owner_key = TRUE,
                                   echo = FALSE,
                                   dc_grid_size = NULL,
                                   cache_profile = NULL,
                                   run_cpd_ids = NULL,
                                   dc_grid = NULL,
                                   nfolds = 5,
                                   fold_type = "op",
                                   op_size = 1, buffer_lag = 0,
                                   fold_stable_const = 1,
                                   loss_output_types = "recv",
                                   lam_set = NULL, nlambda = 30L,
                                   family = "gaussian", thresh = 1e-7) {
  lam_set <- .reliever_resolve_lam_set(
    data, lam_set, family, thresh, nlambda = nlambda
  )
  out <- .reliever_with_reg_fun(
    data = data,
    reg_fun = reg_fun_lasso_crossfit,
    reg_args = list(
      nfolds = nfolds, lam_set = lam_set, family = family, thresh = thresh,
      fold_type = fold_type, op_size = op_size,
      buffer_lag = buffer_lag, fold_stable_const = fold_stable_const,
      loss_output_types = loss_output_types
    ),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
  .reliever_set_input_spec(
    out, .reliever_response_input_spec(data, y = NULL)
  )
}

#' Reliever with cross-fitted mean-square losses
#'
#' Run Reliever with mean-square losses built by interval-level recycled
#' cross-validation (ReCV). The main built-in entry is
#' \code{reliever(X, cpd_family = "mean_crossfit")}; this wrapper exposes all
#' mean-crossfit arguments.
#'
#' @section Selection:
#' Use \code{select_by_run(fit, run_type = "recv", cpn_crit = "loss")} for
#' unpenalized ReCV selection of K. Replace \code{"loss"} with \code{"sic"} in
#' that call to minimize \eqn{ReCV\ loss + \log(n)K}.
#'
#' @inheritParams reliever_lasso
#' @inheritParams reg_fun_mean_crossfit
#' @param data Numeric vector or matrix with observations in rows.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @seealso \code{\link{reliever}()}, \code{\link{reg_fun_mean_crossfit}()},
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
#' fit <- reliever_mean_crossfit(
#'   data = x, cpn_max = 7, dm = 30, cov_rate = 0.8,
#'   nfolds = 2, method = "SN"
#' )
#' selected <- select_by_run(
#'   result = fit, run_type = "recv", cpn_crit = "loss"
#' )
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#'
#' # A conservative alternative is:
#' selected_sic <- select_by_run(
#'   result = fit, run_type = "recv", cpn_crit = "sic"
#' )
#' selected_sic
#' stopifnot(identical(selected_sic$K_hat, 2L))
reliever_mean_crossfit <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                                  method = "SN", cpn_crit = "none",
                                  pen_val = 1, prune_value = 0,
                                  M = 100, wbs_seed = NULL,
                                  wbs_stop_crit = NULL,
                                  detail = FALSE,
                                  cache_backend = "by_loss_block",
                                  owner_key = TRUE,
                                  echo = FALSE,
                                  dc_grid_size = NULL,
                                  cache_profile = NULL,
                                  run_cpd_ids = NULL,
                                  dc_grid = NULL,
                                  nfolds = 5,
                                  fold_type = "op",
                                  op_size = 1, buffer_lag = 0,
                                  fold_stable_const = 1,
                                  loss_output_types = "recv") {
  .reliever_with_reg_fun(
    data = data,
    reg_fun = reg_fun_mean_crossfit,
    reg_args = list(
      nfolds = nfolds, fold_type = fold_type, op_size = op_size,
      buffer_lag = buffer_lag, fold_stable_const = fold_stable_const,
      loss_output_types = loss_output_types
    ),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
}

# Kernel and nonparametric wrappers -------------------------------------------

#' Reliever with KDE negative-log-likelihood solution-path losses
#'
#' Detect distribution changes by fitting one isotropic density kernel over a
#' path of bandwidths. Pass the original observations in rows; the wrapper
#' computes one reusable distance matrix before the changepoint search. A
#' precomputed squared-distance matrix is accepted through
#' \code{input_type = "distance"} and \code{var_dim}. The function returns one
#' candidate path per bandwidth.
#' The main built-in entry is
#' \code{reliever(X, cpd_family = "kde_nll")}; this wrapper exposes all
#' KDE-NLL arguments.
#'
#' @section Selection:
#' Use \code{reliever(X, cpd_family = "kde_nll_crossfit")} for data-driven
#' bandwidth selection by cross-fitted loss. With an independent evaluation
#' sample, \code{select_holdout()} can instead select bandwidth and K.
#' On an ordinary KDE-NLL fit,
#' \code{select_by_run(fit, cpn_crit = "sic")} selects K separately for every
#' fixed bandwidth; it does not select a bandwidth.
#'
#' @inheritParams reliever_lasso
#' @param bandwidth_vec Positive KDE bandwidths \eqn{h}. When \code{NULL}, a
#'   scale-adaptive grid is constructed from 20 log-spaced values \eqn{b}
#'   satisfying
#'   \eqn{b/\sqrt{s_D}\in[1/\sqrt{120},1]}, where \eqn{s_D} is the median squared
#'   distance among observation pairs separated by at most
#'   \eqn{\lceil\sqrt n\rceil} time points. Using local pairs estimates the
#'   within-segment scale without being dominated by large changes.
#'   This reduces to approximately \eqn{[\sqrt{p/60},\sqrt{2p}]} for standardized
#'   \eqn{p}-dimensional Gaussian data. Kernel-specific conversions are:
#'   \itemize{
#'     \item Gaussian: \eqn{h=b}.
#'     \item Radial Laplace: \eqn{h=b/\sqrt{p+1}}.
#'     \item Student with \eqn{\nu>2}: \eqn{h=b\sqrt{(\nu-2)/\nu}}.
#'   }
#'   These give every kernel the same marginal variance at a fixed \eqn{b}.
#'   Student kernels with \code{df <= 2} require an explicit bandwidth vector.
#'   An explicit vector is always used unchanged.
#' @param data Numeric observations in rows, or an \eqn{n} by \eqn{n}
#'   precomputed squared-distance matrix.
#' @param input_type Input representation. \code{"data"} treats rows as raw
#'   observations and infers \code{var_dim}; \code{"distance"} requires a
#'   symmetric squared-distance matrix and explicit \code{var_dim}.
#'   \code{"auto"}, the default, preserves the former distance-matrix call
#'   when a distance-like matrix and \code{var_dim} are supplied. A
#'   distance-like matrix without \code{var_dim} is rejected as ambiguous;
#'   other inputs are treated as raw observations. For an unusually
#'   distance-like raw square matrix, set \code{input_type = "data"}
#'   explicitly.
#' @param var_dim Dimension of the original observations. It is inferred for
#'   raw data and required for a precomputed squared-distance matrix.
#' @param kernel Isotropic density kernel: \code{"gaussian"},
#'   \code{"laplace"}, or \code{"student"}. Here Laplace means the radial
#'   \eqn{L_2} density, not a coordinate-wise product or \eqn{L_1} kernel.
#'   Gaussian uses
#'   \eqn{\exp\{-\|x-y\|^2/(2h^2)\}}; radial Laplace uses
#'   \eqn{\exp\{-\|x-y\|/h\}}; and Student uses the multivariate Student
#'   density with scale matrix \eqn{h^2 I}. Other fixed-kernel choices are
#'   documented under \code{\link{reliever_kde_l2}()}.
#' @param kernel_args Named list of kernel-specific settings. Student accepts
#'   \code{df}, a positive degrees-of-freedom value with default 5. Gaussian
#'   and Laplace currently accept no additional settings.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @seealso \code{\link{reliever}()},
#'   \code{\link{reg_fun_kde_nll_solpath}()}, and
#'   \code{\link{select_by_run}()}.
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 200
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' fit <- reliever_kde_nll(
#'   data = x, cpn_max = 5, dm = 20, cov_rate = 0.6, # just for speed in example
#'   method = "SN", bandwidth_vec = c(0.5, 1, 2)
#' )
#' # Optional SIC selects K separately at every fixed bandwidth.
#' selected <- select_by_run(result = fit, cpn_crit = "sic")
#' selected
#' # At least one candidate bandwidth should recover the true K.
#' stopifnot(any(selected$K_hat == 2L))
reliever_kde_nll <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                             method = "SN", cpn_crit = "none",
                             pen_val = 1, prune_value = 0,
                             M = 100, wbs_seed = NULL,
                             wbs_stop_crit = NULL,
                             detail = FALSE,
                             cache_backend = "by_loss_block",
                             owner_key = TRUE,
                             echo = FALSE,
                             dc_grid_size = NULL,
                             cache_profile = NULL, run_cpd_ids = NULL,
                             dc_grid = NULL,
                             bandwidth_vec = NULL, var_dim = NULL,
                             kernel = "gaussian",
                             kernel_args = list(),
                             input_type = c("auto", "data", "distance")) {
  prepared <- .kde_prepare_fit_input(
    data = data, var_dim = var_dim, input_type = input_type,
    kernel = kernel, kernel_args = kernel_args,
    bandwidth_vec = bandwidth_vec
  )
  out <- .reliever_with_reg_fun(
    data = prepared$data,
    reg_fun = reg_fun_kde_nll_solpath,
    reg_args = list(
      bandwidth_vec = prepared$bandwidth_vec,
      var_dim = prepared$var_dim,
      kernel = prepared$kernel,
      kernel_args = prepared$kernel_args,
      distance_power = prepared$distance_power
    ),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
  .reliever_set_input_spec(out, prepared$input_spec)
}

#' Reliever with cross-fitted KDE negative-log-likelihood losses
#'
#' Detect distribution changes with isotropic KDE losses cross-fitted
#' separately within each interval. Pass the original observations in rows;
#' one bandwidth-independent distance matrix is computed before searching.
#' Alternatively, supply a precomputed squared-distance matrix. Kernel
#' bandwidth is chosen separately within each fitted interval for the
#' \code{recv} loss path. The default retains that path and the fixed-bandwidth
#' cross-fitted paths without choosing K or one global bandwidth.
#' The main built-in entry is
#' \code{reliever(X, cpd_family = "kde_nll_crossfit")}; this wrapper exposes
#' all crossfit KDE-NLL arguments.
#'
#' @section Selection:
#' \itemize{
#'   \item Interval-adaptive ReCV uses \code{select_by_run()} with run type
#'   \code{"recv"} and criterion \code{"loss"}.
#'   \item Homogeneous ReCV jointly selects one bandwidth and K from stored
#'   fixed-bandwidth rows.
#' }
#' Apply \code{select_across_runs()} for the homogeneous workflow.
#'
#' Its run type is \code{"crossfit_homo_hyper"}; request both crossfit outputs
#' when fitting.
#' Replace \code{"loss"} with \code{"sic"} in either selector to minimize
#' \eqn{ReCV\ loss + \log(n)K}.
#'
#' To check the automatic bandwidth range at a chosen K, use a hyperparameter
#' plot with run type \code{"crossfit_homo_hyper"}.
#'
#' @inheritParams reliever_lasso
#' @inheritParams reliever_kde_nll
#' @inheritParams reg_fun_kde_nll_crossfit
#' @param data Numeric observations in rows, or an \eqn{n} by \eqn{n}
#'   precomputed squared-distance matrix. See \code{input_type}.
#' @param var_dim Dimension of the original observations. It is inferred for
#'   raw data and required for a precomputed squared-distance matrix.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @seealso \code{\link{reliever}()},
#'   \code{\link{reg_fun_kde_nll_crossfit}()}, and
#'   \code{\link{select_by_run}()}.
#' @export
#'
#' @examples
#' \donttest{
#' # This cross-fitted KDE example takes more than five seconds.
#' set.seed(2026)
#' n_seg <- 150
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' fit <- reliever_kde_nll_crossfit(
#'   data = x, cpn_max = 5, dm = 15, cov_rate = 0.6, # for speed in example
#'   kernel = "student", kernel_args = list(df = 5),
#'   nfolds = 2,
#'   loss_output_types = c("recv", "crossfit_homo_hyper"),
#'   method = "SN"
#' )
#' selected <- select_by_run(
#'   result = fit, run_type = "recv", cpn_crit = "sic"
#' )
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#' selected_fixed <- select_across_runs(
#'   result = fit, run_type = "crossfit_homo_hyper", cpn_crit = "sic"
#' )
#' selected_fixed
#' stopifnot(identical(selected_fixed$K_hat, 2L))
#' plot(fit, run_type = "recv")
#' plot(
#'   fit, x_axis = "hyperparameter", K = selected$K_hat[[1L]],
#'   run_type = "crossfit_homo_hyper", cpn_crit = "sic", log = "x"
#' )
#' }
reliever_kde_nll_crossfit <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                                      method = "SN", cpn_crit = "none",
                                      pen_val = 1, prune_value = 0,
                                      M = 100, wbs_seed = NULL,
                                      wbs_stop_crit = NULL,
                                      detail = FALSE,
                                      cache_backend = "by_loss_block",
                                      owner_key = TRUE,
                                      echo = FALSE,
                                      dc_grid_size = NULL,
                                      cache_profile = NULL,
                                      run_cpd_ids = NULL,
                                      dc_grid = NULL,
                                      nfolds = 5,
                                      fold_type = "op",
                                      op_size = 1, buffer_lag = 0,
                                      fold_stable_const = 1,
                                      loss_output_types = "recv",
                                      bandwidth_vec = NULL, var_dim = NULL,
                                      kernel = "gaussian",
                                      kernel_args = list(),
                                      input_type = c("auto", "data", "distance")) {
  prepared <- .kde_prepare_fit_input(
    data = data, var_dim = var_dim, input_type = input_type,
    kernel = kernel, kernel_args = kernel_args,
    bandwidth_vec = bandwidth_vec
  )
  out <- .reliever_with_reg_fun(
    data = prepared$data,
    reg_fun = reg_fun_kde_nll_crossfit,
    reg_args = list(
      bandwidth_vec = prepared$bandwidth_vec,
      var_dim = prepared$var_dim,
      kernel = prepared$kernel,
      kernel_args = prepared$kernel_args,
      distance_power = prepared$distance_power,
      nfolds = nfolds,
      fold_type = fold_type, op_size = op_size,
      buffer_lag = buffer_lag, fold_stable_const = fold_stable_const,
      loss_output_types = loss_output_types
    ),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
  .reliever_set_input_spec(out, prepared$input_spec)
}

#' Reliever with fixed-kernel KDE L2 loss
#'
#' Detect distributional changes with an empirical-grid kernel-density
#' \eqn{L_2} loss. For every candidate segment,
#' \code{\link{reg_fun_kde_l2}()} estimates the segment KDE and evaluates the
#' corresponding squared \eqn{L_2} discrepancy on a fixed common grid. With
#' a fixed density \code{kernel}, this wrapper computes one kernel-evaluation matrix
#' and reuses it throughout the search; that matrix is a computational
#' representation of the KDE loss, not the statistical target itself.
#'
#' With \code{kernel = NULL}, \code{data} is an already computed
#' kernel-evaluation matrix; a generic fixed feature matrix is also accepted as
#' a fixed-feature extension. The kernel and bandwidth remain fixed throughout
#' the search.
#' The main built-in entry is
#' \code{reliever(X, cpd_family = "kde_l2", kernel = ..., bandwidth = ...)};
#' this wrapper exposes all fixed-kernel KDE-L2 arguments.
#'
#' @section Selection:
#' Use \code{select_by_run(fit, cpn_crit = "rss_sic")} to select K by
#' minimizing \eqn{n\log(RSS/n)/2 + \log(n)K}. The kernel and bandwidth remain
#' fixed.
#'
#' @inheritParams reliever_lasso
#' @param data With \code{kernel = NULL}, a numeric kernel-evaluation matrix
#'   with observations in rows and common evaluation locations in columns; an
#'   \eqn{n} by \eqn{n} Gram matrix is the ordinary case. A generic fixed
#'   feature matrix is also accepted. For the KDE-L2 interpretation, its
#'   columns are fixed kernel-density evaluations on a common grid. Otherwise,
#'   \code{data} contains the raw numeric observations used to construct one
#'   kernel-evaluation matrix.
#' @param kernel Fixed kernel. The default \code{NULL} means \code{data}
#'   is already a feature or Gram matrix. Built-in choices are:
#'   \itemize{
#'     \item \code{"gaussian"} (or \code{"rbf"}) and \code{"laplace"};
#'     \item \code{"linear"} and \code{"polynomial"};
#'     \item \code{"matern32"} and \code{"matern52"}; and
#'     \item \code{"rational_quadratic"}.
#'   }
#'   Gaussian, radial
#'   Laplace, Matérn-3/2, and Matérn-5/2 admit a density-kernel interpretation
#'   up to an omitted fixed normalizing constant. Rational quadratic does so
#'   in dimension \eqn{p} when \code{alpha > p / 2}. Linear and polynomial
#'   kernels are fixed-feature extensions. A custom function must accept
#'   \code{x}, \code{y}, and the
#'   named values in \code{kernel_args}, plus \code{bandwidth} when it is
#'   supplied, and return a finite
#'   \code{nrow(x)} by \code{nrow(y)} matrix. Its full-data matrix must be
#'   symmetric.
#' @param bandwidth One fixed positive bandwidth for a radial named kernel, or
#'   an optional value passed to a custom kernel. Use \code{NULL} for linear
#'   and polynomial kernels.
#' @param kernel_args Named list of additional fixed settings. Polynomial
#'   accepts \code{degree = 2}, \code{scale = 1}, and \code{offset = 1}.
#'   Rational quadratic accepts \code{alpha = 1}. For a custom function, all
#'   entries are passed unchanged.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @encoding UTF-8
#' @references Padilla, O. H. M., Yu, Y., Wang, D., and Rinaldo, A. (2021).
#'   Optimal nonparametric multivariate change point detection and localization.
#'   \emph{IEEE Transactions on Information Theory}, 68(3), 1922--1944.
#' @seealso \code{\link{reliever}()}, \code{\link{reg_fun_kde_l2}()},
#'   \code{\link{reliever_kde_nll_crossfit}()},
#'   \code{\link{select_by_run}()}
#' @export
#'
#' @examples
#' \donttest{
#' # Constructing and searching the KDE-L2 Gram matrix takes over five seconds.
#' set.seed(2026)
#' n_seg <- 300
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' bandwidth <- sqrt(10)
#' fit <- reliever_kde_l2(
#'   data = x, cpn_max = 7, dm = 30, cov_rate = 0.8, method = "SN",
#'   kernel = "gaussian", bandwidth = bandwidth
#' )
#' selected <- select_by_run(result = fit, cpn_crit = "rss_sic")
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#' }
reliever_kde_l2 <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                            method = "SN", cpn_crit = "none",
                            pen_val = 1, prune_value = 0,
                            M = 100, wbs_seed = NULL,
                            wbs_stop_crit = NULL,
                            detail = FALSE,
                            cache_backend = "by_loss_block",
                            owner_key = TRUE,
                            echo = FALSE,
                            dc_grid_size = NULL,
                            cache_profile = NULL, run_cpd_ids = NULL,
                            dc_grid = NULL,
                            kernel = NULL,
                            bandwidth = NULL,
                            kernel_args = list()) {
  prepared <- .kernel_l2_prepare_input(
    data = data, kernel = kernel, bandwidth = bandwidth,
    kernel_args = kernel_args
  )
  out <- .reliever_with_reg_fun(
    data = prepared$data,
    reg_fun = reg_fun_kde_l2,
    reg_args = list(),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
  .reliever_set_input_spec(out, prepared$input_spec)
}

#' Reliever with univariate nonparametric CDF loss
#'
#' Detect arbitrary changes in a univariate distribution by aggregating
#' empirical-CDF log losses. Unlike mean-change detection, this can respond to
#' changes in location, scale, or distributional shape.
#' The main built-in entry is \code{reliever(x, cpd_family = "nmcd")}; this
#' wrapper exposes all empirical-CDF arguments.
#'
#' @section Selection:
#' Use \code{select_by_run(fit, cpn_crit = "sic")} to select K by minimizing
#' \eqn{CDF\ loss + \log(n)K}. A non-negative numeric criterion supplies
#' another additive penalty per changepoint.
#'
#' @inheritParams reliever_lasso
#' @inheritParams reg_fun_nmcd
#' @param data Numeric vector or one-column matrix.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @references Zou, C., Yin, G., Feng, L., and Wang, Z. (2014).
#'   Nonparametric maximum likelihood approach to multiple change-point
#'   problems. \emph{The Annals of Statistics}, 42(3), 970--1002.
#'   DOI: 10.1214/14-AOS1210.
#' @seealso \code{\link{reliever}()}, \code{\link{reg_fun_nmcd}()},
#'   \code{\link{select_by_run}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- c(
#'   rnorm(n_seg, mean = 0, sd = 0.2),
#'   rnorm(n_seg, mean = 4, sd = 0.2),
#'   rnorm(n_seg, mean = -4, sd = 0.2)
#' )
#' fit <- reliever_nmcd(
#'   data = x, cpn_max = 7, dm = 30, cov_rate = 0.8, method = "SN"
#' )
#' selected <- select_by_run(result = fit, cpn_crit = "sic")
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
reliever_nmcd <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                          method = "SN", cpn_crit = "none",
                          pen_val = 1, prune_value = 0,
                          M = 100, wbs_seed = NULL,
                          wbs_stop_crit = NULL,
                          detail = FALSE,
                          cache_backend = "by_loss_block",
                          owner_key = TRUE,
                          echo = FALSE,
                          dc_grid_size = NULL,
                          cache_profile = NULL, run_cpd_ids = NULL,
                          dc_grid = NULL,
                          w_trunc = 0, sort_X = NULL) {
  if (is.null(sort_X)) {
    sort_X <- sort(as.numeric(data))
  }
  .reliever_with_reg_fun(
    data = data,
    reg_fun = reg_fun_nmcd,
    reg_args = list(w_trunc = w_trunc, sort_X = sort_X),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
}

# Classifier crossfit wrappers ------------------------------------------------

#' Reliever with cross-fitted ranger-classifier losses
#'
#' Detect distribution changes by training ranger classifiers to distinguish
#' observations inside a candidate interval from those outside it. Predicted
#' class probabilities are converted to interval-level ReCV losses. This
#' requires the optional \code{ranger} package.
#' The main built-in entry is \code{reliever()} with
#' \code{cpd_family = "ranger_crossfit"}; this wrapper exposes all
#' ranger-specific arguments.
#'
#' @section Selection:
#' \itemize{
#'   \item Interval-adaptive ReCV uses \code{select_by_run()} with run type
#'   \code{"recv"} and criterion \code{"loss"}.
#'   \item Joint selection of one classifier setting and K applies
#'   \code{select_across_runs()} to run type \code{"crossfit_homo_hyper"}.
#' }
#' Request those paths through \code{loss_output_types} when fitting.
#' Replace \code{"loss"} with \code{"sic"} in either call to add a
#' \eqn{\log(n)K} penalty.
#'
#' @inheritParams reliever_lasso
#' @inheritParams reg_fun_ranger_crossfit
#' @param data Numeric matrix or data frame with observations in rows and
#'   classifier features in columns.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @encoding UTF-8
#' @references Londschien, M., Bühlmann, P., and Kovács, S. (2023). Random
#'   forests for change point detection. \emph{Journal of Machine Learning
#'   Research}, 24(216), 1--45.
#' @seealso \code{\link{reliever}()},
#'   \code{\link{reg_fun_ranger_crossfit}()}, and
#'   \code{\link{select_by_run}()}.
#' @export
#'
#' @examples
#' \donttest{
#' # This optional-model example takes more than five seconds.
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   set.seed(2026)
#'   n_seg <- 300
#'   x <- rbind(
#'     matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#'   )
#'   fit <- suppressWarnings(reliever_ranger_crossfit(
#'     data = x, cpn_max = 7, dm = 30, cov_rate = 0.8, nfolds = 2,
#'     hyper_set = data.frame(num.trees = 10, min.node.size = 5),
#'     ranger_args = list(seed = 2026), method = "SN"
#'   ))
#'   selected <- select_by_run(
#'     result = fit, run_type = "recv", cpn_crit = "loss"
#'   )
#'   selected
#'   stopifnot(identical(selected$K_hat, 2L))
#' }
#' }
reliever_ranger_crossfit <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                                    method = "SN", cpn_crit = "none",
                                    pen_val = 1, prune_value = 0,
                                    M = 100, wbs_seed = NULL,
                                    wbs_stop_crit = NULL,
                                    detail = FALSE,
                                    cache_backend = "by_loss_block",
                                    owner_key = TRUE,
                                    echo = FALSE,
                                    dc_grid_size = NULL,
                                    cache_profile = NULL,
                                    run_cpd_ids = NULL,
                                    dc_grid = NULL,
                                    nfolds = 5,
                                    fold_type = "op",
                                    op_size = 1, buffer_lag = 0,
                                    fold_stable_const = 1,
                                    loss_output_types = "recv",
                                    hyper_set = list(list()),
                                    ranger_args = list()) {
  .reliever_with_reg_fun(
    data = data,
    reg_fun = reg_fun_ranger_crossfit,
    reg_args = list(
      nfolds = nfolds, hyper_set = hyper_set,
      fold_type = fold_type, op_size = op_size,
      buffer_lag = buffer_lag, fold_stable_const = fold_stable_const,
      loss_output_types = loss_output_types,
      ranger_args = ranger_args
    ),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
}

#' Reliever with cross-fitted MLP-classifier losses
#'
#' Detect distribution changes by training a small multilayer perceptron (MLP)
#' to distinguish observations inside a candidate interval from those outside
#' it. Predicted class probabilities are converted to interval-level ReCV
#' losses. This uses the optional \code{nnet} package.
#' The main built-in entry is
#' \code{reliever(X, cpd_family = "mlp_crossfit")}; this wrapper exposes all
#' MLP-specific arguments.
#'
#' @section Selection:
#' \itemize{
#'   \item Interval-adaptive ReCV uses \code{select_by_run()} with run type
#'   \code{"recv"} and criterion \code{"loss"}.
#'   \item Joint selection of one classifier setting and K applies
#'   \code{select_across_runs()} to run type \code{"crossfit_homo_hyper"}.
#' }
#' Request those paths through \code{loss_output_types} when fitting.
#' Replace \code{"loss"} with \code{"sic"} in either call to add a
#' \eqn{\log(n)K} penalty.
#'
#' @inheritParams reliever_lasso
#' @inheritParams reg_fun_mlp_crossfit
#' @param data Numeric matrix with observations in rows and classifier features
#'   in columns.
#' @return A \code{reliever_result} object with the common structure documented
#'   in \code{\link{reliever}()}.
#' @seealso \code{\link{reliever}()},
#'   \code{\link{reg_fun_mlp_crossfit}()}, and
#'   \code{\link{select_by_run}()}.
#' @export
#'
#' @examples
#' \donttest{
#' # This optional-model example takes more than five seconds.
#' if (requireNamespace("nnet", quietly = TRUE)) {
#'   set.seed(2026)
#'   n_seg <- 300
#'   x <- rbind(
#'     matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#'   )
#'   set.seed(2026)
#'   fit <- suppressWarnings(reliever_mlp_crossfit(
#'     data = x, cpn_max = 7, dm = 30, cov_rate = 0.8, nfolds = 2,
#'     hyper_set = data.frame(size = 2, maxit = 30),
#'     method = "SN"
#'   ))
#'   # Numeric criteria supply an additive penalty per changepoint.
#'   selected <- select_by_run(
#'     result = fit, run_type = "recv", cpn_crit = 20
#'   )
#'   selected
#'   stopifnot(identical(selected$K_hat, 2L))
#' }
#' }
reliever_mlp_crossfit <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
                                 method = "SN", cpn_crit = "none",
                                 pen_val = 1, prune_value = 0,
                                 M = 100, wbs_seed = NULL,
                                 wbs_stop_crit = NULL,
                                 detail = FALSE,
                                 cache_backend = "by_loss_block",
                                 owner_key = TRUE,
                                 echo = FALSE,
                                 dc_grid_size = NULL,
                                 cache_profile = NULL,
                                 run_cpd_ids = NULL,
                                 dc_grid = NULL,
                                 nfolds = 5,
                                 fold_type = "op",
                                 op_size = 1, buffer_lag = 0,
                                 fold_stable_const = 1,
                                 loss_output_types = "recv",
                                 hyper_set = data.frame(size = 8),
                                 nnet_args = list()) {
  .reliever_with_reg_fun(
    data = data,
    reg_fun = reg_fun_mlp_crossfit,
    reg_args = list(
      nfolds = nfolds, hyper_set = hyper_set,
      fold_type = fold_type, op_size = op_size,
      buffer_lag = buffer_lag, fold_stable_const = fold_stable_const,
      loss_output_types = loss_output_types,
      nnet_args = nnet_args
    ),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value, M = M,
    wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail,
    cache_backend = cache_backend, owner_key = owner_key, echo = echo,
    dc_grid_size = dc_grid_size,
    cache_profile = cache_profile,
    run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
}
