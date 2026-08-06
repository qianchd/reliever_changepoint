.reliever_parametric_input_spec <- function(data, response_ncol) {
  data <- .reliever_parametric_matrix(data)
  response_ncol <- .reliever_validate_positive_integer(
    response_ncol, "response_ncol"
  )
  if (response_ncol >= ncol(data)) {
    stop("Regression data must contain at least one predictor column.",
         call. = FALSE)
  }
  list(
    type = "response_predictor",
    n_predictors = ncol(data) - response_ncol,
    n_response_columns = response_ncol,
    original_form = "response_first"
  )
}

#' Reliever for covariance changes
#'
#' Fit a single changepoint path using multivariate Gaussian negative twice
#' log-likelihood, with one common mean and interval-specific covariance
#' matrices. If \code{mu = NULL}, the common mean is estimated once from the
#' complete input and then held fixed in every interval.
#'
#' The main built-in entry is \code{reliever()} with
#' \code{cpd_family = "var"}. Use \code{cv.reliever()} with the same family
#' to select K by outer sample-splitting.
#'
#' @inheritParams reliever_lasso
#' @param data Numeric vector or matrix with observations in rows.
#' @param mu Optional fixed common mean. A scalar is applied to every column.
#' @return A \code{reliever_result} object.
#' @seealso \code{\link{reg_fun_var}()}, \code{\link{reliever_meanvar}()},
#'   \code{\link{cv.reliever}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' x <- c(rnorm(40, sd = 0.5), rnorm(40, sd = 2), rnorm(40, sd = 0.5))
#' fit <- reliever_var(
#'   data = x, cpn_max = 3, dm = 10, cov_rate = 0.6, method = "BS"
#' )
#' selected <- select_by_run(fit, cpn_crit = "sic")
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
reliever_var <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
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
                         mu = NULL) {
  data <- .reliever_parametric_matrix(data)
  if (is.null(mu)) {
    mu <- .reliever_stable_col_means(data)
  }
  .reliever_with_reg_fun(
    data = data, reg_fun = reg_fun_var, reg_args = list(mu = mu),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value,
    M = M, wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail, cache_backend = cache_backend, owner_key = owner_key,
    echo = echo, dc_grid_size = dc_grid_size,
    cache_profile = cache_profile, run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
}

#' Reliever for simultaneous mean and covariance changes
#'
#' Fit a single changepoint path using multivariate Gaussian negative twice
#' log-likelihood, estimating both the mean and maximum-likelihood covariance
#' separately in every fitted interval.
#'
#' The main built-in entry is \code{reliever()} with
#' \code{cpd_family = "meanvar"}. Use \code{cv.reliever()} with the same family
#' to select K by outer sample-splitting.
#'
#' @inheritParams reliever_var
#' @return A \code{reliever_result} object.
#' @seealso \code{\link{reg_fun_meanvar}()}, \code{\link{reliever_var}()},
#'   \code{\link{cv.reliever}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' x <- c(
#'   rnorm(40, mean = 0, sd = 0.5),
#'   rnorm(40, mean = 2, sd = 1.5),
#'   rnorm(40, mean = -1, sd = 0.7)
#' )
#' fit <- reliever_meanvar(
#'   data = x, cpn_max = 3, dm = 10, cov_rate = 0.6, method = "BS"
#' )
#' # A numeric criterion is an additive penalty per changepoint.
#' selected <- select_by_run(fit, cpn_crit = 20)
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
reliever_meanvar <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
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
                             dc_grid = NULL) {
  .reliever_with_reg_fun(
    data = data, reg_fun = reg_fun_meanvar, reg_args = list(),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value,
    M = M, wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail, cache_backend = cache_backend, owner_key = owner_key,
    echo = echo, dc_grid_size = dc_grid_size,
    cache_profile = cache_profile, run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
}

#' Reliever for linear-model changes
#'
#' Fit one ordinary least-squares model in each interval and use
#' observation-level squared residuals as the segment loss. The first column
#' of \code{data} is the response and the remaining columns are predictors.
#'
#' The main built-in entry is \code{reliever()} with
#' \code{cpd_family = "lm"}. Use \code{cv.reliever()} with the same inputs to
#' select K by outer sample-splitting.
#'
#' @inheritParams reliever_var
#' @param data Numeric matrix with the response first and predictors afterward.
#' @param intercept Add an intercept to each interval model.
#' @return A \code{reliever_result} object.
#' @seealso \code{\link{reg_fun_lm}()}, \code{\link{reliever_glm}()},
#'   \code{\link{cv.reliever}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n <- 120
#' x <- rnorm(n)
#' beta <- rep(c(1, -1, 0.5), each = 40)
#' y <- beta * x + rnorm(n, sd = 0.4)
#' fit <- reliever_lm(
#'   data = cbind(y, x), cpn_max = 3, dm = 10,
#'   cov_rate = 0.6, method = "BS"
#' )
#' selected <- select_by_run(fit, cpn_crit = "sic")
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
reliever_lm <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
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
                        intercept = TRUE) {
  data <- .reliever_parametric_matrix(data)
  input_spec <- .reliever_parametric_input_spec(data, 1L)
  out <- .reliever_with_reg_fun(
    data = data, reg_fun = reg_fun_lm,
    reg_args = list(intercept = intercept),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value,
    M = M, wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail, cache_backend = cache_backend, owner_key = owner_key,
    echo = echo, dc_grid_size = dc_grid_size,
    cache_profile = cache_profile, run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
  .reliever_set_input_spec(out, input_spec)
}

#' Reliever for generalized-linear-model changes
#'
#' Fit one generalized linear model in each interval and use the selected
#' family's observation-level deviance residuals as segment losses. The leading
#' \code{response_ncol} columns of \code{data} contain the response. Use two
#' response columns for binomial successes and failures.
#'
#' The main built-in entry is \code{reliever()} with
#' \code{cpd_family = "glm"}. Use \code{cv.reliever()} with the same inputs to select K by outer
#' sample-splitting.
#'
#' @inheritParams reliever_var
#' @param data Numeric matrix with the response column or columns first and
#'   predictors afterward.
#' @param family A GLM family object, family function, or the name of a
#'   \code{stats} family function.
#' @param intercept Add an intercept to each interval model.
#' @param response_ncol Number of leading response columns; use 2 for binomial
#'   successes and failures.
#' @return A \code{reliever_result} object.
#' @seealso \code{\link{reg_fun_glm}()}, \code{\link{reliever_lm}()},
#'   \code{\link{cv.reliever}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n <- 120
#' x <- rnorm(n)
#' beta <- rep(c(1.5, -1.5, 0.5), each = 40)
#' trials <- rep(5L, n)
#' success <- rbinom(n, trials, plogis(beta * x))
#' fit <- reliever_glm(
#'   data = cbind(success, trials - success, x),
#'   family = binomial(), response_ncol = 2,
#'   cpn_max = 3, dm = 10, cov_rate = 0.6, method = "BS"
#' )
#' # A numeric criterion is an additive penalty per changepoint.
#' selected <- select_by_run(fit, cpn_crit = 20)
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
reliever_glm <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
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
                         family, intercept = TRUE, response_ncol = 1L) {
  if (missing(family)) {
    stop("family is required for reliever_glm().", call. = FALSE)
  }
  data <- .reliever_parametric_matrix(data)
  input_spec <- .reliever_parametric_input_spec(data, response_ncol)
  out <- .reliever_with_reg_fun(
    data = data, reg_fun = reg_fun_glm,
    reg_args = list(
      family = family, intercept = intercept, response_ncol = response_ncol
    ),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value,
    M = M, wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail, cache_backend = cache_backend, owner_key = owner_key,
    echo = echo, dc_grid_size = dc_grid_size,
    cache_profile = cache_profile, run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
  .reliever_set_input_spec(out, input_spec)
}

#' Reliever for exponential-family distribution changes
#'
#' Fit one distribution in each interval and use observation-level negative
#' twice log-likelihood as the segment loss. Supported families are
#' \code{"binom"}, \code{"multinom"}, \code{"pois"}, \code{"exp"},
#' \code{"geom"}, \code{"diri"}, \code{"gamma"}, \code{"beta"},
#' \code{"chisq"}, and \code{"invgauss"}.
#'
#' The main built-in entry is \code{reliever()} with
#' \code{cpd_family = "em"}. Use \code{cv.reliever()} with the same inputs to select K by outer
#' sample-splitting.
#'
#' @inheritParams reliever_var
#' @param family Exponential-family distribution name.
#' @param size Number of trials, required for binomial and multinomial data.
#' @return A \code{reliever_result} object.
#' @seealso \code{\link{reg_fun_em}()}, \code{\link{cv.reliever}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' x <- c(rexp(40, rate = 1), rexp(40, rate = 4), rexp(40, rate = 1))
#' fit <- reliever_em(
#'   data = x, family = "exp", cpn_max = 3, dm = 10,
#'   cov_rate = 0.6, method = "BS"
#' )
#' selected <- select_by_run(fit, cpn_crit = "sic")
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
reliever_em <- function(data, cpn_max = 3, dm = 50, cov_rate = 0.8,
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
                        family, size = NULL) {
  if (missing(family)) {
    stop("family is required for reliever_em().", call. = FALSE)
  }
  .reliever_with_reg_fun(
    data = data, reg_fun = reg_fun_em,
    reg_args = list(family = family, size = size),
    cpn_max = cpn_max, dm = dm, cov_rate = cov_rate, method = method,
    cpn_crit = cpn_crit, pen_val = pen_val, prune_value = prune_value,
    M = M, wbs_seed = wbs_seed, wbs_stop_crit = wbs_stop_crit,
    detail = detail, cache_backend = cache_backend, owner_key = owner_key,
    echo = echo, dc_grid_size = dc_grid_size,
    cache_profile = cache_profile, run_cpd_ids = run_cpd_ids,
    dc_grid = dc_grid
  )
}
