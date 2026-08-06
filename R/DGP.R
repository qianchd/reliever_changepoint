.dgp_validate_ar1 <- function(rho_ar1) {
  if (!is.numeric(rho_ar1) || length(rho_ar1) != 1L ||
      is.na(rho_ar1) || abs(rho_ar1) >= 1) {
    stop("rho_ar1 must be a single number in (-1, 1).", call. = FALSE)
  }
  rho_ar1
}

.dgp_validate_linear_inputs <- function(n, p, tau, b0, delta, rho) {
  n <- .reliever_validate_positive_integer(n, "n")
  p <- .reliever_validate_positive_integer(p, "p")
  if (!is.numeric(tau) || anyNA(tau) || any(!is.finite(tau)) ||
      any(tau != floor(tau))) {
    stop("tau must contain integer changepoint locations.", call. = FALSE)
  }
  tau <- as.integer(tau)
  if (length(tau) > 0L &&
      (any(tau <= 0L | tau >= n) || is.unsorted(tau, strictly = TRUE))) {
    stop("tau must be strictly increasing with every value in (0, n).",
         call. = FALSE)
  }
  if (!is.numeric(b0) || length(b0) != p || anyNA(b0) ||
      any(!is.finite(b0))) {
    stop("b0 must contain p finite coefficients.", call. = FALSE)
  }
  if (!is.matrix(delta) || !is.numeric(delta) ||
      !identical(dim(delta), c(p, length(tau))) || anyNA(delta) ||
      any(!is.finite(delta))) {
    stop("delta must be a finite p by length(tau) numeric matrix.",
         call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) ||
      !is.finite(rho) || abs(rho) >= 1) {
    stop("rho must be a single number in (-1, 1).", call. = FALSE)
  }
  list(
    n = n,
    p = p,
    tau = tau,
    b0 = as.numeric(b0),
    delta = delta,
    rho = as.numeric(rho)
  )
}

.dgp_ar1_sequence <- function(innovations, rho_ar1) {
  if (rho_ar1 == 0) {
    return(innovations)
  }
  out <- innovations
  scale <- sqrt(1 - rho_ar1^2)
  if (is.null(dim(out))) {
    if (length(out) > 1L) {
      for (i in 2:length(out)) {
        out[i] <- rho_ar1 * out[i - 1] + scale * innovations[i]
      }
    }
    return(out)
  }
  if (NROW(out) > 1L) {
    for (i in 2:NROW(out)) {
      out[i, ] <- rho_ar1 * out[i - 1, ] + scale * innovations[i, ]
    }
  }
  out
}

#' Generate a linear regression sequence with coefficient changepoints
#'
#' Generate observations from a piecewise linear regression model. The response
#' is stored in the first column and predictors are stored in the remaining
#' columns, matching the response-first input convention used by
#' `reg_fun_lasso_solpath()`. For the primary built-in interface, pass the
#' predictor columns and response separately to
#' `reliever(X, y, cpd_family = "lasso")`.
#'
#' @param n Number of observations.
#' @param p Number of predictors.
#' @param tau Increasing integer changepoints. At \code{t}, the coefficient
#'   changes after observation \code{t}.
#' @param b0 Coefficient vector in the first segment, of length \code{p}.
#' @param delta A \code{p} by \code{length(tau)} matrix. Column \code{j} is
#'   added to the current coefficient vector after changepoint \code{tau[j]},
#'   so coefficient changes accumulate across segments.
#' @param sig Noise standard deviation.
#' @param rho Correlation between adjacent predictor coordinates. The
#'   within-observation predictor covariance has entries
#'   \eqn{rho^{|j-k|}}. Use 0 for independent predictor coordinates.
#' @param rho_ar1 AR(1) correlation over time, applied separately to each
#'   predictor coordinate and to the response noise. Predictors retain unit
#'   marginal variance before spatial correlation, and response innovations
#'   have unit variance before either \code{sig} scaling or segment-specific
#'   scaling. \code{rho_ar1 = 0} gives temporally independent observations.
#'
#' @return A list containing \code{data}, the response followed by predictors;
#'   \code{mu}, the noise-free conditional response mean; and \code{ep}, the
#'   response noise added to \code{mu}.
#' @seealso \code{\link{dgp_linear_regression_equal_response}()},
#'   \code{\link{reliever}()}, \code{\link{cv.reliever}()}
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' set.seed(2026)
#' n <- 900
#' p <- 100
#' tau <- c(300, 600)
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' generated <- dgp_linear_regression(n, p, tau, b0, delta, sig = 0.5)
#' dim(generated$data)
#'
#' set.seed(2026)
#' generated_ar <- dgp_linear_regression(
#'   n, p, tau, b0, delta, sig = 0.5, rho_ar1 = 0.3
#' )
#' dim(generated_ar$data)
dgp_linear_regression <- function(n, p, tau, b0, delta, sig = 1, rho = 0,
                                  rho_ar1 = 0) {
  setup <- .dgp_validate_linear_inputs(n, p, tau, b0, delta, rho)
  n <- setup$n
  p <- setup$p
  tau <- setup$tau
  beta <- setup$b0
  delta <- setup$delta
  rho <- setup$rho
  rho_ar1 <- .dgp_validate_ar1(rho_ar1)
  if (!is.numeric(sig) || length(sig) != 1L || is.na(sig) ||
      !is.finite(sig) || sig < 0) {
    stop("sig must be a non-negative finite number.", call. = FALSE)
  }
  x <- .dgp_ar1_sequence(matrix(rnorm(n * p), n, p), rho_ar1)
  if (rho != 0) {
    x <- x %*% chol(toeplitz(rho^c(0:(p - 1))))
  }

  mu <- rep(NA_real_, n)
  tau_ext <- c(0, tau, n)
  for (j in seq_len(length(tau) + 1L)) {
    idx_j <- seq.int(tau_ext[j] + 1L, tau_ext[j + 1L])
    if (j > 1) {
      beta <- beta + delta[, j - 1L]
    }
    mu[idx_j] <- drop(x[idx_j, , drop = FALSE] %*% beta)
  }
  ep <- as.vector(.dgp_ar1_sequence(rnorm(n), rho_ar1))
  ep <- sig * ep
  y <- mu + ep
  data <- cbind(y, x)
  colnames(data) <- NULL
  list(data = data, mu = mu, ep = ep)
}

#' Generate a linear regression sequence with equalized response variance
#'
#' Generate a piecewise linear regression sequence while choosing a separate
#' noise variance in each segment so that the marginal response variance is
#' constant across segments. This avoids making coefficient changes detectable
#' merely through a change in response variance.
#'
#' @inheritParams dgp_linear_regression
#' @param snr_y Positive response-level signal-to-noise control. Larger values
#'   reduce the common target response variance toward the largest segment-wise
#'   signal variance and therefore reduce added noise.
#'
#' @return A list containing \code{data}, the response followed by predictors;
#'   \code{mu}, the noise-free conditional response mean; and \code{ep}, the
#'   segment-scaled response noise.
#' @seealso \code{\link{dgp_linear_regression}()},
#'   \code{\link{reliever}()}, \code{\link{cv.reliever}()},
#'   \code{\link{reliever_lasso}()}
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' set.seed(2026)
#' n <- 900
#' p <- 100
#' tau <- c(300, 600)
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' generated <- dgp_linear_regression_equal_response(
#'   n, p, tau, b0, delta, snr_y = 2
#' )
#' dim(generated$data)
dgp_linear_regression_equal_response <- function(n, p, tau, b0, delta, rho = 0,
                                                 snr_y = 1, rho_ar1 = 0) {
  setup <- .dgp_validate_linear_inputs(n, p, tau, b0, delta, rho)
  n <- setup$n
  p <- setup$p
  tau <- setup$tau
  beta <- setup$b0
  delta <- setup$delta
  rho <- setup$rho
  rho_ar1 <- .dgp_validate_ar1(rho_ar1)
  if (!is.numeric(snr_y) || length(snr_y) != 1L || is.na(snr_y) ||
      !is.finite(snr_y) || snr_y <= 0) {
    stop("snr_y must be a positive finite number.", call. = FALSE)
  }
  x <- .dgp_ar1_sequence(matrix(rnorm(n * p), n, p), rho_ar1)

  if (rho != 0) {
    Sig_x <- toeplitz(rho^c(0:(p - 1)))
    x <- x %*% chol(Sig_x)
  } else {
    Sig_x <- diag(p)
  }

  mu <- rep(NA_real_, n)
  tau_ext <- c(0, tau, n)
  var_xb <- numeric(length(tau) + 1L)
  for (j in seq_len(length(tau) + 1L)) {
    idx_j <- seq.int(tau_ext[j] + 1L, tau_ext[j + 1L])
    if (j > 1) {
      beta <- beta + delta[, j - 1L]
    }
    var_xb[j] <- drop(crossprod(beta, Sig_x %*% beta))
    mu[idx_j] <- drop(x[idx_j, , drop = FALSE] %*% beta)
  }
  ep <- as.vector(.dgp_ar1_sequence(rnorm(n), rho_ar1))
  y <- rep(NA_real_, n)
  var_y <- max(var_xb) * (1 + 1 / snr_y^2)
  for (j in seq_len(length(tau) + 1L)) {
    idx_j <- seq.int(tau_ext[j] + 1L, tau_ext[j + 1L])
    ep[idx_j] <- ep[idx_j] * sqrt(max(0, var_y - var_xb[j]))
    y[idx_j] <- mu[idx_j] + ep[idx_j]
  }
  data <- cbind(y, x)
  colnames(data) <- NULL
  list(data = data, mu = mu, ep = ep)
}

#' Symmetric changepoint-location error
#'
#' Compute the symmetric Hausdorff distance between estimated and true
#' changepoint sets: every estimated point is matched to its nearest true point,
#' every true point is matched to its nearest estimate, and the largest of
#' these distances is returned.
#'
#' @param estimate Estimated changepoint locations.
#' @param truth True changepoint locations.
#'
#' @return A non-negative numeric error. The error is zero when both sets are
#'   empty and infinite when exactly one set is empty.
#' @seealso \code{\link{reliever}()}, \code{\link{reliever_generic}()},
#'   \code{\link{select_by_run}()}
#' @export
#'
#' @examples
#' cp_error(c(298, 603), c(300, 600))
cp_error <- function(estimate, truth) {
  if (length(estimate) == 0L && length(truth) == 0L) {
    return(0)
  }
  if (length(estimate) == 0L || length(truth) == 0L) {
    return(Inf)
  }
  diff_abs <- abs(outer(estimate, truth, "-"))
  ue <- max(apply(diff_abs, 1, min))
  oe <- max(apply(diff_abs, 2, min))
  max(ue, oe)
}
