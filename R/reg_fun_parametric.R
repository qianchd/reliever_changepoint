.reliever_validate_parametric_matrix <- function(x, name = "data") {
  if (!is.numeric(x) || nrow(x) == 0L || ncol(x) == 0L ||
      anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain finite numeric observations.", call. = FALSE)
  }
  invisible(x)
}

.reliever_parametric_matrix <- function(data, name = "data") {
  x <- if (is.vector(data)) matrix(data, ncol = 1L) else as.matrix(data)
  if (!is.numeric(x) || nrow(x) == 0L || ncol(x) == 0L) {
    stop(name, " must be a non-empty numeric vector or matrix.", call. = FALSE)
  }
  x
}

.reliever_prepare_parametric_data <- function(
    data, response_ncol = NULL, intercept = NULL, family = NULL) {
  x <- .reliever_parametric_matrix(data)
  .reliever_validate_parametric_matrix(x)

  if (!is.null(response_ncol)) {
    response_ncol <- .reliever_validate_positive_integer(
      response_ncol, "response_ncol"
    )
    if (response_ncol >= ncol(x)) {
      stop("Regression data must contain at least one predictor column.",
           call. = FALSE)
    }
    if (!is.null(family)) {
      resolved_family <- .reliever_resolve_glm_family(family)
      .reliever_glm_response(
        x[, seq_len(response_ncol), drop = FALSE], resolved_family
      )
    }
  }
  x
}

.reliever_parametric_virtual_info <- function(row_type, loss_kind) {
  list(
    n_loss_outputs = 1L,
    loss_output_meta = data.frame(
      loss_output_id = 1L,
      row_type = row_type,
      loss_kind = loss_kind,
      stringsAsFactors = FALSE
    )
  )
}

.reliever_gaussian_model <- function(x, center) {
  center <- as.numeric(center)
  if (length(center) == 1L && ncol(x) > 1L) {
    center <- rep(center, ncol(x))
  }
  if (length(center) != ncol(x) || anyNA(center) ||
      any(!is.finite(center))) {
    stop("The Gaussian center must be one finite value or one per column.",
         call. = FALSE)
  }
  centered <- sweep(x, 2L, center, "-")
  covariance <- crossprod(centered) / nrow(centered)
  covariance <- matrix(covariance, nrow = ncol(x), ncol = ncol(x))
  chol_covariance <- tryCatch(
    chol(covariance),
    error = function(e) NULL
  )
  list(
    center = center,
    covariance = covariance,
    chol = chol_covariance,
    log_det = if (is.null(chol_covariance)) {
      Inf
    } else {
      2 * sum(log(diag(chol_covariance)))
    }
  )
}

.reliever_gaussian_loss <- function(x, model) {
  if (is.null(model$chol) || !is.finite(model$log_det)) {
    return(rep(Inf, nrow(x)))
  }
  centered <- sweep(x, 2L, model$center, "-")
  standardized <- backsolve(
    model$chol, t(centered), transpose = TRUE
  )
  ncol(x) * log(2 * pi) + model$log_det +
    colSums(standardized^2)
}

.reliever_prepare_parametric_regression <- function(
    X, y, cpd_family, family_args) {
  if (!cpd_family %in% c("lm", "glm")) {
    stop("Parametric regression input requires cpd_family = \"lm\" or \"glm\".",
         call. = FALSE)
  }
  input_spec <- .reliever_response_input_spec(X, y)

  if (identical(cpd_family, "glm") && !is.null(y)) {
    response <- as.matrix(y)
    predictors <- as.matrix(X)
    if (!is.numeric(response) || anyNA(response) ||
        any(!is.finite(response)) || ncol(response) > 2L) {
      stop("GLM y must have one numeric column, or two binomial count columns.",
           call. = FALSE)
    }
    if (nrow(response) != nrow(predictors)) {
      stop("X and y must have the same number of observations.", call. = FALSE)
    }
    if (!is.null(family_args$response_ncol)) {
      supplied_ncol <- .reliever_validate_positive_integer(
        family_args$response_ncol, "response_ncol"
      )
      if (supplied_ncol != ncol(response)) {
        stop("response_ncol must match the number of columns in y.",
             call. = FALSE)
      }
    }
    family_args$response_ncol <- ncol(response)
    data <- cbind(response, predictors)
    input_spec <- list(
      type = "response_predictor",
      n_predictors = ncol(predictors),
      n_response_columns = ncol(response),
      original_form = "separate_xy"
    )
  } else {
    if (identical(cpd_family, "lm") && !is.null(y) &&
        ncol(as.matrix(y)) != 1L) {
      stop("LM y must contain one response column.", call. = FALSE)
    }
    data <- .reliever_prepare_input(
      X, y, cpd_family, response_families = c("lm", "glm")
    )
    response_ncol <- if (identical(cpd_family, "glm") &&
                         !is.null(family_args$response_ncol)) {
      .reliever_validate_positive_integer(
        family_args$response_ncol, "response_ncol"
      )
    } else {
      1L
    }
    if (response_ncol >= ncol(as.matrix(data))) {
      stop(
        toupper(cpd_family),
        " data must contain at least one predictor column.",
        call. = FALSE
      )
    }
    input_spec$n_response_columns <- response_ncol
    input_spec$n_predictors <- ncol(as.matrix(data)) - response_ncol
  }

  list(data = data, input_spec = input_spec, family_args = family_args)
}

#' Covariance-change interval loss
#'
#' Fit one covariance matrix on rows \code{l:r}, using a fixed common mean.
#' Return observation-level negative twice Gaussian log-likelihood on rows
#' \code{l_end:r_end}. If \code{mu = NULL}, use the mean of all supplied data.
#'
#' This is the low-level interval loss used by
#' \code{\link{reliever_var}()}; most users should call that wrapper.
#'
#' @param data Numeric vector or matrix with observations in rows.
#' @param l,r Inclusive fitting-interval endpoints.
#' @param l_end,r_end Inclusive evaluation-interval endpoints.
#' @param save_model Retain fitted parameters in the returned \code{model}.
#' @param is_virtual_run Return loss-output metadata without fitting.
#' @param mu Optional fixed common mean. A scalar is applied to every column.
#' @return A list containing the observation-level \code{loss} matrix and,
#'   when requested, the fitted \code{model}.
#' @seealso \code{\link{reliever_var}()}, \code{\link{reg_fun_meanvar}()},
#'   \code{\link{reg_fun_mean}()}
#' @examples
#' set.seed(2026)
#' x <- matrix(rnorm(60), ncol = 2)
#' out <- reg_fun_var(
#'   x, l = 1, r = 20, l_end = 21, r_end = 30,
#'   mu = c(0, 0), save_model = TRUE
#' )
#' head(out$loss)
#' out$model$covariance
#' @export
reg_fun_var <- function(data, l, r, l_end = l, r_end = r,
                        save_model = FALSE, is_virtual_run = FALSE,
                        mu = NULL) {
  if (is_virtual_run) {
    return(.reliever_parametric_virtual_info("var", "gaussian_nll"))
  }
  x <- .reliever_parametric_matrix(data)
  if (is.null(mu)) {
    mu <- .reliever_stable_col_means(x)
  }
  model <- .reliever_gaussian_model(
    x[l:r, , drop = FALSE], center = mu
  )
  loss <- .reliever_gaussian_loss(
    x[l_end:r_end, , drop = FALSE], model
  )
  list(
    loss = matrix(loss, ncol = 1L),
    model = if (isTRUE(save_model)) {
      model[c("center", "covariance")]
    } else {
      NULL
    }
  )
}

#' Simultaneous mean-and-covariance interval loss
#'
#' Estimate both the mean and maximum-likelihood covariance matrix on rows
#' \code{l:r}, then return observation-level negative twice multivariate
#' Gaussian log-likelihood on rows \code{l_end:r_end}.
#'
#' This is the low-level interval loss used by
#' \code{\link{reliever_meanvar}()}; most users should call that wrapper.
#'
#' @inheritParams reg_fun_var
#' @return A list containing the observation-level \code{loss} matrix and,
#'   when requested, the fitted \code{model}.
#' @seealso \code{\link{reliever_meanvar}()}, \code{\link{reg_fun_var}()}
#' @examples
#' set.seed(2026)
#' x <- matrix(rnorm(60), ncol = 2)
#' out <- reg_fun_meanvar(
#'   x, l = 1, r = 20, l_end = 21, r_end = 30, save_model = TRUE
#' )
#' head(out$loss)
#' out$model$center
#' @export
reg_fun_meanvar <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE, is_virtual_run = FALSE) {
  if (is_virtual_run) {
    return(.reliever_parametric_virtual_info(
      "meanvar", "gaussian_nll"
    ))
  }
  x <- .reliever_parametric_matrix(data)
  fit_x <- x[l:r, , drop = FALSE]
  model <- .reliever_gaussian_model(
    fit_x, center = .reliever_stable_col_means(fit_x)
  )
  loss <- .reliever_gaussian_loss(
    x[l_end:r_end, , drop = FALSE], model
  )
  list(
    loss = matrix(loss, ncol = 1L),
    model = if (isTRUE(save_model)) {
      model[c("center", "covariance")]
    } else {
      NULL
    }
  )
}

.reliever_regression_design <- function(data, response_ncol = 1L,
                                        intercept = TRUE) {
  response_ncol <- .reliever_validate_positive_integer(
    response_ncol, "response_ncol"
  )
  if (response_ncol >= ncol(data)) {
    stop("Regression data must contain at least one predictor column.",
         call. = FALSE)
  }
  x <- data[, -(seq_len(response_ncol)), drop = FALSE]
  if (isTRUE(intercept)) {
    x <- cbind(`(Intercept)` = 1, x)
  }
  list(
    response = data[, seq_len(response_ncol), drop = FALSE],
    design = x
  )
}

#' Linear-model interval loss
#'
#' Fit ordinary least squares on rows \code{l:r} and return squared residuals
#' for rows \code{l_end:r_end}. The first column of \code{data} is the response
#' and all remaining columns are predictors.
#'
#' This is the low-level interval loss used by
#' \code{\link{reliever_lm}()}; most users should call that wrapper.
#'
#' @inheritParams reg_fun_var
#' @param data Numeric matrix with the response first and predictors afterward.
#' @param intercept Add an intercept to the interval model.
#' @return A list containing the observation-level \code{loss} matrix and,
#'   when requested, the fitted \code{model}.
#' @seealso \code{\link{reliever_lm}()}, \code{\link{reg_fun_glm}()}
#' @examples
#' set.seed(2026)
#' x <- matrix(rnorm(60), ncol = 2)
#' y <- 1 + x[, 1] - x[, 2] + rnorm(30, sd = 0.3)
#' out <- reg_fun_lm(
#'   cbind(y, x), l = 1, r = 20, l_end = 21, r_end = 30,
#'   save_model = TRUE
#' )
#' head(out$loss)
#' out$model$coefficients
#' @export
reg_fun_lm <- function(data, l, r, l_end = l, r_end = r,
                       save_model = FALSE, is_virtual_run = FALSE,
                       intercept = TRUE) {
  if (is_virtual_run) {
    return(.reliever_parametric_virtual_info("lm", "rss"))
  }
  x <- .reliever_parametric_matrix(data)
  involved <- sort(unique(c(seq.int(l, r), seq.int(l_end, r_end))))
  parsed <- .reliever_regression_design(
    x[involved, , drop = FALSE],
    response_ncol = 1L,
    intercept = intercept
  )
  y <- as.numeric(parsed$response[, 1L])
  fit_id <- match(seq.int(l, r), involved)
  eval_id <- match(seq.int(l_end, r_end), involved)
  fit <- stats::lm.fit(
    x = parsed$design[fit_id, , drop = FALSE],
    y = y[fit_id]
  )
  coefficients <- as.numeric(fit$coefficients)
  coefficients[is.na(coefficients)] <- 0
  fitted <- drop(
    parsed$design[eval_id, , drop = FALSE] %*% coefficients
  )
  loss <- (y[eval_id] - fitted)^2
  loss[!is.finite(loss)] <- Inf
  list(
    loss = matrix(loss, ncol = 1L),
    model = if (isTRUE(save_model)) {
      list(coefficients = coefficients, intercept = isTRUE(intercept))
    } else {
      NULL
    }
  )
}

.reliever_resolve_glm_family <- function(family) {
  if (is.character(family) && length(family) == 1L && !is.na(family)) {
    family_fun <- tryCatch(
      getExportedValue("stats", family),
      error = function(e) NULL
    )
    if (!is.function(family_fun)) {
      stop("Unknown stats GLM family: ", family, ".", call. = FALSE)
    }
    family <- family_fun()
  } else if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "family") ||
      !is.function(family$linkinv) ||
      !is.function(family$dev.resids)) {
    stop("family must be a valid stats GLM family specification.",
         call. = FALSE)
  }
  family
}

.reliever_glm_response <- function(response, family) {
  if (ncol(response) == 2L) {
    if (!identical(family$family, "binomial")) {
      stop("A two-column response is supported only for binomial GLMs.",
           call. = FALSE)
    }
    if (any(response < 0) || any(response != floor(response))) {
      stop("Binomial successes and failures must be non-negative integers.",
           call. = FALSE)
    }
    weights <- rowSums(response)
    if (any(weights <= 0)) {
      stop("Every binomial response row must contain at least one trial.",
           call. = FALSE)
    }
    return(list(y = response[, 1L] / weights, weights = weights))
  }
  if (ncol(response) != 1L) {
    stop("GLM responses must have one column, or two for binomial counts.",
         call. = FALSE)
  }
  list(y = as.numeric(response[, 1L]), weights = rep(1, nrow(response)))
}

#' Generalized-linear-model interval loss
#'
#' Fit \code{stats::glm.fit()} on rows \code{l:r} and return the selected
#' family's observation-level deviance residuals for rows
#' \code{l_end:r_end}. The leading \code{response_ncol} columns are the
#' response and the remaining columns are predictors. For grouped binomial
#' data, provide successes and failures in the first two columns and use
#' \code{response_ncol = 2}.
#'
#' This is the low-level interval loss used by
#' \code{\link{reliever_glm}()}; most users should call that wrapper.
#'
#' @inheritParams reg_fun_var
#' @param data Numeric matrix with the response column or columns first and
#'   predictors afterward.
#' @param family A GLM family object, family function, or the name of a
#'   \code{stats} family function.
#' @param intercept Add an intercept to the interval model.
#' @param response_ncol Number of leading response columns; use 2 for binomial
#'   successes and failures.
#' @return A list containing the observation-level \code{loss} matrix and,
#'   when requested, the fitted \code{model}.
#' @seealso \code{\link{reliever_glm}()}, \code{\link{reg_fun_lm}()}
#' @examples
#' set.seed(2026)
#' x <- rnorm(30)
#' trials <- rep(5L, 30)
#' success <- rbinom(30, trials, plogis(x))
#' out <- reg_fun_glm(
#'   cbind(success, trials - success, x),
#'   l = 1, r = 20, l_end = 21, r_end = 30,
#'   family = binomial(), response_ncol = 2, save_model = TRUE
#' )
#' head(out$loss)
#' out$model$coefficients
#' @export
reg_fun_glm <- function(data, l, r, l_end = l, r_end = r,
                        save_model = FALSE, is_virtual_run = FALSE,
                        family = stats::gaussian(), intercept = TRUE,
                        response_ncol = 1L) {
  if (is_virtual_run) {
    return(.reliever_parametric_virtual_info("glm", "glm_deviance"))
  }
  family <- .reliever_resolve_glm_family(family)
  x <- .reliever_parametric_matrix(data)
  involved <- sort(unique(c(seq.int(l, r), seq.int(l_end, r_end))))
  parsed <- .reliever_regression_design(
    x[involved, , drop = FALSE],
    response_ncol = response_ncol,
    intercept = intercept
  )
  response <- .reliever_glm_response(parsed$response, family)
  fit_id <- match(seq.int(l, r), involved)
  fit <- suppressWarnings(stats::glm.fit(
    x = parsed$design[fit_id, , drop = FALSE],
    y = response$y[fit_id],
    weights = response$weights[fit_id],
    family = family
  ))
  coefficients <- as.numeric(fit$coefficients)
  eval_id <- match(seq.int(l_end, r_end), involved)
  if (any(!is.finite(coefficients))) {
    loss <- rep(Inf, length(eval_id))
  } else {
    eta <- drop(parsed$design[eval_id, , drop = FALSE] %*% coefficients)
    fitted <- family$linkinv(eta)
    loss <- family$dev.resids(
      response$y[eval_id], fitted, response$weights[eval_id]
    )
    loss[!is.finite(loss)] <- Inf
  }
  list(
    loss = matrix(loss, ncol = 1L),
    model = if (isTRUE(save_model)) {
      list(
        coefficients = coefficients,
        family = family,
        intercept = isTRUE(intercept),
        response_ncol = as.integer(response_ncol),
        converged = isTRUE(fit$converged)
      )
    } else {
      NULL
    }
  )
}

.reliever_match_em_family <- function(family) {
  aliases <- c(
    binomial = "binom",
    multinomial = "multinom",
    poisson = "pois",
    exponential = "exp",
    geometric = "geom",
    dirichlet = "diri",
    chi_square = "chisq",
    inverse_gaussian = "invgauss"
  )
  if (!is.character(family) || length(family) != 1L ||
      is.na(family) || !nzchar(family)) {
    stop("family must be one supported exponential-family name.",
         call. = FALSE)
  }
  if (family %in% names(aliases)) {
    family <- unname(aliases[[family]])
  }
  match.arg(
    family,
    c(
      "binom", "multinom", "pois", "exp", "geom", "diri",
      "gamma", "beta", "chisq", "invgauss"
    )
  )
}

.reliever_positive_optim <- function(start, objective) {
  fit <- suppressWarnings(stats::optim(
    par = log(start),
    fn = function(log_param) objective(exp(log_param)),
    method = "BFGS"
  ))
  if (!is.finite(fit$value)) {
    stop("Distribution parameter estimation did not produce a finite fit.",
         call. = FALSE)
  }
  exp(fit$par)
}

.reliever_fit_em <- function(x, family, size) {
  eps <- sqrt(.Machine$double.eps)
  scalar <- as.numeric(x[, 1L])
  if (family == "binom") {
    size <- .reliever_validate_positive_integer(size, "size")
    if (ncol(x) != 1L || any(scalar != floor(scalar)) ||
        any(scalar < 0 | scalar > size)) {
      stop("Binomial data must be integer counts between 0 and size.",
           call. = FALSE)
    }
    param <- min(1 - eps, max(eps, mean(scalar) / size))
  } else if (family == "multinom") {
    size <- .reliever_validate_positive_integer(size, "size")
    if (any(x != floor(x)) || any(x < 0) ||
        any(rowSums(x) != size)) {
      stop("Each multinomial row must contain non-negative integer counts summing to size.",
           call. = FALSE)
    }
    param <- colSums(x) / sum(x)
    param <- pmax(param, eps)
    param <- param / sum(param)
  } else if (family == "pois") {
    if (ncol(x) != 1L || any(scalar != floor(scalar)) || any(scalar < 0)) {
      stop("Poisson data must contain non-negative integers.",
           call. = FALSE)
    }
    param <- max(eps, mean(scalar))
  } else if (family == "exp") {
    if (ncol(x) != 1L || any(scalar < 0)) {
      stop("Exponential data must be non-negative.", call. = FALSE)
    }
    param <- 1 / max(eps, mean(scalar))
  } else if (family == "geom") {
    if (ncol(x) != 1L || any(scalar != floor(scalar)) || any(scalar < 0)) {
      stop("Geometric data must contain non-negative integers.",
           call. = FALSE)
    }
    param <- min(1 - eps, max(eps, 1 / (1 + mean(scalar))))
  } else if (family == "diri") {
    if (ncol(x) < 2L || any(x <= 0) ||
        any(abs(rowSums(x) - 1) > 1e-7)) {
      stop("Dirichlet rows must be positive vectors summing to one.",
           call. = FALSE)
    }
    mean_x <- colMeans(x)
    component_var <- apply(x, 2L, stats::var)
    concentration <- stats::median(
      mean_x * (1 - mean_x) / pmax(component_var, eps) - 1
    )
    start <- pmax(eps, mean_x * max(1, concentration))
    param <- .reliever_positive_optim(start, function(alpha) {
      -sum(
        lgamma(sum(alpha)) - sum(lgamma(alpha)) +
          rowSums(sweep(log(x), 2L, alpha - 1, "*"))
      )
    })
  } else if (family == "gamma") {
    if (ncol(x) != 1L || any(scalar <= 0)) {
      stop("Gamma data must be positive.", call. = FALSE)
    }
    variance <- max(eps, stats::var(scalar))
    start <- c(mean(scalar)^2 / variance, mean(scalar) / variance)
    param <- .reliever_positive_optim(start, function(theta) {
      -sum(stats::dgamma(
        scalar, shape = theta[1L], rate = theta[2L], log = TRUE
      ))
    })
  } else if (family == "beta") {
    if (ncol(x) != 1L || any(scalar <= 0 | scalar >= 1)) {
      stop("Beta data must lie strictly between zero and one.",
           call. = FALSE)
    }
    mean_x <- mean(scalar)
    concentration <- max(
      2, mean_x * (1 - mean_x) / max(eps, stats::var(scalar)) - 1
    )
    start <- c(mean_x, 1 - mean_x) * concentration
    param <- .reliever_positive_optim(start, function(theta) {
      -sum(stats::dbeta(
        scalar, shape1 = theta[1L], shape2 = theta[2L], log = TRUE
      ))
    })
  } else if (family == "chisq") {
    if (ncol(x) != 1L || any(scalar <= 0)) {
      stop("Chi-square data must be positive.", call. = FALSE)
    }
    param <- .reliever_positive_optim(max(eps, mean(scalar)), function(df) {
      -sum(stats::dchisq(scalar, df = df, log = TRUE))
    })
  } else {
    if (ncol(x) != 1L || any(scalar <= 0)) {
      stop("Inverse-Gaussian data must be positive.", call. = FALSE)
    }
    mu <- mean(scalar)
    denominator <- sum((scalar - mu)^2 / (mu^2 * scalar))
    param <- c(mu = mu, lambda = nrow(x) / max(eps, denominator))
  }
  list(family = family, param = as.numeric(param), size = size)
}

.reliever_em_loss <- function(x, model) {
  family <- model$family
  param <- model$param
  scalar <- as.numeric(x[, 1L])
  log_likelihood <- if (family == "binom") {
    stats::dbinom(scalar, size = model$size, prob = param, log = TRUE)
  } else if (family == "multinom") {
    vapply(seq_len(nrow(x)), function(i) {
      stats::dmultinom(
        x[i, ], size = model$size, prob = param, log = TRUE
      )
    }, numeric(1L))
  } else if (family == "pois") {
    stats::dpois(scalar, lambda = param, log = TRUE)
  } else if (family == "exp") {
    stats::dexp(scalar, rate = param, log = TRUE)
  } else if (family == "geom") {
    stats::dgeom(scalar, prob = param, log = TRUE)
  } else if (family == "diri") {
    lgamma(sum(param)) - sum(lgamma(param)) +
      rowSums(sweep(log(x), 2L, param - 1, "*"))
  } else if (family == "gamma") {
    stats::dgamma(
      scalar, shape = param[1L], rate = param[2L], log = TRUE
    )
  } else if (family == "beta") {
    stats::dbeta(
      scalar, shape1 = param[1L], shape2 = param[2L], log = TRUE
    )
  } else if (family == "chisq") {
    stats::dchisq(scalar, df = param, log = TRUE)
  } else {
    mu <- param[1L]
    lambda <- param[2L]
    0.5 * (
      log(lambda / (2 * pi)) - 3 * log(scalar) -
        lambda * (scalar - mu)^2 / (mu^2 * scalar)
    )
  }
  loss <- -2 * log_likelihood
  loss[!is.finite(loss)] <- Inf
  loss
}

#' Exponential-family interval loss
#'
#' Estimate one distribution on rows \code{l:r} and return observation-level
#' negative twice log-likelihood on rows \code{l_end:r_end}. Supported family
#' names are grouped as follows:
#' \itemize{
#'   \item \code{"binom"}, \code{"multinom"}, and \code{"pois"};
#'   \item \code{"exp"}, \code{"geom"}, and \code{"diri"};
#'   \item \code{"gamma"}, \code{"beta"}, \code{"chisq"}, and
#'   \code{"invgauss"}.
#' }
#' Full distribution names are also accepted; see \code{family} below.
#'
#' This is the low-level interval loss used by
#' \code{\link{reliever_em}()}; most users should call that wrapper.
#'
#' @inheritParams reg_fun_var
#' @param family Exponential-family distribution name.
#' @param size Number of trials, required for binomial and multinomial data.
#' @return A list containing the observation-level \code{loss} matrix and,
#'   when requested, the fitted \code{model}.
#' @seealso \code{\link{reliever_em}()}, \code{\link{reg_fun_glm}()}
#' @examples
#' set.seed(2026)
#' x <- rexp(30, rate = 2)
#' out <- reg_fun_em(
#'   x, l = 1, r = 20, l_end = 21, r_end = 30,
#'   family = "exp", save_model = TRUE
#' )
#' head(out$loss)
#' out$model$param
#' @export
reg_fun_em <- function(data, l, r, l_end = l, r_end = r,
                       save_model = FALSE, is_virtual_run = FALSE,
                       family, size = NULL) {
  if (missing(family)) {
    stop("family is required for reg_fun_em().", call. = FALSE)
  }
  family <- .reliever_match_em_family(family)
  if (is_virtual_run) {
    return(.reliever_parametric_virtual_info(
      paste0("em_", family), "negative_log_likelihood"
    ))
  }
  x <- .reliever_parametric_matrix(data)
  model <- .reliever_fit_em(x[l:r, , drop = FALSE], family, size)
  loss <- .reliever_em_loss(x[l_end:r_end, , drop = FALSE], model)
  list(
    loss = matrix(loss, ncol = 1L),
    model = if (isTRUE(save_model)) model else NULL
  )
}
