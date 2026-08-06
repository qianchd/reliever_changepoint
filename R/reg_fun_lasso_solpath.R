.reliever_validate_lam_set <- function(lam_set) {
  if (!is.numeric(lam_set) || length(lam_set) == 0L ||
      anyNA(lam_set) || any(!is.finite(lam_set)) || any(lam_set <= 0)) {
    stop("lam_set must contain finite positive numeric values.", call. = FALSE)
  }
  sort(unique(as.numeric(lam_set)), decreasing = TRUE)
}

.reliever_validate_glmnet_family <- function(family) {
  supported <- c("gaussian", "binomial", "poisson")
  if (!is.character(family) || length(family) != 1L ||
      is.na(family) || !nzchar(family) || !family %in% supported) {
    stop(
      "family must be one of \"gaussian\", \"binomial\", or \"poisson\".",
      call. = FALSE
    )
  }
  family
}

.reliever_validate_lasso_response <- function(response, family,
                                              name = "response") {
  family <- .reliever_validate_glmnet_family(family)
  if (!is.numeric(response) || anyNA(response) ||
      any(!is.finite(response))) {
    stop(name, " must contain finite numeric values.", call. = FALSE)
  }
  if (identical(family, "binomial") && any(!response %in% c(0, 1))) {
    stop(name, " must contain only 0 and 1 for family = \"binomial\".",
         call. = FALSE)
  }
  if (identical(family, "poisson") && any(response < 0)) {
    stop(name, " must contain non-negative values for family = \"poisson\".",
         call. = FALSE)
  }
  invisible(response)
}

.reliever_lasso_data_matrix <- function(data) {
  data <- as.matrix(data)
  if (!is.numeric(data) || nrow(data) < 2L || ncol(data) < 3L) {
    stop(
      "Lasso losses require a numeric response-plus-predictor matrix with at least two observations and two predictors.",
      call. = FALSE
    )
  }
  data
}

.reliever_auto_lam_set <- function(data, family = "gaussian", thresh = 1e-7,
                                   nlambda = 30L) {
  data <- .reliever_lasso_data_matrix(data)
  family <- .reliever_validate_glmnet_family(family)
  .reliever_validate_lasso_response(data[, 1L], family)
  nlambda <- .reliever_validate_positive_integer(nlambda, "nlambda")
  if (nlambda < 2L) {
    stop("nlambda must be a single integer >= 2.", call. = FALSE)
  }

  model <- glmnet::glmnet(
    x = data[, -1, drop = FALSE],
    y = data[, 1],
    family = family,
    nlambda = nlambda,
    intercept = FALSE,
    standardize = FALSE,
    control = list(thresh = thresh)
  )
  base_path <- .reliever_validate_lam_set(model$lambda * sqrt(nrow(data)))
  if (length(base_path) < 2L) {
    stop(
      "glmnet returned fewer than two lambda values; supply lam_set explicitly.",
      call. = FALSE
    )
  }

  ratios <- base_path[-length(base_path)] / base_path[-1L]
  ratio <- stats::median(ratios[is.finite(ratios) & ratios > 1])
  if (!is.finite(ratio) || ratio <= 1) {
    stop(
      "The glmnet lambda path could not be extended; supply lam_set explicitly.",
      call. = FALSE
    )
  }
  extension <- seq_len(2L)
  .reliever_validate_lam_set(c(
    base_path[1L] * ratio^extension,
    base_path,
    base_path[length(base_path)] / ratio^extension
  ))
}

.reliever_resolve_lam_set <- function(data, lam_set, family, thresh,
                                      nlambda = 30L) {
  .reliever_validate_glmnet_family(family)
  if (is.null(lam_set)) {
    return(.reliever_auto_lam_set(
      data, family = family, thresh = thresh, nlambda = nlambda
    ))
  }
  .reliever_validate_lam_set(lam_set)
}

#' Lasso solution-path interval loss
#'
#' Fit one glmnet lasso solution path on rows \code{l:r}, then return the
#' observation-level prediction loss on \code{l_end:r_end} for every lambda.
#' This is the low-level loss function used by \code{reliever()} when
#' \code{cpd_family = "lasso"}, and by outer-CV lasso workflows.
#' The returned columns define candidates; an
#' in-sample minimum must not be used to select lambda. Use cross-fitted loss,
#' outer CV, or independent hold-out loss for that choice. For an exact fitted interval,
#' training loss ordinarily decreases as lambda decreases, so an interior
#' minimum is neither expected nor a valid grid check. A fixed-K held-out or
#' cross-fitted loss profile should instead attain its minimum away from both
#' ends of a sufficiently wide lambda grid.
#' The built-in model uses \code{intercept = FALSE} and
#' \code{standardize = FALSE}; predictors should already be on their intended
#' scale.
#'
#' @param data Numeric matrix with response in the first column and predictors
#'   in the remaining columns. At least two predictor columns are required by
#'   \code{glmnet}.
#' @param l,r Rows used to fit the lasso path.
#' @param l_end,r_end Rows for which prediction losses are returned.
#' @param save_model Return the fitted glmnet model.
#' @param is_virtual_run Query flag used by \code{reliever_generic()}. When
#'   \code{TRUE}, return lambda metadata without fitting glmnet.
#' @param lam_set Finite positive values on a sample-size-independent scale. On
#'   an interval containing \eqn{m} observations, the lambda passed to glmnet is
#'   \code{lam_set / sqrt(m)}. This low-level function requires an explicit
#'   path; \code{reliever(X, y, cpd_family = "lasso")} constructs one
#'   automatically when omitted.
#' @param family Response family: \code{"gaussian"}, \code{"binomial"}, or
#'   \code{"poisson"}. Gaussian returns squared prediction loss; binomial and
#'   Poisson return negative log-likelihood loss. Binomial responses must be
#'   0 or 1; Poisson responses may contain any finite non-negative values.
#' @param thresh Convergence threshold passed to \code{glmnet::glmnet()}.
#'
#' @return With \code{is_virtual_run = TRUE}, metadata linking each loss output
#'   to its lambda. Otherwise, a list with:
#'   \describe{
#'     \item{\code{loss}}{An observation-by-lambda matrix.}
#'     \item{\code{model}}{The fitted glmnet path when
#'     \code{save_model = TRUE}; otherwise \code{NULL}.}
#'   }
#' @seealso \code{\link{reliever}()} and
#'   \code{\link{reg_fun_lasso_crossfit}()}.
#' @export
#'
#' @examples
#' set.seed(2026)
#' n <- 900
#' p <- 100
#' tau <- c(300, 600)
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' data <- dgp_linear_regression(n, p, tau, b0, delta, sig = 1)$data
#' lam_set <- 0.1 * 1.23^(30:0)
#' out <- reg_fun_lasso_solpath(data, 1, 300, 1, 360, lam_set = lam_set)
#' dim(out$loss)
reg_fun_lasso_solpath <- function(data, l, r, l_end = l, r_end = r,
                          save_model = FALSE, is_virtual_run = FALSE,
                          lam_set = NULL,
                          family = "gaussian", thresh = 1e-7) {
  data <- .reliever_lasso_data_matrix(data)
  family <- .reliever_validate_glmnet_family(family)
  if (is.null(lam_set)) {
    stop(
      paste0(
        "lam_set is required by reg_fun_lasso_solpath(); use ",
        "reliever(..., cpd_family = \"lasso\") for an automatic path."
      ),
      call. = FALSE
    )
  }
  lam_set <- .reliever_validate_lam_set(lam_set)
  if (is_virtual_run) {
    return(list(
      n_loss_outputs = length(lam_set),
      loss_output_meta = data.frame(
        loss_output_id = seq_along(lam_set),
        row_type = "lasso",
        hyper_id = seq_along(lam_set),
        hyper_value = I(as.list(lam_set)),
        hyper_name = "lambda",
        loss_kind = if (identical(family, "gaussian")) {
          "rss"
        } else {
          "negative_log_likelihood"
        }
      )
    ))
  }
  involved <- sort(unique(c(seq.int(l, r), seq.int(l_end, r_end))))
  .reliever_validate_lasso_response(data[involved, 1L], family)
  y <- data[l:r, 1]
  x <- as.matrix(data[l:r, -1, drop = FALSE])
  model <- glmnet::glmnet(x, y, family = family,
                          lambda = lam_set / sqrt(r - l + 1L),
                          intercept = FALSE, standardize = FALSE,
                          control = list(thresh = thresh))

  x <- as.matrix(data[l_end:r_end, -1, drop = FALSE])
  y <- data[l_end:r_end, 1]
  loss <- .reliever_glmnet_loss(model, x, y, family)
  list(
    loss = loss,
    model = if (isTRUE(save_model)) model else NULL
  )
}
