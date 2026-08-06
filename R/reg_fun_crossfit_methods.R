# Mean-square crossfit reg_fun ------------------------------------------------

.loss_mean_crossfit <- function(data, train_id, eval_id, hyper_set, ...) {
  x_train <- data[train_id, , drop = FALSE]
  center <- .reliever_stable_col_means(x_train)
  x_eval <- data[eval_id, , drop = FALSE]
  matrix(rowMeans(sweep(x_eval, 2L, center, "-")^2), ncol = 1L)
}

#' Cross-fitted mean-square interval loss
#'
#' Construct interval-level cross-fitted and in-sample mean-square losses. In
#' each fold, the segment mean is estimated from the training observations and
#' squared error is returned for held-out observations. Most users should call
#' \code{reliever(X, cpd_family = "mean_crossfit")}; this function is the
#' lower-level loss function for \code{reliever_generic()} and custom
#' extensions.
#'
#' @param data Numeric vector or matrix with observations in rows.
#' @inheritParams reg_fun_crossfit_template
#'
#' @return A named list. With \code{is_virtual_run = TRUE}, it contains the
#'   integer scalar \code{n_loss_outputs} and data frame
#'   \code{loss_output_meta}. Otherwise, it contains a numeric
#'   observation-by-output \code{loss} matrix and \code{model = NULL}. The
#'   default output is \code{recv}, the
#'   cross-fitted mean-square loss; \code{incv} can additionally return its
#'   in-sample counterpart. See
#'   \code{\link{reg_fun_crossfit_template}()} for the full return convention.
#'
#' @seealso \code{\link{reliever_mean_crossfit}()},
#'   \code{\link{reliever_mean}()}, \code{\link{reg_fun_crossfit_template}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- c(
#'   rnorm(n_seg, mean = 0, sd = 0.5),
#'   rnorm(n_seg, mean = 4, sd = 0.5),
#'   rnorm(n_seg, mean = -4, sd = 0.5)
#' )
#' out <- reg_fun_mean_crossfit(
#'   data = x, l = 241, r = 660, l_end = 211, r_end = 690, nfolds = 2
#' )
#' colSums(out$loss)
reg_fun_mean_crossfit <- function(data, l, r, l_end = l, r_end = r,
                            nfolds = 5,
                            fold_type = c("op", "blk", "stable_blk"),
                            op_size = 1,
                            buffer_lag = 0,
                            fold_stable_const = 1,
                            loss_output_types = "recv",
                            save_model = FALSE, is_virtual_run = FALSE) {
  external_train_n <- attr(
    data, ".reliever_external_train_n", exact = TRUE
  )
  if (is_virtual_run) {
    meta <- .crossfit_loss_output_meta(
      list(list(method = "mean")),
      loss_output_types = loss_output_types
    )
    meta$loss_kind <- ifelse(
      meta$row_type == "incv", "rss", "crossfit_loss"
    )
    return(list(
      n_loss_outputs = nrow(meta),
      loss_output_meta = meta
    ))
  }
  data <- if (is.vector(data)) matrix(data, ncol = 1L) else as.matrix(data)
  if (!is.null(external_train_n)) {
    attr(data, ".reliever_external_train_n") <- external_train_n
  }
  out <- reg_fun_crossfit_template(
    data = data,
    l = l,
    r = r,
    l_end = l_end,
    r_end = r_end,
    loss_fun = .loss_mean_crossfit,
    hyper_set = list(list(method = "mean")),
    loss_output_types = loss_output_types,
    nfolds = nfolds,
    fold_type = match.arg(fold_type),
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = save_model,
    is_virtual_run = FALSE
  )
  out
}

# KDE negative-log-likelihood crossfit reg_fun --------------------------------

.kde_bandwidth_vec <- function(hyper_set) {
  if (is.data.frame(hyper_set) && "bandwidth" %in% names(hyper_set)) {
    return(hyper_set$bandwidth)
  }
  if (is.matrix(hyper_set) && "bandwidth" %in% colnames(hyper_set)) {
    return(hyper_set[, "bandwidth"])
  }
  as.numeric(unlist(hyper_set, use.names = FALSE))
}

.loss_kde_nll_crossfit <- function(data, train_id, eval_id, hyper_set,
                             var_dim, kernel, kernel_args,
                             distance_power, ...) {
  bandwidth_vec <- .kde_bandwidth_vec(hyper_set)
  dist_sub <- data[eval_id, train_id, drop = FALSE]
  .kde_nll_loss(
    dist_sub, bandwidth_vec, var_dim,
    kernel = kernel, kernel_args = kernel_args,
    distance_power = distance_power
  )
}

.kde_nll_crossfit_external <- function(data, l, r, l_end, r_end,
                                        bandwidth_vec, var_dim,
                                        kernel, kernel_args,
                                        distance_power,
                                        nfolds, fold_type, op_size,
                                        buffer_lag, fold_stable_const,
                                        loss_output_types) {
  n_train <- attr(data, ".reliever_external_train_n", exact = TRUE)
  if (!is.numeric(n_train) || length(n_train) != 1L || is.na(n_train) ||
      n_train < 1L || n_train != floor(n_train)) {
    stop("External KDE evaluation has an invalid training-row marker.",
         call. = FALSE)
  }
  n_train <- as.integer(n_train)
  if (!is.matrix(data) || ncol(data) != n_train ||
      nrow(data) <= n_train) {
    stop(
      "External KDE-NLL crossfit evaluation requires a training-distance ",
      "matrix followed by evaluation-by-training distance rows.",
      call. = FALSE
    )
  }
  if (l < 1L || r > n_train || l > r ||
      l_end <= n_train || r_end > nrow(data) || l_end > r_end) {
    stop("External KDE-NLL crossfit evaluation has invalid interval endpoints.",
         call. = FALSE)
  }
  train_data <- data[seq_len(n_train), , drop = FALSE]

  training_fit <- reg_fun_crossfit_template(
    data = train_data,
    l = l,
    r = r,
    l_end = l,
    r_end = r,
    loss_fun = .loss_kde_nll_crossfit,
    hyper_set = bandwidth_vec,
    hyper_name = "bandwidth",
    loss_output_types = c("recv", "crossfit_homo_hyper"),
    nfolds = nfolds,
    fold_type = fold_type,
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = FALSE,
    is_virtual_run = FALSE,
    var_dim = var_dim,
    kernel = kernel,
    kernel_args = kernel_args,
    distance_power = distance_power
  )

  external_fixed <- .loss_kde_nll_crossfit(
    data = data,
    train_id = l:r,
    eval_id = l_end:r_end,
    hyper_set = bandwidth_vec,
    var_dim = var_dim,
    kernel = kernel,
    kernel_args = kernel_args,
    distance_power = distance_power
  )
  fixed_cf_columns <- which(
    colnames(training_fit$loss) == "crossfit_homo_hyper"
  )
  .crossfit_loss_outputs(
    external_fixed,
    external_fixed,
    selection_loss = training_fit$loss[, fixed_cf_columns, drop = FALSE],
    loss_output_types = loss_output_types
  )
}

#' Cross-fitted KDE negative-log-likelihood interval loss
#'
#' Construct interval-level cross-fitted KDE negative-log-likelihood losses for
#' a path of bandwidths and one fixed isotropic density kernel. The low-level
#' input is a squared-distance matrix; each fold estimates density from its
#' training rows and evaluates held-out rows. Most users should pass the
#' original observations to
#' \code{reliever(X, cpd_family = "kde_nll_crossfit")}. Supplying a
#' precomputed squared-distance matrix together with \code{var_dim} remains
#' supported.
#'
#' @inheritParams reg_fun_kde_nll_solpath
#' @inheritParams reg_fun_crossfit_template
#' @param data Numeric \eqn{n} by \eqn{n} matrix of pairwise squared distances.
#' @param l,r Integer scalars giving the first and last rows used to fit the
#'   interval KDE.
#' @param l_end,r_end Integer scalars giving the first and last rows whose
#'   negative log densities are returned.
#' @param save_model A logical scalar accepted for compatibility with the
#'   interval-loss protocol. KDE models are not retained.
#' @param is_virtual_run A logical scalar used as a metadata-query flag. When
#'   \code{TRUE}, return only the number and metadata of loss outputs.
#' @param bandwidth_vec \code{NULL} or a numeric vector of positive KDE
#'   bandwidths. \code{NULL} uses the scale-adaptive default grid.
#' @param var_dim \code{NULL} or a positive integer scalar giving the dimension
#'   of the original observations. It is required when it cannot be inferred.
#' @param kernel A character scalar selecting \code{"gaussian"},
#'   \code{"laplace"}, or \code{"student"}.
#' @param kernel_args A named list of kernel-specific settings.
#' @param distance_power An integer scalar equal to 1 or 2 identifying whether
#'   \code{data} stores distances or squared distances.
#'
#' @return A named list. With \code{is_virtual_run = TRUE}, it contains the
#'   integer scalar \code{n_loss_outputs} and data frame
#'   \code{loss_output_meta}. Otherwise, it contains the requested numeric
#'   observation-by-output \code{loss} matrix and \code{model = NULL}. See
#'   \code{\link{reg_fun_crossfit_template}()} for the full return convention
#'   and column definitions.
#'
#' @seealso \code{\link{reliever_kde_nll_crossfit}()} and
#'   \code{\link{reg_fun_crossfit_template}()}.
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
#' dist_sq <- as.matrix(dist(x))^2
#' out <- reg_fun_kde_nll_crossfit(
#'   data = dist_sq, l = 241, r = 660, l_end = 211, r_end = 690,
#'   var_dim = 5, nfolds = 2
#' )
#' colSums(out$loss)
reg_fun_kde_nll_crossfit <- function(data, l, r, l_end = l, r_end = r,
                               nfolds = 5,
                               fold_type = c("op", "blk", "stable_blk"),
                               op_size = 1,
                               buffer_lag = 0,
                               fold_stable_const = 1,
                               loss_output_types = "recv",
                               save_model = FALSE, is_virtual_run = FALSE,
                               bandwidth_vec = NULL,
                               var_dim = NULL,
                               kernel = "gaussian",
                               kernel_args = list(),
                               distance_power = 2L) {
  distance_power <- .kde_validate_distance_power(distance_power)
  fold_type <- match.arg(fold_type)
  external_train_n <- attr(
    data, ".reliever_external_train_n", exact = TRUE
  )
  bandwidth_data <- .kde_bandwidth_reference(
    data, bandwidth_vec, distance_power
  )
  bandwidth_vec <- .kde_resolve_bandwidth_vec(
    bandwidth_vec, var_dim, bandwidth_data, kernel = kernel,
    kernel_args = kernel_args
  )
  kernel_settings <- .kde_kernel_settings(kernel, kernel_args)
  if (!is_virtual_run && !is.null(external_train_n)) {
    return(.kde_nll_crossfit_external(
      data = data,
      l = l,
      r = r,
      l_end = l_end,
      r_end = r_end,
      bandwidth_vec = bandwidth_vec,
      var_dim = var_dim,
      kernel = kernel_settings$kernel,
      kernel_args = kernel_settings$kernel_args,
      distance_power = distance_power,
      nfolds = nfolds,
      fold_type = fold_type,
      op_size = op_size,
      buffer_lag = buffer_lag,
      fold_stable_const = fold_stable_const,
      loss_output_types = loss_output_types
    ))
  }
  out <- reg_fun_crossfit_template(
    data = data,
    l = l,
    r = r,
    l_end = l_end,
    r_end = r_end,
    loss_fun = .loss_kde_nll_crossfit,
    hyper_set = bandwidth_vec,
    hyper_name = "bandwidth",
    loss_output_types = loss_output_types,
    nfolds = nfolds,
    fold_type = fold_type,
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = save_model,
    is_virtual_run = is_virtual_run,
    var_dim = var_dim,
    kernel = kernel_settings$kernel,
    kernel_args = kernel_settings$kernel_args,
    distance_power = distance_power
  )
  if (is_virtual_run) {
    missing_kind <- is.na(out$loss_output_meta$loss_kind)
    out$loss_output_meta$loss_kind[missing_kind] <-
      "negative_log_likelihood"
  }
  out
}

# Probability-classifier crossfit reg_fun -------------------------------------

.validate_classifier_prior_control <- function(control, engine) {
  unsafe <- if (identical(engine, "ranger")) {
    c("class.weights", "case.weights", "inbag")
  } else {
    "weights"
  }
  supplied <- intersect(names(control), unsafe)
  class_specific_sampling <- identical(engine, "ranger") &&
    "sample.fraction" %in% names(control) &&
    length(control$sample.fraction) > 1L
  if (length(supplied) > 0L || class_specific_sampling) {
    stop(
      engine,
      " class or case reweighting is not supported: classifier ",
      "probabilities must use the empirical training class prior. Use ",
      "reg_fun_clf_crossfit_template() with probabilities recalibrated to ",
      "that prior for a weighted classifier.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.validate_classifier_reserved_args <- function(control, engine, reserved) {
  supplied <- intersect(names(control), reserved)
  if (length(supplied) > 0L) {
    stop(
      engine, " arguments ", paste(supplied, collapse = ", "),
      " are managed by the crossfit wrapper and cannot be supplied.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.classifier_hyper_settings <- function(hyper_set, engine) {
  if (is.data.frame(hyper_set) || is.matrix(hyper_set)) {
    setting_names <- if (is.data.frame(hyper_set)) {
      names(hyper_set)
    } else {
      colnames(hyper_set)
    }
    if (nrow(hyper_set) < 1L || ncol(hyper_set) < 1L ||
        is.null(setting_names) || anyNA(setting_names) ||
        any(!nzchar(setting_names)) || anyDuplicated(setting_names)) {
      stop(
        engine,
        " hyper_set must have at least one row and uniquely named columns.",
        call. = FALSE
      )
    }
    return(lapply(seq_len(nrow(hyper_set)), function(i) {
      .crossfit_hyper_as_list(hyper_set[i, , drop = FALSE])
    }))
  }
  if (!is.list(hyper_set) || length(hyper_set) < 1L ||
      !all(vapply(hyper_set, is.list, logical(1L)))) {
    stop(
      engine,
      " hyper_set must be a data frame or matrix with named columns, ",
      "or a list of named argument lists.",
      call. = FALSE
    )
  }
  for (setting in hyper_set) {
    setting_names <- names(setting)
    if (length(setting) > 0L &&
        (is.null(setting_names) || anyNA(setting_names) ||
         any(!nzchar(setting_names)) || anyDuplicated(setting_names))) {
      stop(
        engine,
        " hyper_set list elements must be named argument lists.",
        call. = FALSE
      )
    }
  }
  hyper_set
}

.validate_classifier_configuration <- function(hyper_set, fixed_args,
                                               engine, fixed_name,
                                               reserved) {
  settings <- .classifier_hyper_settings(hyper_set, engine)
  if (!is.list(fixed_args)) {
    stop(fixed_name, " must be a named list.", call. = FALSE)
  }
  fixed_names <- names(fixed_args)
  if (length(fixed_args) > 0L &&
      (is.null(fixed_names) || anyNA(fixed_names) ||
       any(!nzchar(fixed_names)) || anyDuplicated(fixed_names))) {
    stop(fixed_name, " must be a named list with unique names.",
         call. = FALSE)
  }
  for (setting in settings) {
    overlap <- intersect(names(setting), fixed_names)
    if (length(overlap) > 0L) {
      stop(
        engine, " arguments ", paste(overlap, collapse = ", "),
        " must be supplied in only one of hyper_set and ", fixed_name, ".",
        call. = FALSE
      )
    }
    control <- utils::modifyList(fixed_args, setting)
    .validate_classifier_prior_control(control, engine)
    .validate_classifier_reserved_args(control, engine, reserved)
  }
  invisible(TRUE)
}

.prob_ranger_crossfit <- function(data, y, train_id, eval_id, hyper, hyper_id,
                            ranger_args = list(), ...) {
  ranger_control <- utils::modifyList(ranger_args, hyper)
  .validate_classifier_prior_control(ranger_control, "ranger")
  .validate_classifier_reserved_args(
    ranger_control, "ranger",
    c("formula", "data", "write.forest", "probability")
  )
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("reg_fun_ranger_crossfit() requires the ranger package.", call. = FALSE)
  }
  y_fac <- factor(y, levels = c(0L, 1L))
  dat <- data.frame(Y = y_fac, X = as.matrix(data))
  args <- utils::modifyList(
    list(
      formula = Y ~ .,
      data = dat[train_id, , drop = FALSE],
      write.forest = TRUE,
      probability = TRUE,
      num.threads = 1
    ),
    ranger_control
  )
  model <- do.call(ranger::ranger, args)
  pred_data <- data.frame(X = as.matrix(data[eval_id, , drop = FALSE]))
  pred <- predict(model, data = pred_data)$predictions
  if ("1" %in% colnames(pred)) {
    pred[, "1"]
  } else {
    pred[, ncol(pred)]
  }
}

#' Cross-fitted ranger-classifier interval loss
#'
#' Use \code{ranger::ranger()} as the probability model inside
#' \code{reg_fun_clf_crossfit_template()}. The classifier distinguishes
#' observations inside the fitted interval from observations outside it; its
#' probabilities are converted to density-ratio losses. Most users should call
#' \code{reliever(X, cpd_family = "ranger_crossfit")}.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   predictor features in columns.
#' @inheritParams reg_fun_clf_crossfit_template
#' @param hyper_set A data frame, matrix, or list of candidate ranger settings.
#'   Supply a data frame or matrix
#'   with one row per setting and uniquely named columns, or a list whose
#'   elements are named argument lists.
#' @param ranger_args A named list of fixed \code{ranger::ranger()} arguments
#'   shared by every candidate in \code{hyper_set}. The built-in conversion
#'   rejects class weights, case weights, and class-specific sampling because
#'   it requires probabilities under the empirical training class proportion.
#'   The wrapper manages \code{formula}, \code{data}, \code{write.forest}, and
#'   \code{probability}; do not include them in either argument list. Other
#'   controls, including \code{num.threads}, may be supplied. Supply each
#'   argument in only one of \code{hyper_set} and \code{ranger_args}.
#'   Weighted classifiers require the generic classifier template with
#'   recalibrated probabilities.
#'
#' @return A named list. With \code{is_virtual_run = TRUE}, it contains the
#'   integer scalar \code{n_loss_outputs} and data frame
#'   \code{loss_output_meta}. Otherwise, it contains the requested numeric
#'   observation-by-output \code{loss} matrix and \code{model = NULL}. See
#'   \code{\link{reg_fun_crossfit_template}()} for the full return convention
#'   and column definitions.
#'
#' @encoding UTF-8
#' @references Londschien, M., Bühlmann, P., and Kovács, S. (2023). Random
#'   forests for change point detection. \emph{Journal of Machine Learning
#'   Research}, 24(216), 1--45.
#' @seealso \code{\link{reliever_ranger_crossfit}()} and
#'   \code{\link{reg_fun_clf_crossfit_template}()}.
#' @export
#'
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   set.seed(2026)
#'   n_seg <- 300
#'   x <- rbind(
#'     matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#'   )
#'   out <- suppressWarnings(reg_fun_ranger_crossfit(
#'     data = x, l = 241, r = 660, l_end = 211, r_end = 690,
#'     hyper_set = data.frame(num.trees = 10, min.node.size = 5),
#'     ranger_args = list(seed = 2026),
#'     nfolds = 2
#'   ))
#'   colSums(out$loss)
#' }
reg_fun_ranger_crossfit <- function(data, l, r, l_end = l, r_end = r,
                              nfolds = 5,
                              fold_type = "op",
                              op_size = 1,
                              buffer_lag = 0,
                              fold_stable_const = 1,
                              loss_output_types = "recv",
                              save_model = FALSE, is_virtual_run = FALSE,
                              hyper_set = list(list()),
                              ranger_args = list()) {
  .validate_classifier_configuration(
    hyper_set = hyper_set,
    fixed_args = ranger_args,
    engine = "ranger",
    fixed_name = "ranger_args",
    reserved = c("formula", "data", "write.forest", "probability")
  )
  reg_fun_clf_crossfit_template(
    data = data,
    l = l,
    r = r,
    l_end = l_end,
    r_end = r_end,
    prob_fun = .prob_ranger_crossfit,
    hyper_set = hyper_set,
    hyper_name = "model_setting",
    loss_output_types = loss_output_types,
    nfolds = nfolds,
    fold_type = fold_type,
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = save_model,
    is_virtual_run = is_virtual_run,
    loss = "density_ratio_nll",
    ranger_args = ranger_args
  )
}

.prob_mlp_crossfit <- function(data, y, train_id, eval_id, hyper, hyper_id,
                         nnet_args = list(), ...) {
  nnet_control <- utils::modifyList(
    list(
      size = 8,
      entropy = TRUE,
      trace = FALSE,
      maxit = 100
    ),
    nnet_args
  )
  nnet_control <- utils::modifyList(nnet_control, hyper)
  .validate_classifier_prior_control(nnet_control, "nnet")
  .validate_classifier_reserved_args(nnet_control, "nnet", c("x", "y"))
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("reg_fun_mlp_crossfit() requires the nnet package.", call. = FALSE)
  }
  args <- c(
    list(
      x = as.matrix(data[train_id, , drop = FALSE]),
      y = y[train_id]
    ),
    nnet_control
  )
  model <- do.call(nnet::nnet, args)
  as.numeric(stats::predict(
    model,
    newdata = as.matrix(data[eval_id, , drop = FALSE]),
    type = "raw"
  ))
}

#' Cross-fitted MLP-classifier interval loss
#'
#' Use \code{nnet::nnet()} as the probability model inside
#' \code{reg_fun_clf_crossfit_template()}. The classifier distinguishes
#' observations inside the fitted interval from observations outside it; its
#' probabilities are converted to density-ratio losses. This function requires
#' the optional \code{nnet} package. Most users should call
#' \code{reliever(X, cpd_family = "mlp_crossfit")}.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   predictor features in columns.
#' @inheritParams reg_fun_clf_crossfit_template
#' @param hyper_set A data frame, matrix, or list of candidate MLP settings.
#'   Supply a data frame or matrix with
#'   one row per setting and uniquely named columns, or a list whose elements
#'   are named argument lists.
#' @param nnet_args A named list of fixed \code{nnet::nnet()} arguments shared by
#'   every candidate in \code{hyper_set}. Observation \code{weights} are
#'   rejected because the built-in probability conversion requires the
#'   empirical training class proportion. The wrapper manages \code{x} and
#'   \code{y}; do not include them in either argument list. Weighted models
#'   require the generic classifier template with recalibrated probabilities.
#'   Supply each argument in only one of \code{hyper_set} and
#'   \code{nnet_args}.
#'
#' @return A named list. With \code{is_virtual_run = TRUE}, it contains the
#'   integer scalar \code{n_loss_outputs} and data frame
#'   \code{loss_output_meta}. Otherwise, it contains the requested numeric
#'   observation-by-output \code{loss} matrix and \code{model = NULL}. See
#'   \code{\link{reg_fun_crossfit_template}()} for the full return convention
#'   and column definitions.
#'
#' @seealso \code{\link{reliever_mlp_crossfit}()} and
#'   \code{\link{reg_fun_clf_crossfit_template}()}.
#' @export
#'
#' @examples
#' if (requireNamespace("nnet", quietly = TRUE)) {
#'   set.seed(2026)
#'   n_seg <- 300
#'   x <- rbind(
#'     matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'     matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#'   )
#'   out <- suppressWarnings(reg_fun_mlp_crossfit(
#'     data = x, l = 241, r = 660, l_end = 211, r_end = 690,
#'     hyper_set = data.frame(size = 2, maxit = 30),
#'     nfolds = 2
#'   ))
#'   colSums(out$loss)
#' }
reg_fun_mlp_crossfit <- function(data, l, r, l_end = l, r_end = r,
                           nfolds = 5,
                           fold_type = "op",
                           op_size = 1,
                           buffer_lag = 0,
                           fold_stable_const = 1,
                           loss_output_types = "recv",
                           save_model = FALSE, is_virtual_run = FALSE,
                           hyper_set = data.frame(size = 8),
                           nnet_args = list()) {
  .validate_classifier_configuration(
    hyper_set = hyper_set,
    fixed_args = nnet_args,
    engine = "nnet",
    fixed_name = "nnet_args",
    reserved = c("x", "y")
  )
  reg_fun_clf_crossfit_template(
    data = data,
    l = l,
    r = r,
    l_end = l_end,
    r_end = r_end,
    prob_fun = .prob_mlp_crossfit,
    hyper_set = hyper_set,
    hyper_name = "model_setting",
    loss_output_types = loss_output_types,
    nfolds = nfolds,
    fold_type = fold_type,
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = save_model,
    is_virtual_run = is_virtual_run,
    loss = "density_ratio_nll",
    nnet_args = nnet_args
  )
}
