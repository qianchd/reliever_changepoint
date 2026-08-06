# Result row and data setup helpers -------------------------------------------

.segment_result_run_meta <- function(result, run_ids = NULL) {
  run_meta <- .model_select_metadata(result)
  if (is.null(run_ids)) {
    return(run_meta)
  }
  run_ids <- unique(.reliever_validate_positive_integer_vector(
    run_ids, "run_ids"
  ))
  if (!all(run_ids %in% run_meta$run_id)) {
    stop("run_ids contains values that are not present in result$run_meta.",
         call. = FALSE)
  }
  out <- run_meta[match(run_ids, run_meta$run_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.segment_evaluation_spec <- function(result, data, run_meta, reg_fun,
                                     arg_overrides) {
  stored <- attr(result, ".reliever_loss_spec", exact = TRUE)
  args <- if (is.null(stored$args)) list() else stored$args

  if (is.null(reg_fun)) {
    reg_fun <- stored$reg_fun
  }
  if (!is.function(reg_fun)) {
    stop(
      "reg_fun is unavailable in result; supply a reg_fun for a manually constructed result.",
      call. = FALSE
    )
  }

  if (length(arg_overrides) > 0L) {
    override_names <- names(arg_overrides)
    if (is.null(override_names) || any(!nzchar(override_names))) {
      stop("Additional reg_fun arguments must be named.", call. = FALSE)
    }
    for (name in override_names) {
      args[name] <- arg_overrides[name]
    }
  }

  virtual_info <- .reliever_reg_fun_virtual_info(reg_fun, data, args)
  if (max(run_meta$loss_output_id) > virtual_info$n_loss_outputs) {
    stop(
      "The evaluation reg_fun returns fewer loss outputs than result$run_meta requires.",
      call. = FALSE
    )
  }
  expected_meta <- .reliever_run_meta(
    virtual_info$loss_output_meta, run_meta$loss_output_id
  )
  expected_meta$run_id <- run_meta$run_id
  .reliever_validate_run_meta_structure(
    expected_meta,
    run_meta,
    paste0(
      "The evaluation reg_fun must expose the same run structure and ",
      "hyperparameter path as the fitted result."
    )
  )

  list(
    reg_fun = reg_fun,
    args = args,
    n_loss_outputs = virtual_info$n_loss_outputs
  )
}

.segment_subset_rows <- function(x, idx) {
  if (is.null(dim(x))) {
    x[idx]
  } else {
    x[idx, , drop = FALSE]
  }
}

.segment_default_data_stack_fun <- function(data, eval_data) {
  if (is.vector(data) && is.null(dim(data))) {
    return(c(data, eval_data))
  }
  if (is.data.frame(data)) {
    return(rbind(data, eval_data))
  }
  rbind(as.matrix(data), as.matrix(eval_data))
}

.segment_response_matrix <- function(data, response, expected_predictors,
                                     data_name, response_name) {
  if (is.data.frame(data)) {
    numeric_columns <- vapply(data, is.numeric, logical(1L))
    if (!all(numeric_columns)) {
      stop(data_name, " must contain only numeric columns.", call. = FALSE)
    }
  }
  data <- as.matrix(data)
  if (!is.numeric(data) || length(dim(data)) != 2L) {
    stop(data_name, " must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (is.null(response)) {
    if (ncol(data) == expected_predictors) {
      stop(
        data_name, " contains ", expected_predictors,
        " predictor columns but no response. Supply ", response_name,
        ", or prepend the response as the first column of ", data_name, ".",
        call. = FALSE
      )
    }
    if (ncol(data) != expected_predictors + 1L) {
      stop(
        data_name, " must contain ", expected_predictors + 1L,
        " columns when ", response_name,
        " is omitted: one response followed by ", expected_predictors,
        " predictors.",
        call. = FALSE
      )
    }
    return(data)
  }

  if (!is.numeric(response) || is.matrix(response) ||
      length(response) != nrow(data) || anyNA(response) ||
      any(!is.finite(response))) {
    stop(
      response_name, " must be a finite numeric vector with one value per ",
      data_name, " row.",
      call. = FALSE
    )
  }
  if (ncol(data) != expected_predictors) {
    stop(
      data_name, " must contain exactly ", expected_predictors,
      " predictor columns when ", response_name, " is supplied.",
      call. = FALSE
    )
  }
  cbind(as.numeric(response), data)
}

.segment_prepare_input_pair <- function(result, data, eval_data,
                                        y = NULL, eval_y = NULL) {
  input_spec <- result$settings$input_spec
  response_based <- is.list(input_spec) &&
    identical(input_spec$type, "response_predictor")

  if (response_based) {
    if (xor(is.null(y), is.null(eval_y))) {
      stop(
        "Supply both y and eval_y for separate predictor/response input, ",
        "or omit both and prepend each response to data and eval_data.",
        call. = FALSE
      )
    }
    expected_predictors <- input_spec$n_predictors
    if (!is.numeric(expected_predictors) ||
        length(expected_predictors) != 1L ||
        is.na(expected_predictors) || expected_predictors < 1L ||
        expected_predictors != floor(expected_predictors)) {
      stop("result$settings$input_spec has an invalid n_predictors value.",
           call. = FALSE)
    }
    expected_predictors <- as.integer(expected_predictors)

    train_data <- .segment_response_matrix(
      data, y, expected_predictors, "data", "y"
    )
    evaluation_data <- .segment_response_matrix(
      eval_data, eval_y, expected_predictors, "eval_data", "eval_y"
    )
    family <- result$settings$loss_args$family
    if (!is.null(family)) {
      .reliever_validate_lasso_response(
        train_data[, 1L], family,
        if (is.null(y)) "the response column of data" else "y"
      )
      .reliever_validate_lasso_response(
        evaluation_data[, 1L], family,
        if (is.null(eval_y)) {
          "the response column of eval_data"
        } else {
          "eval_y"
        }
      )
    }
    return(list(data = train_data, eval_data = evaluation_data))
  }

  if (!is.null(y) || !is.null(eval_y)) {
    stop(
      "y and eval_y are supported only for a built-in response-based fit.",
      call. = FALSE
    )
  }

  if (is.list(input_spec) && identical(input_spec$type, "kde_nll")) {
    distance_power <- input_spec$distance_power
    if (!is.numeric(distance_power) || length(distance_power) != 1L ||
        is.na(distance_power) || !distance_power %in% c(1, 2)) {
      stop("result$settings$input_spec has an invalid distance_power value.",
           call. = FALSE)
    }
    if (identical(input_spec$original_form, "raw_data")) {
      train_x <- .kernel_data_matrix(data, "data")
      eval_x <- .kernel_data_matrix(eval_data, "eval_data")
      if (ncol(train_x) != input_spec$n_features ||
          ncol(eval_x) != input_spec$n_features) {
        stop(
          "data and eval_data must have the same number of columns as the ",
          "raw data used for fitting.",
          call. = FALSE
        )
      }
      train_distance <- .kernel_dist_sq(train_x)
      eval_distance <- .kernel_dist_sq(eval_x, train_x)
    } else if (identical(
      input_spec$original_form, "squared_distance"
    )) {
      train_distance <- as.matrix(data)
      eval_distance <- as.matrix(eval_data)
      .kde_validate_dist_sq(train_distance)
      .kde_validate_dist_sq(eval_distance, square = FALSE)
      if (!.kde_is_distance_like(train_distance) ||
          ncol(eval_distance) != nrow(train_distance)) {
        stop(
          "For a precomputed KDE-NLL fit, data must be the training squared-",
          "distance matrix and eval_data must have one row per evaluation ",
          "observation and one column per training observation.",
          call. = FALSE
        )
      }
    } else {
      stop("result$settings$input_spec has an invalid KDE input form.",
           call. = FALSE)
    }
    if (distance_power == 1) {
      train_distance <- sqrt(train_distance)
      eval_distance <- sqrt(eval_distance)
    }
    return(list(data = train_distance, eval_data = eval_distance))
  }

  if (is.list(input_spec) &&
      identical(input_spec$type, "kernel_features") &&
      identical(input_spec$original_form, "raw_data")) {
    train_x <- .kernel_data_matrix(data, "data")
    eval_x <- .kernel_data_matrix(eval_data, "eval_data")
    if (ncol(train_x) != input_spec$n_features ||
        ncol(eval_x) != input_spec$n_features) {
      stop(
        "data and eval_data must have the same number of columns as the ",
        "raw data used for fitting.",
        call. = FALSE
      )
    }
    return(list(
      data = .kernel_l2_matrix(
        train_x,
        kernel = input_spec$kernel,
        bandwidth = input_spec$bandwidth,
        kernel_args = input_spec$kernel_args
      ),
      eval_data = .kernel_l2_matrix(
        eval_x, train_x,
        kernel = input_spec$kernel,
        bandwidth = input_spec$bandwidth,
        kernel_args = input_spec$kernel_args
      )
    ))
  }

  if (is.list(input_spec) &&
      identical(input_spec$type, "kernel_features") &&
      identical(input_spec$original_form, "precomputed_features")) {
    if (is.vector(data) && is.null(dim(data)) &&
        is.vector(eval_data) && is.null(dim(eval_data))) {
      return(list(data = data, eval_data = eval_data))
    }
    data <- .kernel_data_matrix(data, "data")
    eval_data <- .kernel_data_matrix(eval_data, "eval_data")
    if (ncol(data) != input_spec$n_features ||
        ncol(eval_data) != input_spec$n_features) {
      stop(
        "data and eval_data must have the same feature columns as the ",
        "precomputed input used for fitting.",
           call. = FALSE)
    }
    return(list(data = data, eval_data = eval_data))
  }

  list(data = data, eval_data = eval_data)
}

.segment_eval_data <- function(data, eval_data, eval_index, data_stack_fun) {
  n <- .reliever_nobs(data)
  eval_n <- .reliever_nobs(eval_data)
  eval_index_int <- suppressWarnings(as.integer(eval_index))
  if (!is.numeric(eval_index) || length(eval_index) != eval_n ||
      anyNA(eval_index) || anyNA(eval_index_int) ||
      any(eval_index != eval_index_int) ||
      any(eval_index_int < 1L | eval_index_int > n)) {
    stop("eval_index must contain one integer original sample index in [1, n] per eval_data row.",
         call. = FALSE)
  }

  ord <- order(eval_index_int)
  eval_index <- eval_index_int[ord]
  eval_data <- .segment_subset_rows(eval_data, ord)
  stacked_data <- data_stack_fun(data, eval_data)
  stacked_n <- .reliever_nobs(stacked_data)
  if (stacked_n != n + eval_n) {
    stop("data_stack_fun(data, eval_data) must return exactly n + n_eval rows, with training rows first and eval rows after them.",
         call. = FALSE)
  }
  attr(stacked_data, ".reliever_external_train_n") <- as.integer(n)

  list(
    data = stacked_data,
    train_data = data,
    eval_index = eval_index,
    eval_rows = n + seq_len(eval_n),
    n = n
  )
}

.segment_eval_one <- function(cache, reg_fun, eval_setup, n_loss_outputs,
                               left, right, save_model, para_list) {
  key <- paste(left, right, sep = ":")
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  eval_rows <- eval_setup$eval_rows[
    eval_setup$eval_index >= left & eval_setup$eval_index <= right
  ]
  if (length(eval_rows) > 0L) {
    tem <- do.call(
      reg_fun,
      c(
        list(
          data = eval_setup$data,
          l = left,
          r = right,
          l_end = min(eval_rows),
          r_end = max(eval_rows),
          save_model = save_model
        ),
        para_list
      )
    )
    loss_eval <- colSums(.reliever_loss_matrix(tem[["loss"]]))
  } else {
    loss_eval <- rep(0, n_loss_outputs)
    tem <- NULL
    if (save_model) {
      tem <- do.call(
        reg_fun,
        c(
          list(
            data = eval_setup$train_data,
            l = left,
            r = right,
            l_end = left,
            r_end = right,
            save_model = TRUE
          ),
          para_list
        )
      )
    }
  }

  if (length(loss_eval) < n_loss_outputs) {
    stop("reg_fun returned fewer loss outputs than the reliever result requires.",
         call. = FALSE)
  }
  out <- list(
    loss = loss_eval,
    model = if (save_model && !is.null(tem)) tem[["model"]] else NULL
  )
  assign(key, out, envir = cache)
  if (save_model) {
    keys <- if (exists(".__keys__", envir = cache, inherits = FALSE)) {
      get(".__keys__", envir = cache, inherits = FALSE)
    } else {
      character()
    }
    assign(".__keys__", c(keys, key), envir = cache)
  }
  out
}

# Apply segment evaluation over a reliever path --------------------------------

.segment_apply_evaluation <- function(candidates,
                                      data, n, run_meta, reg_fun,
                                      n_loss_outputs,
                                      save_model, eval_data, eval_index,
                                      data_stack_fun, para_list) {
  eval_setup <- .segment_eval_data(
    data = data,
    eval_data = eval_data,
    eval_index = eval_index,
    data_stack_fun = data_stack_fun
  )
  eval_setup$data <- .reliever_prepare_builtin_reg_fun_data(
    eval_setup$data, reg_fun, para_list
  )
  if (.reliever_uses_builtin_data_preparation(reg_fun)) {
    eval_setup$train_data <- .segment_subset_rows(
      eval_setup$data, seq_len(eval_setup$n)
    )
  }
  cache <- new.env(parent = emptyenv())

  eval_loss <- rep(NA_real_, nrow(candidates))
  for (id_candidate in seq_len(nrow(candidates))) {
    loss_output_id <- run_meta$loss_output_id[
      match(candidates$run_id[id_candidate], run_meta$run_id)
    ]
    cpd <- as.integer(candidates$cpd[[id_candidate]])
    bounds <- cbind(c(1L, cpd + 1L), c(cpd, n))
    eval_loss[id_candidate] <- 0
    for (id_segment in seq_len(nrow(bounds))) {
      segment_eval <- .segment_eval_one(
        cache = cache,
        reg_fun = reg_fun,
        eval_setup = eval_setup,
        n_loss_outputs = n_loss_outputs,
        left = bounds[id_segment, 1L],
        right = bounds[id_segment, 2L],
        save_model = save_model,
        para_list = para_list
      )
      eval_loss[id_candidate] <-
        eval_loss[id_candidate] + segment_eval$loss[loss_output_id]
    }
  }

  models <- NULL
  if (save_model) {
    keys <- get(".__keys__", envir = cache, inherits = FALSE)
    models <- lapply(keys, function(key) {
      get(key, envir = cache, inherits = FALSE)$model
    })
    names(models) <- keys
  }
  candidates$eval_loss <- eval_loss
  list(candidates = candidates, models = models)
}

# Public segment evaluation API ------------------------------------------------

#' Evaluate candidate segmentations on external data
#'
#' Compute segment loss on external evaluation data after candidate
#' segmentations have been fitted by \code{reliever()}. This is a
#' lower-level helper used by \code{select_holdout()}; most users doing
#' hold-out model selection can call \code{select_holdout()} directly.
#'
#' For each candidate changepoint set, the sample is split into final segments.
#' Within each segment, the model is fit on the corresponding rows of
#' \code{data}, and loss is evaluated on rows of \code{eval_data} whose
#' \code{eval_index} values fall in that segment.
#'
#' For built-in lasso-crossfit and KDE-NLL-crossfit results, inner tuning uses
#' only each final training segment. The selected model is refitted on that
#' segment and scored on \code{eval_data}; paired cross-fitted and in-sample
#' outputs therefore share the corresponding external score. Other crossfit
#' templates require a model-specific external scorer.
#'
#' If a WBS-family fit supplied \code{wbs_stop_crit}, evaluation is restricted
#' to the segmentations mapped from those thresholds and the no-change
#' baseline. Omit \code{wbs_stop_crit} when the full K path is needed.
#'
#' @param result A fitted \code{reliever_result} object.
#' @param data Training observations in the representation expected by the
#'   fitted loss. For a built-in lasso fit, supply either the predictor matrix
#'   together with \code{y}, or a response-first matrix and omit \code{y}.
#'   A built-in KDE fit made from raw observations accepts those observations
#'   again and reconstructs the required training pairwise matrix. For a
#'   precomputed KDE-NLL fit, supply the training squared-distance matrix; for
#'   a precomputed KDE-L2 fit, supply the training feature or Gram matrix.
#' @param eval_data Independent observations on which segment losses are
#'   evaluated. For a built-in lasso fit, supply either the evaluation
#'   predictor matrix together with \code{eval_y}, or a response-first matrix
#'   and omit \code{eval_y}. A raw-data KDE fit accepts raw evaluation
#'   observations. A precomputed KDE-NLL fit instead requires an
#'   evaluation-by-training squared-distance matrix, and a precomputed KDE-L2
#'   fit requires an evaluation-by-feature or evaluation-by-training matrix.
#'   It cannot be \code{NULL}.
#' @param eval_index Integer time index for each row of \code{eval_data}. This
#'   tells the function which estimated segment contains that row. For one
#'   independent evaluation observation at every original time point, the
#'   default \code{NULL} uses \code{seq_len(n)}, where \eqn{n} is the number
#'   of training observations. If evaluation
#'   observations are available only at times 10, 20, and 30, use
#'   \code{c(10, 20, 30)}. Repeated time indices are allowed.
#' @param y,eval_y Optional training and evaluation responses for a built-in
#'   lasso or lasso-crossfit result. Supply both when \code{data} and
#'   \code{eval_data} contain predictors only. Omit both when their first
#'   columns already contain the responses. Do not supply these arguments for
#'   custom or non-response loss families.
#' @param reg_fun Optional interval-loss function used to refit and score each
#'   segment. Results returned by \code{reliever()}, its model-specific
#'   wrappers, \code{reliever_generic()}, and \code{twostep()} retain the
#'   original function automatically, so ordinary calls should leave this as
#'   \code{NULL}. Supply it only for a manually constructed result or an
#'   intentional compatible override.
#' @param data_stack_fun Optional function that stacks \code{data} and
#'   \code{eval_data} before scoring. The default row-binds vectors, matrices,
#'   and data frames. Custom representations may supply another function; its
#'   result must put all training rows before all evaluation rows.
#' @param save_model Logical scalar. If \code{TRUE}, return the segment model
#'   objects produced by \code{reg_fun}.
#' @param run_ids Optional identifiers from \code{result$run_meta}. The default
#'   evaluates every fitted run. Use this to evaluate only selected model
#'   settings after fitting; \code{run_cpd_ids} is instead the fitting-time
#'   argument used by \code{reliever_generic()} to choose which loss outputs
#'   become runs. Supply only one of \code{run_ids} and \code{run_type}.
#' @param run_type Optional character values matched against
#'   \code{result$run_meta$row_type}, such as \code{"recv"} or
#'   \code{"crossfit_homo_hyper"}. Supply only one of \code{run_type} and
#'   \code{run_ids}.
#' @param ... Named loss-function arguments for a manually constructed result,
#'   additional evaluation-only arguments, or explicit overrides of retained
#'   arguments. An override must preserve the selected loss-output identifiers
#'   and hyperparameter path.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{candidates}}{The original candidate segmentations and their
#'     total external \code{eval_loss}.}
#'     \item{\code{models}}{Optional fitted models, named \code{"left:right"}
#'     by their training segment.}
#'   }
#' @seealso \code{\link{select_holdout}()},
#'   \code{\link{reliever}()}, \code{\link{reliever_generic}()},
#'   \code{\link{reg_fun_lasso_solpath}()}
#' @export
#'
#' @examples
#' \donttest{
#' # Constructing and evaluating the KDE-L2 path takes more than five seconds.
#' set.seed(2026)
#' n_seg <- 300
#' x <- rbind(
#'   matrix(rnorm(n_seg * 5, mean = 0, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = 4, sd = 0.5), n_seg, 5),
#'   matrix(rnorm(n_seg * 5, mean = -4, sd = 0.5), n_seg, 5)
#' )
#' holdout <- x + matrix(rnorm(length(x), sd = 0.2), nrow = nrow(x))
#' bandwidth <- sqrt(10)
#' res <- reliever(
#'   X = x, cpd_family = "kde_l2",
#'   kernel = "gaussian", bandwidth = bandwidth,
#'   cpn_max = 7, dm = 30, cov_rate = 0.7, method = "SN"
#' )
#' eval_res <- evaluate_reliever_segments(
#'   result = res, data = x, eval_data = holdout
#' )
#' eval_res$candidates
#' }
evaluate_reliever_segments <- function(result, data, eval_data,
                                       eval_index = NULL,
                                       y = NULL, eval_y = NULL,
                                       save_model = FALSE,
                                       run_ids = NULL,
                                       reg_fun = NULL,
                                       data_stack_fun = NULL, ...,
                                       run_type = NULL) {
  if (!is.list(result) || is.null(result$cpd_path$candidates)) {
    stop("result must contain cpd_path$candidates.",
         call. = FALSE)
  }
  if (is.null(data_stack_fun)) {
    data_stack_fun <- .segment_default_data_stack_fun
  }
  if (missing(eval_data) || is.null(eval_data)) {
    stop("eval_data is required for hold-out segment evaluation.",
         call. = FALSE)
  }
  prepared_input <- .segment_prepare_input_pair(
    result, data, eval_data, y = y, eval_y = eval_y
  )
  data <- prepared_input$data
  eval_data <- prepared_input$eval_data
  if (!is.function(data_stack_fun)) {
    stop("data_stack_fun must be a function.", call. = FALSE)
  }
  if (!is.logical(save_model) || length(save_model) != 1L ||
      is.na(save_model)) {
    stop("save_model must be TRUE or FALSE.", call. = FALSE)
  }
  n <- .reliever_nobs(data)
  eval_n <- .reliever_nobs(eval_data)
  if (is.null(eval_index)) {
    if (eval_n != n) {
      stop(
        "eval_index is required when eval_data does not have one row for ",
        "every original observation.",
        call. = FALSE
      )
    }
    eval_index <- seq_len(n)
  }
  fitted_n <- result$settings$n
  if (!is.null(fitted_n) && !identical(as.integer(n), as.integer(fitted_n))) {
    stop("data must have the same number of observations as the fitted result.",
         call. = FALSE)
  }
  selected_run_ids <- NULL
  if (!is.null(run_ids) || !is.null(run_type)) {
    selected_run_ids <- .model_select_run_ids(
      result, run_ids = run_ids, run_type = run_type
    )
  }
  run_meta <- .segment_result_run_meta(result, selected_run_ids)
  eval_spec <- .segment_evaluation_spec(
    result = result,
    data = data,
    run_meta = run_meta,
    reg_fun = reg_fun,
    arg_overrides = list(...)
  )
  candidates <- .reliever_active_candidates(result)
  candidates <- candidates[
    candidates$run_id %in% run_meta$run_id, , drop = FALSE
  ]
  if (nrow(candidates) == 0L) {
    stop("run_ids or run_type do not match any reliever result rows.",
         call. = FALSE)
  }

  .segment_apply_evaluation(
    candidates = candidates,
    data = data,
    n = n,
    run_meta = run_meta,
    reg_fun = eval_spec$reg_fun,
    n_loss_outputs = eval_spec$n_loss_outputs,
    save_model = save_model,
    eval_data = eval_data,
    eval_index = eval_index,
    data_stack_fun = data_stack_fun,
    para_list = eval_spec$args
  )
}

#' Select a fitted candidate path by hold-out loss
#'
#' Refit every requested candidate segmentation on \code{data}, score it on
#' independent \code{eval_data}, and return the candidate with the smallest
#' hold-out loss.
#'
#' With \code{K = NULL}, the function selects among the requested K/path
#' candidates; fixed-setting paths jointly select K and that model setting.
#' For a crossfit result, use \code{run_type = "recv"} for interval-adaptive
#' tuning or \code{run_type = "crossfit_homo_hyper"} for one homogeneous setting.
#' Supplying \code{K} fixes the changepoint number.
#'
#' Use \code{run_ids} or \code{run_type} to restrict the paths. Compared paths
#' must use the same observation-level scoring rule, units, and normalization.
#' Built-in external scoring is available for ordinary losses,
#' lasso-crossfit, and KDE-NLL-crossfit. Other crossfit templates require a
#' model-specific external scorer.
#'
#' If a WBS-family fit supplied \code{wbs_stop_crit}, only the segmentations
#' mapped from those thresholds and the no-change baseline are compared.
#'
#' @param result A fitted \code{reliever_result} object.
#' @param data Training observations. For a built-in lasso fit, this may be
#'   either the predictor matrix used in \code{reliever(X, y, ...)} when
#'   \code{y} is supplied here, or a response-first matrix when \code{y} is
#'   omitted. Built-in kernel fits accept the same raw or precomputed
#'   representation described in \code{evaluate_reliever_segments()}.
#' @param eval_data Independent observations used to compare candidate
#'   segmentations. For a built-in lasso fit, use an evaluation predictor
#'   matrix with \code{eval_y}, or a response-first matrix without
#'   \code{eval_y}. Kernel inputs and their rectangular precomputed forms are
#'   described in \code{evaluate_reliever_segments()}.
#' @param eval_index Time index for each row of \code{eval_data}. For one
#'   evaluation observation at every original time point, the default
#'   \code{NULL} uses \code{seq_len(n)}, where \eqn{n} is the number of
#'   training observations. Supply explicit indices for sparse, repeated, or
#'   otherwise non-one-to-one evaluation observations.
#' @param y,eval_y Optional training and evaluation responses for a built-in
#'   lasso or lasso-crossfit fit. Supply both with predictor-only
#'   \code{data}/\code{eval_data}; omit both with response-first matrices.
#'   Their values must satisfy the response family stored in the fit.
#' @param data_stack_fun Function used to combine \code{data} and
#'   \code{eval_data}; see the segment-evaluation help page.
#' @param K Optional fixed changepoint number. When omitted, all requested K
#'   candidates are compared; fixed-setting paths can also select their model
#'   setting.
#' @param run_ids Optional identifiers from \code{result$run_meta}. The default
#'   compares every fitted run; supply a subset to compare only selected model
#'   settings. Supply only one of \code{run_ids} and \code{run_type}.
#' @param run_type Optional character values matched against
#'   \code{result$run_meta$row_type}, such as \code{"recv"} or
#'   \code{"crossfit_homo_hyper"}. This is the preferred way to restrict hold-out
#'   comparison to a statistically defined kind of stored loss path. Supply
#'   only one of \code{run_ids} and \code{run_type}.
#' @param reg_fun Optional compatible loss-function override. The default
#'   reuses the function retained in \code{result}; see
#'   \code{evaluate_reliever_segments()}.
#' @param ... Named loss-function arguments for a manually constructed result,
#'   additional evaluation-only arguments, or explicit overrides of retained
#'   arguments. An override must preserve the selected loss-output identifiers
#'   and hyperparameter path.
#'
#' @return One row containing the selected \code{K_hat}, \code{cpd_hat},
#'   hold-out loss in \code{score}, and \code{hyper_value} when available.
#'   If the fit contains several statistically distinct loss-row types,
#'   \code{row_type} identifies the selected source.
#'   Internal \code{run_id} and \code{candidate_id} columns remain available
#'   for path tracing but are hidden by the default print method; use
#'   \code{print(x, details = TRUE)} to show them.
#'   Use \code{evaluate_reliever_segments()} when losses for every candidate
#'   are needed.
#'
#' @seealso \code{\link{reliever}()},
#'   \code{\link{evaluate_reliever_segments}()}, and
#'   \code{\link{select_by_run}()}.
#' @export
#'
#' @examples
#' \donttest{
#' # Fitting and evaluating the high-dimensional lasso path takes over five seconds.
#' set.seed(2026)
#' n <- 900
#' p <- 100
#' tau <- c(300, 600)
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' data <- dgp_linear_regression(n, p, tau, b0, delta, sig = 1)$data
#' holdout <- dgp_linear_regression(n, p, tau, b0, delta, sig = 1)$data
#' fit <- reliever(
#'   X = data[, -1, drop = FALSE], y = data[, 1],
#'   cpd_family = "lasso",
#'   cpn_max = 7, dm = 30, cov_rate = 0.6, # just for speed in example
#'   method = "SN"
#' )
#' selected <- select_holdout(
#'   result = fit,
#'   data = data[, -1, drop = FALSE], y = data[, 1],
#'   eval_data = holdout[, -1, drop = FALSE], eval_y = holdout[, 1]
#' )
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#' selected_at_k <- select_holdout(
#'   result = fit, data = data,
#'   eval_data = holdout, K = 2
#' )
#' selected_at_k
#' stopifnot(all(selected_at_k$K_hat == 2L))
#' }
select_holdout <- function(result, data, eval_data, eval_index = NULL,
                           y = NULL, eval_y = NULL,
                           K = NULL, run_ids = NULL,
                           reg_fun = NULL, data_stack_fun = NULL,
                           ..., run_type = NULL) {
  if (missing(eval_data) || is.null(eval_data)) {
    stop("eval_data is required for hold-out selection.",
         call. = FALSE)
  }
  selected_run_ids <- .model_select_run_ids(
    result, run_ids = run_ids, run_type = run_type
  )
  evaluated <- evaluate_reliever_segments(
    result = result,
    data = data,
    eval_data = eval_data,
    eval_index = eval_index,
    y = y,
    eval_y = eval_y,
    reg_fun = reg_fun,
    data_stack_fun = data_stack_fun,
    save_model = FALSE,
    run_ids = selected_run_ids,
    ...
  )
  candidates <- evaluated$candidates
  if (!is.null(K)) {
    K <- .reliever_validate_positive_integer(K, "K", allow_zero = TRUE)
    candidates <- candidates[candidates$K == K, , drop = FALSE]
  }
  finite_id <- which(is.finite(candidates$eval_loss))
  if (length(finite_id) == 0L) {
    stop("No finite hold-out candidate is available for selection.",
         call. = FALSE)
  }
  best_id <- finite_id[order(
    candidates$eval_loss[finite_id],
    candidates$K[finite_id],
    candidates$run_id[finite_id],
    candidates$candidate_id[finite_id]
  )[1L]]
  selected <- candidates[best_id, , drop = FALSE]
  out <- data.frame(
    rule = "holdout",
    run_id = as.integer(selected$run_id[1L]),
    candidate_id = as.integer(selected$candidate_id[1L]),
    K_hat = as.integer(selected$K[1L]),
    score = as.numeric(selected$eval_loss[1L]),
    stringsAsFactors = FALSE
  )
  out$cpd_hat <- I(list(as.integer(selected$cpd[[1L]])))
  .model_select_add_metadata(
    out[c("rule", "K_hat", "cpd_hat", "score", "run_id", "candidate_id")],
    result
  )
}
