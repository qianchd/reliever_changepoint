.cv_reliever_default_data_folding_fun <- function(data, train_id, eval_id) {
  list(
    train_data = .segment_subset_rows(data, train_id),
    eval_data = .segment_subset_rows(data, eval_id)
  )
}

.cv_reliever_eval_index <- function(train_id, eval_id) {
  n_train <- length(train_id)
  # This preserves comparisons with original boundaries train_id[t_k].
  pmin.int(n_train, findInterval(eval_id, train_id) + 1L)
}

.cv_reliever_fold_grid_size <- function(dc_grid_size, n_train, n) {
  max(1L, as.integer(ceiling(
    as.numeric(dc_grid_size) * n_train / n
  )))
}

.cv_reliever_reject_args <- function(args, caller, call = NULL) {
  supplied <- names(args)
  if (!is.null(call) && "dc_grid" %in% names(as.list(call)[-1L])) {
    supplied <- c(supplied, "dc_grid")
  }
  unsupported <- intersect(
    unique(supplied), c("dc_grid", "cache_profile", "cpn_crit")
  )
  if (length(unsupported) == 0L) {
    return(invisible(NULL))
  }
  hint <- if ("dc_grid" %in% unsupported) {
    " Use dc_grid_size instead."
  } else {
    ""
  }
  stop(
    caller, "() does not support ",
    paste(unsupported, collapse = ", "), ".",
    hint,
    call. = FALSE
  )
}

.cv_reliever_selector_rows <- function(result, candidates) {
  select_by <- result$cpd_path$select_by
  candidate_key <- paste(candidates$run_id, candidates$candidate_id)

  if (identical(select_by, "K")) {
    selectors <- candidates[c("run_id", "candidate_id", "K")]
    selectors$selector_id <- selectors$K + 1L
    selectors$select_value <- as.numeric(selectors$K)
  } else {
    selectors <- result$cpd_path$selector_map
    if (is.null(selectors)) {
      stop("A non-K path must provide cpd_path$selector_map.", call. = FALSE)
    }
    selectors <- selectors[
      selectors$run_id %in% candidates$run_id,
      c("run_id", "select_value", "candidate_id"),
      drop = FALSE
    ]
    selectors$selector_id <- as.integer(stats::ave(
      seq_len(nrow(selectors)), selectors$run_id, FUN = seq_along
    ))
    candidate_id <- match(
      paste(selectors$run_id, selectors$candidate_id), candidate_key
    )
    if (anyNA(candidate_id)) {
      stop("A selector refers to an unavailable changepoint candidate.",
           call. = FALSE)
    }
    selectors$K <- candidates$K[candidate_id]
  }

  if ("eval_loss" %in% names(candidates)) {
    candidate_id <- match(
      paste(selectors$run_id, selectors$candidate_id), candidate_key
    )
    selectors$eval_loss <- candidates$eval_loss[candidate_id]
  }
  selectors <- selectors[
    order(selectors$run_id, selectors$selector_id),
    c(
      "run_id", "selector_id", "select_value", "candidate_id", "K",
      intersect("eval_loss", names(selectors))
    ),
    drop = FALSE
  ]
  rownames(selectors) <- NULL
  selectors
}

.cv_reliever_loss_table <- function(fold_loss, nfolds) {
  group_id <- interaction(
    fold_loss$loss_output_id, fold_loss$selector_id,
    drop = TRUE, lex.order = TRUE
  )
  rows <- lapply(split(fold_loss, group_id), function(x) {
    if (length(unique(x$select_value)) != 1L) {
      stop("A selector position must have the same value in every fold.",
           call. = FALSE)
    }
    fold_mean <- x$eval_loss / x$n_eval
    data.frame(
      loss_output_id = x$loss_output_id[1L],
      selector_id = x$selector_id[1L],
      select_value = x$select_value[1L],
      cv_loss = sum(x$eval_loss),
      cv_mean = sum(x$eval_loss) / sum(x$n_eval),
      cv_se = if (nrow(x) > 1L) {
        stats::sd(fold_mean) / sqrt(nrow(x))
      } else {
        NA_real_
      },
      n_folds = length(unique(x$fold_id)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$complete <- out$n_folds == nfolds
  out <- out[order(out$loss_output_id, out$selector_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.cv_reliever_outer <- function(data, search_fun, search_args,
                               dm, nfolds, op_size, dc_grid_size,
                               data_folding_fun, data_stack_fun, detail, echo) {
  if (!is.function(data_folding_fun)) {
    stop("data_folding_fun must be a function.", call. = FALSE)
  }
  if (!is.function(data_stack_fun)) {
    stop("data_stack_fun must be a function.", call. = FALSE)
  }
  detail_mode <- .reliever_detail_mode(detail)
  n <- .reliever_nobs(data)
  dm <- .reliever_validate_positive_integer(dm, "dm")
  nfolds <- .reliever_validate_positive_integer(nfolds, "nfolds")
  if (nfolds < 2L || nfolds > n) {
    stop("nfolds must be between 2 and the number of observations.",
         call. = FALSE)
  }
  op_size <- .reliever_validate_positive_integer(op_size, "op_size")
  if (!is.null(dc_grid_size)) {
    dc_grid_size <- .reliever_validate_positive_integer(
      dc_grid_size, "dc_grid_size"
    )
  }
  fold_id <- rep(rep(seq_len(nfolds), each = op_size), length.out = n)
  if (any(tabulate(fold_id, nbins = nfolds) == 0L)) {
    stop("nfolds and op_size must leave every outer fold nonempty.",
         call. = FALSE)
  }

  fold_losses <- vector("list", nfolds)
  fold_dm <- integer(nfolds)
  expected_run_meta <- NULL
  expected_select_by <- NULL
  cache_backend_warning_seen <- FALSE
  run_search <- function(args) {
    withCallingHandlers(
      do.call(search_fun, args),
      reliever_cache_backend_warning = function(w) {
        if (cache_backend_warning_seen) {
          invokeRestart("muffleWarning")
        }
        cache_backend_warning_seen <<- TRUE
      }
    )
  }

  for (fold in seq_len(nfolds)) {
    eval_id <- which(fold_id == fold)
    train_id <- which(fold_id != fold)
    fold_dm[fold] <- max(1L, ceiling(dm * length(train_id) / n))
    fold_search_args <- list(
      data = NULL,
      dm = fold_dm[fold],
      detail = FALSE
    )
    if (!is.null(dc_grid_size)) {
      fold_grid_size <- .cv_reliever_fold_grid_size(
        dc_grid_size = dc_grid_size,
        n_train = length(train_id),
        n = n
      )
      fold_search_args$dc_grid_size <- fold_grid_size
    }
    fold_parts <- data_folding_fun(data, train_id, eval_id)
    if (!is.list(fold_parts) ||
        is.null(fold_parts$train_data) || is.null(fold_parts$eval_data)) {
      stop("data_folding_fun must return train_data and eval_data.", call. = FALSE)
    }
    if (.reliever_nobs(fold_parts$train_data) != length(train_id) ||
        .reliever_nobs(fold_parts$eval_data) != length(eval_id)) {
      stop("data_folding_fun returned data with unexpected row counts.", call. = FALSE)
    }
    if (echo) {
      message(sprintf("Outer CV: fold %d/%d", fold, nfolds))
    }

    fold_search_args$data <- fold_parts$train_data
    fold_fit <- run_search(c(fold_search_args, search_args))
    if (is.null(expected_run_meta)) {
      expected_run_meta <- fold_fit$run_meta
      expected_select_by <- fold_fit$cpd_path$select_by
    } else {
      .reliever_validate_run_meta_structure(
        fold_fit$run_meta,
        expected_run_meta,
        "The search must expose the same run structure in every outer fold.",
        compare_hyper_value = FALSE
      )
      if (!identical(fold_fit$cpd_path$select_by, expected_select_by)) {
        stop("The search must use the same path selector in every outer fold.",
             call. = FALSE)
      }
    }

    evaluated <- evaluate_reliever_segments(
      result = fold_fit,
      data = fold_parts$train_data,
      eval_data = fold_parts$eval_data,
      eval_index = .cv_reliever_eval_index(train_id, eval_id),
      data_stack_fun = data_stack_fun
    )
    fold_selectors <- .cv_reliever_selector_rows(
      fold_fit, evaluated$candidates
    )
    fold_selectors$loss_output_id <- fold_fit$run_meta$loss_output_id[
      match(fold_selectors$run_id, fold_fit$run_meta$run_id)
    ]
    if (anyDuplicated(fold_selectors[c("loss_output_id", "selector_id")])) {
      stop("Each fold must provide one candidate per selector position.",
           call. = FALSE)
    }
    fold_losses[[fold]] <- data.frame(
      fold_id = fold,
      loss_output_id = fold_selectors$loss_output_id,
      selector_id = fold_selectors$selector_id,
      select_value = fold_selectors$select_value,
      K = fold_selectors$K,
      eval_loss = fold_selectors$eval_loss,
      n_eval = length(eval_id),
      stringsAsFactors = FALSE
    )
  }

  fold_loss <- do.call(rbind, fold_losses)
  cv_loss_raw <- .cv_reliever_loss_table(fold_loss, nfolds)
  cv_loss_raw <- cv_loss_raw[cv_loss_raw$complete, , drop = FALSE]
  if (!any(is.finite(cv_loss_raw$cv_loss))) {
    stop("No finite candidate is shared by all folds.", call. = FALSE)
  }

  if (echo) {
    message("Outer CV: fitting changepoint candidates on all observations")
  }
  full_search_args <- list(data = data, dm = dm, detail = detail)
  if (!is.null(dc_grid_size)) {
    full_search_args$dc_grid_size <- dc_grid_size
  }
  full_fit <- run_search(c(full_search_args, search_args))
  .reliever_validate_run_meta_structure(
    full_fit$run_meta,
    expected_run_meta,
    paste0(
      "The search must expose the same run structure in folds and the ",
      "full-data fit."
    ),
    compare_hyper_value = FALSE
  )
  if (!identical(full_fit$cpd_path$select_by, expected_select_by)) {
    stop("The search must use the same path selector in folds and the full-data fit.",
         call. = FALSE)
  }

  full_selectors <- .cv_reliever_selector_rows(
    full_fit, .reliever_active_candidates(full_fit)
  )
  full_selectors$loss_output_id <- full_fit$run_meta$loss_output_id[
    match(full_selectors$run_id, full_fit$run_meta$run_id)
  ]
  selector_key <- paste(
    full_selectors$loss_output_id, full_selectors$selector_id
  )
  full_id <- match(
    paste(cv_loss_raw$loss_output_id, cv_loss_raw$selector_id), selector_key
  )
  available_in_full <- !is.na(full_id)
  cv_loss_raw <- cv_loss_raw[available_in_full, , drop = FALSE]
  full_id <- full_id[available_in_full]
  if (nrow(cv_loss_raw) == 0L ||
      !any(is.finite(cv_loss_raw$cv_loss))) {
    stop("No finite candidate is shared by all folds and the full-data path.",
         call. = FALSE)
  }
  if (any(full_selectors$select_value[full_id] != cv_loss_raw$select_value)) {
    stop("Fold selectors do not match the full-data path.", call. = FALSE)
  }
  cv_loss_raw$K <- full_selectors$K[full_id]
  cv_loss_raw$candidate_id <- full_selectors$candidate_id[full_id]

  eligible <- which(is.finite(cv_loss_raw$cv_loss))
  best <- eligible[order(
    cv_loss_raw$cv_loss[eligible],
    cv_loss_raw$K[eligible],
    match(
      cv_loss_raw$loss_output_id[eligible],
      full_fit$run_meta$loss_output_id
    ),
    cv_loss_raw$selector_id[eligible]
  )[1L]]

  meta_id <- match(
    cv_loss_raw$loss_output_id, full_fit$run_meta$loss_output_id
  )
  cv_loss_internal <- cbind(
    full_fit$run_meta[meta_id, , drop = FALSE],
    cv_loss_raw[c("select_value", "K", "cv_loss", "cv_mean", "cv_se")]
  )
  rownames(cv_loss_internal) <- NULL

  selected_run <- full_selectors$run_id[full_id[best]]
  final_candidates <- full_fit$cpd_path$candidates
  final_id <- which(
    final_candidates$run_id == selected_run &
      final_candidates$candidate_id == cv_loss_raw$candidate_id[best]
  )
  if (length(final_id) != 1L) {
    stop("The selected candidate is unavailable in the full-data fit.",
         call. = FALSE)
  }

  selected_loss <- cv_loss_internal[best, , drop = FALSE]
  summary <- data.frame(
    rule = "outer_cv",
    K_hat = as.integer(selected_loss$K),
    stringsAsFactors = FALSE
  )
  summary$cpd_hat <- I(list(
    as.integer(final_candidates$cpd[[final_id]])
  ))
  summary_columns <- "rule"
  valid_row_types <- if ("row_type" %in% names(cv_loss_internal)) {
    unique(as.character(cv_loss_internal$row_type))
  } else {
    character()
  }
  valid_row_types <- valid_row_types[
    !is.na(valid_row_types) & nzchar(valid_row_types)
  ]
  if (length(valid_row_types) > 1L) {
    summary$row_type <- as.character(selected_loss$row_type)
    summary_columns <- c(summary_columns, "row_type")
  }
  if ("hyper_value" %in% names(selected_loss) &&
      .reliever_column_has_value(selected_loss$hyper_value)) {
    summary$hyper_value <- .reliever_simplify_meta_column(
      selected_loss$hyper_value
    )
    summary_columns <- c(summary_columns, "hyper_value")
  }
  if (!identical(expected_select_by, "K")) {
    summary[[expected_select_by]] <- selected_loss$select_value
    summary_columns <- c(summary_columns, expected_select_by)
  }
  summary$cv_mean <- selected_loss$cv_mean
  summary$cv_se <- selected_loss$cv_se
  summary <- summary[c(
    summary_columns, "K_hat", "cpd_hat", "cv_mean", "cv_se"
  )]
  rownames(summary) <- NULL
  class(summary) <- c("reliever_summary", "data.frame")
  hyper_name <- .reliever_hyper_name(full_fit$run_meta, selected_run)
  if (!is.null(hyper_name) && "hyper_value" %in% names(summary)) {
    attr(summary, "hyper_name") <- hyper_name
  }

  cv_loss <- data.frame(
    K = as.integer(cv_loss_internal$K),
    cv_mean = as.numeric(cv_loss_internal$cv_mean),
    cv_se = as.numeric(cv_loss_internal$cv_se),
    stringsAsFactors = FALSE
  )
  cv_columns <- character()
  if (length(valid_row_types) > 1L) {
    cv_loss$row_type <- as.character(cv_loss_internal$row_type)
    cv_columns <- c(cv_columns, "row_type")
  }
  if ("hyper_value" %in% names(cv_loss_internal) &&
      .reliever_column_has_value(cv_loss_internal$hyper_value)) {
    cv_loss$hyper_value <- .reliever_simplify_meta_column(
      cv_loss_internal$hyper_value
    )
    cv_columns <- c(cv_columns, "hyper_value")
  }
  if (!identical(expected_select_by, "K")) {
    cv_loss[[expected_select_by]] <- cv_loss_internal$select_value
    cv_columns <- c(cv_columns, expected_select_by)
  }
  cv_columns <- c(cv_columns, "K", "cv_mean", "cv_se")
  if (length(unique(cv_loss_internal$run_id)) > 1L) {
    cv_loss$run_id <- as.integer(cv_loss_internal$run_id)
    cv_columns <- c(cv_columns, "run_id")
  }
  cv_loss <- cv_loss[cv_columns]

  search_settings <- full_fit$settings[setdiff(
    names(full_fit$settings),
    c("n", "cpn_crit", "cpn_penalty", "cpn_criterion")
  )]

  out <- list(
    summary = summary,
    cv_loss = cv_loss,
    full_data_fit = full_fit,
    settings = c(
      list(
        n = n,
        selection = "outer_cv",
        nfolds = nfolds,
        op_size = op_size,
        fold_size = as.integer(tabulate(fold_id, nbins = nfolds)),
        fold_dm = fold_dm
      ),
      search_settings
    )
  )
  if (detail_mode != "none") {
    out$diagnostics <- list(
      fold_id = as.integer(fold_id),
      fold_loss = fold_loss
    )
  }
  class(out) <- c("cv_reliever_result", "list")
  out
}

#' Outer cross-validation for generic Reliever models
#'
#' CPSS-style outer cross-validation for a user-supplied interval-loss
#' function. Each outer training fold runs \code{\link{reliever_generic}()};
#' held-out loss selects a model setting and candidate path, and a final
#' full-data fit supplies the reported K and changepoints. This final
#' \code{reliever_generic()} fit is performed automatically and returned in
#' \code{full_data_fit}; users do not need to call \code{reliever_generic()}
#' again after model selection. By default all declared loss outputs are
#' compared; \code{run_cpd_ids} selects a subset.
#'
#' Use \code{data_folding_fun} and \code{data_stack_fun} for custom data representations.
#' Compare only outputs with the same held-out scoring rule, units, and
#' normalization; select them with \code{run_cpd_ids}.
#'
#' Built-in lasso-crossfit and KDE-NLL-crossfit losses provide model-specific
#' external scorers. Other crossfit templates require a model-specific
#' external scorer before they can be used here.
#'
#' @param data,reg_fun Custom data object and interval-loss function; see
#'   \code{\link{reliever_generic}()}.
#' @param cpn_max,dm,cov_rate,method Common path controls; see
#'   \code{\link{cv.reliever}()}.
#' @param nfolds,op_size Common outer-CV controls; see
#'   \code{\link{cv.reliever}()}.
#' @param pen_val,prune_value,M,wbs_seed,wbs_stop_crit Common penalty and
#'   search controls; see \code{\link{cv.reliever}()}.
#' @param cache_backend,owner_key,detail,echo Common output and cache controls;
#'   see \code{\link{cv.reliever}()}.
#' @param dc_grid_size Candidate-grid control; see
#'   \code{\link{cv.reliever}()}.
#' @param run_cpd_ids Optional declared loss-output identifiers to compare.
#'   The default compares every output supplied by \code{reg_fun}.
#' @param data_folding_fun Advanced function called as
#'   \code{data_folding_fun(data, train_id, eval_id)}. It must return a list containing
#'   \code{train_data} and \code{eval_data}, with rows kept in the respective
#'   \code{train_id} and \code{eval_id} orders. The default subsets row-oriented
#'   vectors, matrices, and data frames.
#' @param data_stack_fun Advanced function that combines a fold's training and
#'   held-out data into the object expected by \code{reg_fun}. The default
#'   row-binds vectors, matrices, and data frames.
#' @param ... Additional arguments passed to \code{reg_fun}. Outer CV accepts
#'   \code{dc_grid_size}, but not \code{dc_grid}, \code{cache_profile}, or
#'   \code{cpn_crit}.
#'
#' @return A \code{cv_reliever_result} with:
#'   \describe{
#'     \item{\code{summary}}{The selected K, changepoints, loss output when
#'     applicable, and mean outer-CV loss with its standard error.
#'     \code{summary()} returns this row.}
#'     \item{\code{cv_loss}}{Every compared loss output and candidate path,
#'     with mean held-out loss and standard error.}
#'     \item{\code{full_data_fit}}{The automatically fitted final full-data
#'     \code{\link{reliever_generic}()} result and its complete candidate
#'     paths. No additional external refit is required.}
#'     \item{\code{settings}}{The outer-fold and full-data search settings.}
#'     \item{\code{diagnostics}}{Fold assignments and losses when requested.}
#'   }
#'
#' @references
#' Qian, C., Wang, G., and Zou, C. (2025). Reliever: Relieving the burden of
#' costly model fits for changepoint detection. \emph{Journal of Machine
#' Learning Research}, 26(203), 1--57.
#'
#' Qian, C., Wang, G., Wang, Z., and Zou, C. (2024). Changepoint detection in
#' complex models: Cross-fitting is needed. arXiv:2411.07874.
#'
#' Zou, C., Wang, G., and Li, R. (2020). Consistent selection of
#'   the number of change-points via sample-splitting. \emph{The Annals of
#'   Statistics}, 48(1), 413--439.
#'
#' @seealso \code{\link{cv.reliever}()},
#'   \code{\link{reliever_generic}()}, and
#'   \code{\link{reg_fun_crossfit_template}()}.
#' @export
#'
#' @examples
#' set.seed(2026)
#' x <- c(rnorm(40), rnorm(40, 3), rnorm(40, -3))
#' absolute_location_loss <- function(data, l, r, l_end = l, r_end = r,
#'                                    save_model = FALSE,
#'                                    is_virtual_run = FALSE) {
#'   if (is_virtual_run) {
#'     return(1L)
#'   }
#'   center <- stats::median(data[l:r])
#'   list(
#'     loss = matrix(abs(data[l_end:r_end] - center), ncol = 1L),
#'     model = if (save_model) list(center = center) else NULL
#'   )
#' }
#'
#' fit <- cv.reliever_generic(
#'   data = x,
#'   reg_fun = absolute_location_loss,
#'   cpn_max = 3, dm = 10, cov_rate = 0.6,
#'   method = "SN", nfolds = 3,
#'   dc_grid_size = 5
#' )
#' selected <- summary(fit)
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
cv.reliever_generic <- function(data, reg_fun,
                                cpn_max = 3, dm = 50, cov_rate = 0.8,
                                method = "SN",
                                nfolds = 5, op_size = 1,
                                pen_val = 1, prune_value = 0,
                                M = 100, wbs_seed = NULL,
                                wbs_stop_crit = NULL,
                                run_cpd_ids = NULL,
                                data_folding_fun = NULL,
                                data_stack_fun = NULL,
                                cache_backend = "by_loss_block",
                                owner_key = TRUE,
                                detail = FALSE, echo = FALSE,
                                dc_grid_size = NULL, ...) {
  method <- match.arg(
    method,
    c("SN", "WBS", "WBS_recursive", "SeedBS", "BS", "PELT", "OP")
  )
  cache_backend <- match.arg(
    cache_backend, c("by_loss_block", "by_cost_mat")
  )
  if (is.null(data_folding_fun)) {
    data_folding_fun <- .cv_reliever_default_data_folding_fun
  }
  if (is.null(data_stack_fun)) {
    data_stack_fun <- .segment_default_data_stack_fun
  }
  para_list <- list(...)
  .reliever_reject_renamed_cpn_max(para_list, "cv.reliever_generic")
  .cv_reliever_reject_args(
    para_list, "cv.reliever_generic", call = sys.call()
  )
  .cv_reliever_outer(
    data = data,
    search_fun = reliever_generic,
    search_args = c(
      list(
        reg_fun = reg_fun,
        cpn_max = cpn_max,
        cov_rate = cov_rate,
        method = method,
        cpn_crit = "none",
        pen_val = pen_val,
        prune_value = prune_value,
        M = M,
        wbs_seed = wbs_seed,
        wbs_stop_crit = wbs_stop_crit,
        run_cpd_ids = run_cpd_ids,
        cache_backend = cache_backend,
        owner_key = owner_key,
        echo = FALSE
      ),
      para_list
    ),
    dm = dm,
    nfolds = nfolds,
    op_size = op_size,
    dc_grid_size = dc_grid_size,
    data_folding_fun = data_folding_fun,
    data_stack_fun = data_stack_fun,
    detail = detail,
    echo = echo
  )
}

#' Cross-validated changepoint detection with built-in loss families
#'
#' The sample-efficient cross-fitted version of the CPSS sample-splitting
#' selector of Zou, Wang, and Li (2020). Built-in single-path losses cover
#' changes in the mean, covariance, mean and covariance, linear models,
#' generalized linear models, common exponential-family distributions,
#' fixed-kernel KDE-L2 features, and univariate distributions through NMCD.
#' Lasso additionally supplies a normalized lambda path. Each outer fold
#' reruns the changepoint search on its training rows and scores the fitted
#' segments on held-out rows. Held-out loss selects K and, for lasso, a
#' normalized lambda setting. A final full-data fit supplies the reported
#' changepoints.
#'
#' As part of the same call, \code{cv.reliever()} automatically runs the final
#' \code{\link{reliever}()} fit on all observations. The returned
#' \code{full_data_fit} contains this fit, and \code{summary()} reports its
#' outer-CV-selected K and changepoints. Do not call \code{reliever()} again
#' after \code{cv.reliever()} merely to refit the selected result.
#'
#' Training rows retain their order, and fold-specific changepoints are mapped
#' to the original time axis before held-out observations are scored.
#'
#' Candidates are matched by K, \code{wbs_stop_crit}, or \code{pen_val}, as
#' appropriate. Only candidates available in every training fold and the final
#' full-data path are compared.
#' WBS-family and PELT/OP comparisons also contain an implicit
#' \code{Inf}-valued no-change baseline, which may be selected.
#'
#' The default mean family uses mean-square loss. The \code{"var"} and
#' \code{"meanvar"} families use multivariate Gaussian negative twice
#' log-likelihood. The \code{"lm"} family uses squared residuals,
#' \code{"glm"} uses family deviance residuals, and \code{"em"} uses negative
#' twice log-likelihood for one selected distribution. These are single-path
#' families: outer CV selects K, not a model or distribution family.
#' Fixed-kernel \code{"kde_l2"} evaluates held-out kernel feature rows against
#' segment centers, while \code{"nmcd"} evaluates held-out empirical-CDF loss
#' using fold-specific training references. For KDE-L2, the kernel and
#' bandwidth are fixed before outer CV; this interface selects K rather than
#' tuning the bandwidth.
#'
#' The lasso family evaluates glmnet prediction loss over a normalized
#' \code{lam_set}; a fit with \eqn{m} rows passes
#' \code{lam_set/ sqrt(m)} to glmnet. When \code{lam_set = NULL}, one full-data
#' path, with length requested through \code{nlambda}, is generated and reused
#' in all folds. An explicit \code{lam_set} takes precedence over
#' \code{nlambda}. Lasso fits use no intercept or automatic predictor
#' standardization.
#'
#' For interval-level cross-fitting, use
#' \code{reliever()} with \code{cpd_family = "lasso_crossfit"}, followed by
#' \code{select_by_run(..., run_type = "recv", cpn_crit = "loss")}.
#'
#' @param cpn_max,cov_rate,method Common path controls; see
#'   \code{\link{reliever}()}.
#' @param pen_val,prune_value,M,wbs_seed,wbs_stop_crit Common penalty and
#'   search controls; see \code{\link{reliever}()}.
#' @param cache_backend,owner_key,echo Common cache and output controls; see
#'   \code{\link{reliever}()}.
#' @param X For \code{"mean"}, \code{"var"}, \code{"meanvar"}, \code{"em"},
#'   and \code{"nmcd"}, a numeric vector or matrix with observations in rows;
#'   NMCD requires univariate input. For \code{"kde_l2"}, a precomputed
#'   kernel-feature matrix, or raw observations together with a fixed
#'   \code{kernel} and \code{bandwidth}. For
#'   \code{"lm"}, \code{"glm"}, and \code{"lasso"}, a predictor matrix when
#'   \code{y} is supplied, or a matrix whose leading column or columns contain
#'   the response when \code{y = NULL}. Lasso requires at least two predictor
#'   columns.
#' @param y Optional response for \code{"lm"}, \code{"glm"}, and
#'   \code{"lasso"}. A binomial GLM may use a two-column matrix of successes and
#'   failures. If omitted, the first column of \code{X} is used as the response
#'   by default. It is ignored with a warning for non-regression families.
#' @param cpd_family Built-in outer-CV loss family: \code{"mean"},
#'   \code{"var"}, \code{"meanvar"}, \code{"lm"}, \code{"glm"},
#'   \code{"em"}, \code{"lasso"}, \code{"kde_l2"}, or \code{"nmcd"}.
#' @param dm Minimum segment length for the full-data fit. It is rescaled
#'   proportionally inside each training fold.
#' @param nfolds Number of global outer folds. The sample size need not be
#'   divisible by \code{nfolds}; remainder observations are distributed across
#'   folds.
#' @param op_size Number of consecutive observations assigned to one outer fold
#'   before assignment cycles to the next fold.
#' @param ratio Computational update threshold in \code{(0, 1]} used only for
#'   ungridded mean loss. It affects speed, not the criterion. An explicitly
#'   supplied value is ignored with a warning for other cases.
#' @param detail Detail level. \code{TRUE} or \code{"cache"} retains fold
#'   assignments, losses, and final-fit cache information. Use \code{FALSE} or
#'   \code{"none"} for a compact result.
#' @param ... Model-specific arguments passed to the selected \code{reg_fun};
#'   see that function's help page.
#' @param dc_grid_size Preferred regular-grid spacing for outer CV. The
#'   full-data grid uses this spacing and ends at \eqn{n}; each training-fold
#'   spacing is scaled in proportion to its sample size. The default
#'   \code{NULL} searches every boundary. See \code{\link{reliever}()} for the
#'   candidate-grid motivation and reference.
#'
#' @return An object of class \code{cv_reliever_result} with:
#'   \describe{
#'     \item{\code{summary}}{The selected K, changepoints, model setting when
#'     applicable, and mean CV loss with its standard error. The
#'     \code{summary()} method returns this row.}
#'     \item{\code{cv_loss}}{All compared candidates with their mean held-out
#'     loss and standard error. Use \code{plot()} to visualize this table.}
#'     \item{\code{full_data_fit}}{The automatically fitted final full-data
#'     \code{\link{reliever}()} result and its candidate path. No additional
#'     external \code{reliever()} call is required.}
#'     \item{\code{settings}}{Fold and full-data search settings.}
#'     \item{\code{diagnostics}}{Fold assignments and losses when requested.}
#'   }
#' @references
#' Qian, C., Wang, G., and Zou, C. (2025). Reliever: Relieving the burden of
#' costly model fits for changepoint detection. \emph{Journal of Machine
#' Learning Research}, 26(203), 1--57.
#'
#' Qian, C., Wang, G., Wang, Z., and Zou, C. (2024). Changepoint detection in
#' complex models: Cross-fitting is needed. arXiv:2411.07874.
#'
#' Zou, C., Wang, G., and Li, R. (2020). Consistent selection of
#'   the number of change-points via sample-splitting. \emph{The Annals of
#'   Statistics}, 48(1), 413--439.
#' @seealso \code{\link{reliever}()}, \code{\link{cv.reliever_generic}()},
#'   and \code{\link{select_by_run}()}.
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
#' fit <- cv.reliever(
#'   X = x,
#'   cpn_max = 5, dm = 30, cov_rate = 0.7,
#'   method = "SN", nfolds = 5,
#'   dc_grid_size = 10
#' )
#' selected <- summary(fit)
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#' # fit$full_data_fit is already the final reliever fit on all observations.
#' plot(fit)
#'
#' # CPSS outer CV also selects K for NMCD and fixed-kernel KDE-L2.
#' set.seed(2026)
#' n_nmcd <- 90
#' x_nmcd <- c(
#'   rnorm(n_nmcd, -4, 0.35),
#'   rnorm(n_nmcd, 4, 0.35),
#'   rnorm(n_nmcd, 0, 0.35)
#' )
#' fit_nmcd <- cv.reliever(
#'   X = x_nmcd, cpd_family = "nmcd",
#'   cpn_max = 4, dm = 18, cov_rate = 0.7,
#'   method = "SN", nfolds = 5, dc_grid_size = 10
#' )
#' stopifnot(identical(summary(fit_nmcd)$K_hat, 2L))
#'
#' set.seed(2026)
#' n_kde <- 50
#' x_kde <- c(
#'   rnorm(n_kde, -3, 0.45),
#'   rnorm(n_kde, 3, 0.45),
#'   rnorm(n_kde, 0, 0.45)
#' )
#' fit_kde_l2 <- cv.reliever(
#'   X = x_kde, cpd_family = "kde_l2",
#'   kernel = "gaussian", bandwidth = 0.8,
#'   cpn_max = 4, dm = 10, cov_rate = 0.6,
#'   method = "SN", nfolds = 5, dc_grid_size = 5
#' )
#' stopifnot(identical(summary(fit_kde_l2)$K_hat, 2L))
#'
#' \donttest{
#' # This high-dimensional lasso example takes more than five seconds.
#' # Lasso outer CV jointly selects one automatic-path setting and K.
#' p <- 20
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' set.seed(2026)
#' reg_data <- dgp_linear_regression(
#'   n = 450, p = p, tau = c(150, 300),
#'   b0 = b0, delta = delta, sig = 1
#' )$data
#' fit_lasso <- cv.reliever(
#'   X = reg_data[, -1, drop = FALSE], y = reg_data[, 1],
#'   cpd_family = "lasso",
#'   cpn_max = 5, dm = 15, cov_rate = 0.6, method = "SN",
#'   nfolds = 3, nlambda = 20
#' )
#' plot(fit_lasso)
#' selected_lasso <- summary(fit_lasso)
#' selected_lasso
#' stopifnot(identical(selected_lasso$K_hat, 2L))
#' }
#'
#' # The other single-path parametric families use the same outer-CV selector.
#' set.seed(2027)
#' n <- 100
#' z <- c(rnorm(50, sd = 0.5), rnorm(50, sd = 2))
#' fit_var <- cv.reliever(
#'   X = z, cpd_family = "var", cpn_max = 2, dm = 8,
#'   cov_rate = 0.6, method = "BS", nfolds = 2
#' )
#' plot(fit_var)
#' fit_meanvar <- cv.reliever(
#'   X = z, cpd_family = "meanvar", cpn_max = 2, dm = 8,
#'   cov_rate = 0.6, method = "BS", nfolds = 2
#' )
#' plot(fit_meanvar)
#'
#' predictor <- rnorm(n)
#' coefficient <- rep(c(2, -2), each = 50)
#' response <- coefficient * predictor + rnorm(n, sd = 0.25)
#' fit_lm <- cv.reliever(
#'   X = matrix(predictor, ncol = 1), y = response, cpd_family = "lm",
#'   cpn_max = 2, dm = 8, cov_rate = 0.6, method = "BS", nfolds = 2
#' )
#' plot(fit_lm)
#'
#' trials <- rep(20L, n)
#' success <- rbinom(n, trials, plogis(2.5 * coefficient * predictor))
#' fit_glm <- cv.reliever(
#'   X = matrix(predictor, ncol = 1),
#'   y = cbind(success, trials - success), cpd_family = "glm",
#'   family = binomial(), cpn_max = 2, dm = 8,
#'   cov_rate = 0.6, method = "BS", nfolds = 2
#' )
#' plot(fit_glm)
#'
#' exponential <- c(rexp(50, 1), rexp(50, 5))
#' fit_em <- cv.reliever(
#'   X = exponential, cpd_family = "em", family = "exp",
#'   cpn_max = 2, dm = 8, cov_rate = 0.6, method = "BS", nfolds = 2
#' )
#' plot(fit_em)
#'
#' selected_k <- vapply(
#'   list(fit_var, fit_meanvar, fit_lm, fit_glm, fit_em),
#'   function(object) object$summary$K_hat,
#'   integer(1)
#' )
#' selected_k
#' stopifnot(identical(selected_k, rep(1L, 5L)))
cv.reliever <- function(X, y = NULL, cpd_family = "mean",
                        cpn_max = 3, dm = 50, cov_rate = 0.8,
                        method = "SN", nfolds = 5, op_size = 1,
                        pen_val = 1, prune_value = 0,
                        M = 100, wbs_seed = NULL,
                        wbs_stop_crit = NULL,
                        ratio = 0.9,
                        cache_backend = "by_loss_block",
                        owner_key = TRUE,
                        detail = FALSE, echo = FALSE,
                        dc_grid_size = NULL, ...) {
  cpd_family <- .reliever_match_cpd_family(
    cpd_family,
    c(
      "mean", "var", "meanvar", "lm", "glm", "em", "lasso",
      "kde_l2", "nmcd"
    ),
    "cv.reliever"
  )
  method <- match.arg(
    method,
    c("SN", "WBS", "WBS_recursive", "SeedBS", "BS", "PELT", "OP")
  )
  cache_backend <- match.arg(
    cache_backend, c("by_loss_block", "by_cost_mat")
  )
  family_args <- list(...)
  .reliever_reject_renamed_cpn_max(family_args, "cv.reliever")
  .cv_reliever_reject_args(family_args, "cv.reliever", call = sys.call())

  input_spec <- if (identical(cpd_family, "lasso")) {
    .reliever_response_input_spec(X, y)
  } else {
    NULL
  }
  if (cpd_family %in% c("lm", "glm")) {
    prepared <- .reliever_prepare_parametric_regression(
      X, y, cpd_family, family_args
    )
    data <- prepared$data
    input_spec <- prepared$input_spec
    family_args <- prepared$family_args
  } else {
    data <- .reliever_prepare_input(
      X, y, cpd_family, response_families = "lasso"
    )
  }
  if (identical(cpd_family, "kde_l2")) {
    prepared <- .kernel_l2_prepare_input(
      data = data,
      kernel = family_args$kernel,
      bandwidth = family_args$bandwidth,
      kernel_args = if (is.null(family_args$kernel_args)) {
        list()
      } else {
        family_args$kernel_args
      }
    )
    data <- prepared$data
    input_spec <- prepared$input_spec
    family_args$kernel <- NULL
    family_args$bandwidth <- NULL
    family_args$kernel_args <- NULL
  }
  common_args <- list(
    data = data,
    cpn_max = cpn_max,
    dm = dm,
    cov_rate = cov_rate,
    method = method,
    nfolds = nfolds,
    op_size = op_size,
    pen_val = pen_val,
    prune_value = prune_value,
    M = M,
    wbs_seed = wbs_seed,
    wbs_stop_crit = wbs_stop_crit,
    cache_backend = cache_backend,
    owner_key = owner_key,
    detail = detail,
    echo = echo,
    dc_grid_size = dc_grid_size
  )
  if (!missing(ratio) &&
      (!identical(cpd_family, "mean") || !is.null(dc_grid_size))) {
    warning(
      "ratio is ignored unless cpd_family = \"mean\" and dc_grid_size = NULL.",
      call. = FALSE
    )
  }

  if (identical(cpd_family, "mean")) {
    if (length(family_args) > 0L) {
      supplied <- names(family_args)
      supplied[is.na(supplied) | !nzchar(supplied)] <- "<unnamed>"
      stop(
        "Unused arguments for cpd_family = \"mean\": ",
        paste(supplied, collapse = ", "),
        call. = FALSE
      )
    }
    if (!is.null(dc_grid_size)) {
      out <- do.call(
        cv.reliever_generic,
        c(common_args, list(reg_fun = reg_fun_mean))
      )
    } else {
      out <- .cv_reliever_outer(
        data = data,
        search_fun = reliever_mean,
        search_args = list(
          cpn_max = cpn_max,
          cov_rate = cov_rate,
          method = method,
          pen_val = pen_val,
          prune_value = prune_value,
          M = M,
          wbs_seed = wbs_seed,
          wbs_stop_crit = wbs_stop_crit,
          cpn_crit = "none",
          cache_backend = cache_backend,
          owner_key = owner_key,
          echo = FALSE,
          ratio = ratio
        ),
        dm = dm,
        nfolds = nfolds,
        op_size = op_size,
        dc_grid_size = NULL,
        data_folding_fun = .cv_reliever_default_data_folding_fun,
        data_stack_fun = .segment_default_data_stack_fun,
        detail = detail,
        echo = echo
      )
    }
  } else if (identical(cpd_family, "lasso")) {
    family <- if (is.null(family_args$family)) {
      "gaussian"
    } else {
      family_args$family
    }
    thresh <- if (is.null(family_args$thresh)) {
      1e-7
    } else {
      family_args$thresh
    }
    nlambda <- if (is.null(family_args$nlambda)) {
      30L
    } else {
      family_args$nlambda
    }
    family_args$nlambda <- NULL
    family_args$lam_set <- .reliever_resolve_lam_set(
      data, family_args$lam_set, family, thresh, nlambda = nlambda
    )
    out <- do.call(
      cv.reliever_generic,
      c(common_args, list(reg_fun = reg_fun_lasso_solpath), family_args)
    )
  } else {
    reg_fun <- switch(
      cpd_family,
      var = reg_fun_var,
      meanvar = reg_fun_meanvar,
      lm = reg_fun_lm,
      glm = reg_fun_glm,
      em = reg_fun_em,
      kde_l2 = reg_fun_kde_l2,
      nmcd = reg_fun_nmcd
    )
    if (identical(cpd_family, "var") && is.null(family_args$mu)) {
      family_args$mu <- .reliever_stable_col_means(
        .reliever_parametric_matrix(data)
      )
    }
    if (cpd_family %in% c("glm", "em") &&
        is.null(family_args$family)) {
      stop(
        "family is required for cpd_family = \"", cpd_family, "\".",
        call. = FALSE
      )
    }
    out <- do.call(
      cv.reliever_generic,
      c(common_args, list(reg_fun = reg_fun), family_args)
    )
  }

  out$settings <- c(list(cpd_family = cpd_family), out$settings)
  if (!is.null(input_spec)) {
    out <- .reliever_set_input_spec(out, input_spec)
    out$full_data_fit <- .reliever_set_input_spec(
      out$full_data_fit, input_spec
    )
  }
  out
}

#' Plot held-out losses from cross-validated Reliever
#'
#' Plot mean held-out loss against K for every fitted model setting. The point
#' with the smallest held-out loss among the plotted rows is circled. With the
#' default, unfiltered K path this is the outer-CV selection. For solution-path
#' models,
#' \code{x_axis = "hyperparameter"} shows held-out loss across numeric or
#' categorical tuning values; categorical values follow their declared run
#' order. Supplying \code{K} holds the changepoint number fixed while checking
#' whether a lambda or bandwidth grid spans its minimum. For a numeric grid,
#' an endpoint minimum triggers a suggestion to consider a wider grid.
#'
#' @param x A result returned by \code{cv.reliever()} or
#'   \code{cv.reliever_generic()}.
#' @param x_axis Horizontal-axis choice:
#'   \itemize{
#'     \item \code{"K"} for changepoint number;
#'     \item \code{"hyperparameter"} for values from \code{x$cv_loss}; or
#'     \item \code{"search_value"} for WBS, PELT, or OP controls.
#'   }
#' @param K Optional non-negative changepoint number used only for a fixed-K
#'   hyperparameter plot. With \code{K = NULL}, the best K is selected within
#'   each hyperparameter before plotting.
#' @param run_ids Optional positive integer run identifiers to display. For a
#'   crossfit loss, a hyperparameter plot defaults to its
#'   \code{"crossfit_homo_hyper"} runs. Supply only one of \code{run_ids} and
#'   \code{run_type}.
#' @param run_type Optional character values matched against the
#'   \code{row_type} metadata of the full-data fit. Supply only one of
#'   \code{run_type} and \code{run_ids}.
#' @param show_se Draw one-standard-error bars when available.
#' @param col,lty,pch Graphical parameters for curves and points.
#' @param xlab,ylab,main Optional axis labels and title.
#' @param ... Additional arguments passed to \code{graphics::plot.default()}.
#'
#' @return The rows of \code{x$cv_loss} used for plotting, invisibly, with a
#'   logical \code{selected} column identifying the minimum among those rows.
#'   When the full-data run metadata supplies \code{hyper_name}, tuning-path
#'   plots use it instead of the generic word "hyperparameter".
#' @seealso \code{\link{summary.cv_reliever_result}()},
#'   \code{\link{plot.reliever_result}()},
#'   \code{\link{cv.reliever}()}
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
#' fit <- cv.reliever(
#'   X = x, cpn_max = 20, dm = 30, cov_rate = 0.8,
#'   method = "WBS", M = 200,
#'   wbs_stop_crit = c(5, 10, 15, 20), nfolds = 3
#' )
#' stopifnot(identical(summary(fit)$K_hat, 2L))
#' plot(fit, x_axis = "search_value")
plot.cv_reliever_result <- function(x, x_axis = c("K", "hyperparameter",
                                                  "search_value"),
                                    K = NULL, run_ids = NULL, show_se = TRUE,
                                    col = NULL, lty = 1, pch = 19,
                                    xlab = NULL, ylab = NULL, main = NULL,
                                    ..., run_type = NULL) {
  x_axis <- match.arg(x_axis)
  if (!is.logical(show_se) || length(show_se) != 1L || is.na(show_se)) {
    stop("show_se must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(K)) {
    if (!is.numeric(K) || length(K) != 1L || is.na(K) ||
        K < 0 || K != floor(K)) {
      stop("K must be a non-negative integer.", call. = FALSE)
    }
    K <- as.integer(K)
  }
  if (x_axis != "hyperparameter" && !is.null(K)) {
    stop("K is used only when x_axis = \"hyperparameter\".", call. = FALSE)
  }
  if (!is.data.frame(x$cv_loss) || nrow(x$cv_loss) == 0L) {
    stop("x must contain a non-empty cv_loss table.", call. = FALSE)
  }

  plot_data <- x$cv_loss
  if (!"run_id" %in% names(plot_data)) {
    plot_data$run_id <- 1L
  }
  available <- sort(unique(as.integer(plot_data$run_id)))
  if (is.null(run_ids) && is.null(run_type) &&
      x_axis == "hyperparameter") {
    meta <- x$full_data_fit$run_meta
    if (!is.null(meta) && all(c("run_id", "row_type") %in% names(meta))) {
      fixed_cf <- as.integer(meta$run_id[
        !is.na(meta$row_type) & meta$row_type == "crossfit_homo_hyper"
      ])
      fixed_cf <- intersect(fixed_cf, available)
      if (length(fixed_cf) > 0L) {
        run_ids <- fixed_cf
      }
    }
  }
  if (!is.null(run_type)) {
    run_ids <- .model_select_run_ids(
      x$full_data_fit, run_ids = run_ids, run_type = run_type
    )
  }
  if (!is.null(run_ids)) {
    if (!is.numeric(run_ids) || anyNA(run_ids) ||
        any(run_ids != floor(run_ids)) || any(run_ids < 1L)) {
      stop("run_ids must contain positive integer run identifiers.",
           call. = FALSE)
    }
    run_ids <- unique(as.integer(run_ids))
    if (!all(run_ids %in% available)) {
      stop("run_ids contains values that are not present in x$cv_loss.",
           call. = FALSE)
    }
  } else {
    run_ids <- available
  }
  plot_data <- plot_data[plot_data$run_id %in% run_ids, , drop = FALSE]

  if (x_axis == "search_value") {
    select_by <- x$full_data_fit$cpd_path$select_by
    if (is.null(select_by) || !select_by %in% names(plot_data)) {
      stop(
        "This result is indexed by K and has no WBS/PELT/OP search values.",
        call. = FALSE
      )
    }
    plot_data <- plot_data[is.finite(plot_data[[select_by]]), , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      stop("No finite search values are available to plot.", call. = FALSE)
    }
    x_value <- plot_data[[select_by]]
    group <- plot_data$run_id
    axis_label <- .reliever_search_axis_label(select_by)
    xlab <- if (is.null(xlab)) axis_label else xlab
    main <- if (is.null(main)) "Cross-validated search-value profile" else main
  } else if (x_axis == "hyperparameter") {
    if (!"hyper_value" %in% names(plot_data)) {
      stop("x$cv_loss does not contain hyper_value.", call. = FALSE)
    }
    if (is.null(K)) {
      row_split <- split(seq_len(nrow(plot_data)), plot_data$run_id)
      keep <- vapply(row_split, function(id) {
        finite <- id[is.finite(plot_data$cv_mean[id])]
        if (!length(finite)) return(NA_integer_)
        finite[order(plot_data$cv_mean[finite], plot_data$K[finite])[1L]]
      }, integer(1L))
      plot_data <- plot_data[keep[!is.na(keep)], , drop = FALSE]
    } else {
      plot_data <- plot_data[plot_data$K == K, , drop = FALSE]
      if (nrow(plot_data) == 0L) {
        stop("No requested run contains the selected K.", call. = FALSE)
      }
      if (anyDuplicated(plot_data$run_id)) {
        stop(
          "K does not identify one CV candidate per run; use x_axis = \"K\".",
          call. = FALSE
        )
      }
    }
    hyper_axis <- .reliever_plot_hyper_axis(plot_data$hyper_value)
    if (anyDuplicated(hyper_axis$display)) {
      stop(
        "Hyperparameter values are duplicated; select one run type with ",
        "run_type or run_ids.",
        call. = FALSE
      )
    }
    hyper_order <- order(hyper_axis$position)
    plot_data <- plot_data[hyper_order, , drop = FALSE]
    plot_data$hyper_value <- hyper_axis$display[hyper_order]
    hyper_position <- hyper_axis$position[hyper_order]
    hyper_name <- .reliever_hyper_name(
      x$full_data_fit$run_meta, plot_data$run_id
    )
    hyper_label <- .reliever_hyper_label(hyper_name)
    x_value <- hyper_position
    if (!hyper_axis$is_numeric) {
      attr(x_value, "tick_at") <- hyper_position
      attr(x_value, "tick_label") <- plot_data$hyper_value
    }
    group <- rep.int(1L, nrow(plot_data))
    xlab <- if (is.null(xlab)) hyper_label else xlab
    main <- if (is.null(main)) {
      paste("Cross-validated", tolower(hyper_label), "profile")
    } else {
      main
    }
  } else {
    x_value <- plot_data$K
    group <- plot_data$run_id
    xlab <- if (is.null(xlab)) "Number of changepoints (K)" else xlab
    main <- if (is.null(main)) "Cross-validated changepoint path" else main
  }

  finite <- which(is.finite(plot_data$cv_mean))
  if (!length(finite)) {
    stop("No finite cross-validation loss is available to plot.",
         call. = FALSE)
  }
  best <- finite[order(plot_data$cv_mean[finite], plot_data$K[finite],
                       plot_data$run_id[finite])[1L]]
  plot_data$selected <- seq_len(nrow(plot_data)) == best
  if (x_axis == "hyperparameter" && hyper_axis$is_numeric &&
      !is.null(K) && nrow(plot_data) > 1L &&
      best %in% c(1L, nrow(plot_data))) {
    boundary <- if (best == 1L) "lower" else "upper"
    warning(
      "The fixed-K CV loss is minimized at the ", boundary,
      " endpoint of the ", tolower(hyper_label),
      " grid; consider a wider candidate range.",
      call. = FALSE
    )
  }
  labels <- .reliever_plot_run_labels(x$full_data_fit, unique(group))
  .reliever_draw_path(
    data = plot_data,
    x_value = x_value,
    y_value = plot_data$cv_mean,
    group = group,
    selected = plot_data$selected,
    se = if (show_se) plot_data$cv_se else NULL,
    labels = labels,
    col = col,
    lty = lty,
    pch = pch,
    xlab = xlab,
    ylab = if (is.null(ylab)) "Mean held-out loss" else ylab,
    main = main,
    ...
  )
}

#' Print a cross-validated Reliever result
#'
#' Display the final full-data changepoint estimate selected by outer
#' cross-validation. The printed columns are the same compact columns returned
#' by \code{summary(x)}; fold-level candidate losses remain in
#' \code{x$cv_loss}. The built-in loss family, search method, and number of
#' folds are printed above the result when available.
#'
#' @param x A result returned by \code{cv.reliever()} or
#'   \code{cv.reliever_generic()}.
#' @param max_rows Maximum number of result rows printed. Cross-validated
#'   results normally contain one row; the default is 10.
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#'
#' @return \code{x}, invisibly.
#'   Changepoint locations are printed in full; \code{x$summary} retains them
#'   as a list-column for programmatic use.
#' @export
print.cv_reliever_result <- function(x, max_rows = 10L, ...) {
  max_rows <- .reliever_validate_max_rows(max_rows)
  cat("Cross-validated Reliever result\n")
  context <- character()
  if (!is.null(x$settings$cpd_family)) {
    context <- c(context, paste("Family:", x$settings$cpd_family))
  }
  if (!is.null(x$settings$method)) {
    context <- c(context, paste("Method:", x$settings$method))
  }
  if (!is.null(x$settings$nfolds)) {
    context <- c(context, paste("Folds:", x$settings$nfolds))
  }
  if (length(context)) {
    cat(paste(context, collapse = " | "), "\n", sep = "")
  }
  .reliever_print_table(summary(x), max_rows = max_rows, ...)
  invisible(x)
}

#' Summarize a cross-validated Reliever result
#'
#' Return the final full-data segmentation chosen by outer cross-validation.
#' This is the same compact one-row data frame stored in
#' \code{object$summary}.
#'
#' @param object A result returned by \code{cv.reliever()} or a related
#'   cross-validation wrapper.
#' @param ... Unused.
#'
#' @return A one-row data frame. The compact view contains
#'   \code{rule = "outer_cv"}, \code{K_hat}, the changepoint locations in
#'   list-column \code{cpd_hat}, the selected hyperparameter, WBS stopping
#'   threshold, or PELT/OP penalty when relevant, and \code{cv_mean} with
#'   \code{cv_se}. It includes \code{row_type} when statistically distinct
#'   output types were compared. The complete comparison is in
#'   \code{object$cv_loss}; its optional \code{run_id} maps to
#'   \code{object$full_data_fit$run_meta}. A selected \code{hyper_value} is
#'   printed under its model-specific \code{hyper_name}, such as
#'   \code{lambda}, when that metadata is available.
#' @seealso \code{\link{cv.reliever}()},
#'   \code{\link{cv.reliever_generic}()},
#'   \code{\link{plot.cv_reliever_result}()}
#' @export
summary.cv_reliever_result <- function(object, ...) {
  object$summary
}
