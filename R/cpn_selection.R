# Penalty presets --------------------------------------------------------------

.cpn_penalty <- function(penalty = "sic", n) {
  n <- .reliever_validate_positive_integer(n, "n")

  if (is.numeric(penalty)) {
    if (length(penalty) != 1L || is.na(penalty) || !is.finite(penalty) ||
        penalty < 0) {
      stop("penalty must be a non-negative numeric scalar.", call. = FALSE)
    }
    return(list(
      value = as.numeric(penalty),
      label = as.character(penalty),
      criterion = "loss"
    ))
  }

  if (!is.character(penalty) || length(penalty) != 1L || is.na(penalty)) {
    stop("penalty must be a character preset or a non-negative numeric scalar.",
         call. = FALSE)
  }

  key <- tolower(penalty)
  criterion <- if (startsWith(key, "rss_")) "rss_log_loss" else "loss"
  value <- switch(
    key,
    none = 0,
    loss = 0,
    aic = 2,
    hqc = 2 * log(log(max(n, 3L))),
    sic = log(n),
    rss_aic = 2,
    rss_hqc = 2 * log(log(max(n, 3L))),
    rss_sic = log(n),
    stop(
      "Unknown penalty preset. Supported presets are \"none\", \"loss\", ",
      "\"aic\", \"hqc\", \"sic\", \"rss_aic\", ",
      "\"rss_hqc\", and \"rss_sic\".",
      call. = FALSE
    )
  )

  list(value = as.numeric(value), label = key, criterion = criterion)
}

.reliever_cpn_is_none <- function(cpn_crit) {
  is.character(cpn_crit) &&
    length(cpn_crit) == 1L &&
    !is.na(cpn_crit) &&
    identical(tolower(cpn_crit), "none")
}

.selection_result_n <- function(result, n = NULL) {
  if (!is.null(n)) {
    return(.reliever_validate_positive_integer(n, "n"))
  }
  n_original <- result$settings$n
  if (!is.null(n_original)) {
    return(.reliever_validate_positive_integer(n_original, "result$settings$n"))
  }
  search_n <- result$settings$search_n
  if (!is.null(search_n)) {
    return(.reliever_validate_positive_integer(
      search_n, "result$settings$search_n"
    ))
  }
  stop("n must be supplied when result does not contain result$settings$n.",
       call. = FALSE)
}

# Selection value helpers ------------------------------------------------------

.selection_criterion <- function(loss, K, n, penalty,
                                 criterion = "loss") {
  if (criterion == "loss") {
    return(loss + penalty * K)
  }
  if (criterion != "rss_log_loss") {
    stop("Unknown internal selection criterion.", call. = FALSE)
  }
  if (any(is.finite(loss) & loss < 0)) {
    stop("RSS criteria require non-negative loss values.", call. = FALSE)
  }
  positive_loss <- pmax(loss, .Machine$double.xmin)
  n / 2 * log(positive_loss / n) + penalty * K
}

.model_select_run_ids <- function(result, run_ids = NULL, run_type = NULL) {
  if (!is.list(result) || is.null(result$cpd_path$candidates)) {
    stop("result must contain cpd_path$candidates.", call. = FALSE)
  }
  available <- sort(unique(as.integer(result$cpd_path$candidates$run_id)))
  if (!is.null(run_ids) && !is.null(run_type)) {
    stop("Supply only one of run_ids and run_type.", call. = FALSE)
  }
  if (!is.null(run_type)) {
    if (!is.character(run_type) || length(run_type) == 0L ||
        anyNA(run_type) || any(!nzchar(run_type))) {
      stop("run_type must contain non-empty character values.",
           call. = FALSE)
    }
    meta <- result$run_meta
    if (is.null(meta) || !"row_type" %in% names(meta)) {
      stop("result$run_meta must contain row_type when run_type is supplied.",
           call. = FALSE)
    }
    unmatched <- setdiff(unique(run_type), unique(as.character(meta$row_type)))
    if (length(unmatched) > 0L) {
      available_types <- sort(unique(as.character(meta$row_type)))
      available_types <- available_types[
        !is.na(available_types) & nzchar(available_types)
      ]
      stop(
        "run_type does not match result$run_meta$row_type: ",
        paste(unmatched, collapse = ", "),
        ". Available values are: ",
        paste(available_types, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    run_ids <- unique(as.integer(meta$run_id[meta$row_type %in% run_type]))
  }
  if (is.null(run_ids)) {
    return(available)
  }
  if (!is.numeric(run_ids) || any(is.na(run_ids)) ||
      any(run_ids != floor(run_ids)) || any(run_ids < 1)) {
    stop("run_ids must be positive integer run identifiers.", call. = FALSE)
  }
  run_ids <- unique(as.integer(run_ids))
  if (!all(run_ids %in% available)) {
    stop("run_ids contains values that are not present in result$cpd_path$candidates.",
         call. = FALSE)
  }
  run_ids
}

.reliever_warn_incomparable_loss_kinds <- function(run_meta) {
  if (is.null(run_meta) || !"loss_kind" %in% names(run_meta)) {
    return(invisible(run_meta))
  }
  loss_kind <- unique(as.character(run_meta$loss_kind))
  loss_kind <- loss_kind[!is.na(loss_kind) & nzchar(loss_kind)]
  if (length(loss_kind) > 1L) {
    warning(warningCondition(
      paste0(
        "Selected runs have different loss_kind values (",
        paste(sort(loss_kind), collapse = ", "),
        "); losses may not be comparable. Selection will continue."
      ),
      call = NULL,
      class = "reliever_loss_kind_comparability_warning"
    ))
  }
  invisible(run_meta)
}

.model_select_metadata <- function(result) {
  if (is.null(result$run_meta)) {
    stop("result must contain run_meta.", call. = FALSE)
  }
  meta <- as.data.frame(result$run_meta)
  if (!all(c("run_id", "loss_output_id") %in% names(meta))) {
    stop("result$run_meta must contain run_id and loss_output_id.",
         call. = FALSE)
  }
  meta <- meta[, c("run_id", "loss_output_id",
                   setdiff(names(meta), c("run_id", "loss_output_id"))),
               drop = FALSE]
  rownames(meta) <- NULL
  meta
}

.reliever_active_candidates <- function(result) {
  candidates <- result$cpd_path$candidates
  selector_map <- result$cpd_path$selector_map
  if (is.null(selector_map)) {
    return(candidates)
  }
  active_key <- unique(paste(selector_map$run_id, selector_map$candidate_id))
  candidate_key <- paste(candidates$run_id, candidates$candidate_id)
  candidates[candidate_key %in% active_key, , drop = FALSE]
}

.model_select_candidates <- function(result, run_ids = NULL,
                                     run_type = NULL, cpn_crit, n = NULL) {
  run_ids <- .model_select_run_ids(result, run_ids, run_type)
  n <- .selection_result_n(result, n)
  if (is.null(cpn_crit)) {
    stop("cpn_crit must be supplied.", call. = FALSE)
  }
  penalty_info <- .cpn_penalty(cpn_crit, n)
  if (identical(penalty_info$label, "none")) {
    stop(
      "cpn_crit = \"none\" does not select a model; use \"loss\" for unpenalized selection.",
      call. = FALSE
    )
  }
  candidates <- .reliever_active_candidates(result)
  candidates <- candidates[candidates$run_id %in% run_ids, , drop = FALSE]
  candidates$score <- as.numeric(.selection_criterion(
    loss = candidates$loss,
    K = candidates$K,
    n = n,
    penalty = penalty_info$value,
    criterion = penalty_info$criterion
  ))
  list(
    candidates = candidates,
    run_ids = run_ids,
    penalty_info = penalty_info
  )
}

.model_select_one <- function(candidates) {
  finite_id <- which(is.finite(candidates$score))
  if (length(finite_id) == 0L) {
    return(NULL)
  }
  best_id <- finite_id[order(
    candidates$score[finite_id],
    candidates$K[finite_id],
    candidates$run_id[finite_id],
    candidates$candidate_id[finite_id]
  )[1L]]
  candidates[best_id, , drop = FALSE]
}

.model_select_one_or_missing <- function(candidates) {
  selected <- .model_select_one(candidates)
  if (!is.null(selected)) {
    return(selected)
  }
  selected <- candidates[1L, , drop = FALSE]
  selected$candidate_id <- NA_integer_
  selected$K <- NA_integer_
  selected$loss <- NA_real_
  selected$cpd <- I(list(NA_integer_))
  selected$score <- NA_real_
  selected
}

.model_select_format <- function(selected, penalty_info) {
  if (is.null(selected) || nrow(selected) == 0L) {
    return(data.frame())
  }
  out <- data.frame(
    rule = penalty_info$label,
    K_hat = as.integer(selected$K),
    stringsAsFactors = FALSE
  )
  out$cpd_hat <- I(lapply(selected$cpd, as.integer))
  out$score <- as.numeric(selected$score)
  out$run_id <- as.integer(selected$run_id)
  out$candidate_id <- as.integer(selected$candidate_id)
  rownames(out) <- NULL
  out
}

.model_select_add_metadata <- function(selected, result) {
  meta <- .model_select_metadata(result)
  meta_id <- match(selected$run_id, meta$run_id)
  hyper_name <- NULL
  if ("row_type" %in% names(meta)) {
    available_types <- unique(as.character(meta$row_type))
    available_types <- available_types[!is.na(available_types) &
                                         nzchar(available_types)]
    if (length(available_types) > 1L) {
      selected$row_type <- as.character(meta$row_type[meta_id])
    }
  }
  if ("hyper_value" %in% names(meta)) {
    hyper_value <- meta$hyper_value[meta_id]
    if (.reliever_column_has_value(hyper_value)) {
      selected$hyper_value <- .reliever_simplify_meta_column(hyper_value)
      hyper_name <- .reliever_hyper_name(meta, selected$run_id)
    }
  }
  column_order <- c(
    "rule", "row_type", "hyper_value", "K_hat", "cpd_hat", "score",
    "run_id", "candidate_id"
  )
  selected <- selected[intersect(column_order, names(selected))]
  if (!is.null(hyper_name)) {
    attr(selected, "hyper_name") <- hyper_name
  }
  class(selected) <- c("reliever_model_selection", "data.frame")
  selected
}

#' Print a Reliever model-selection result
#'
#' Display the statistical result without the internal path identifiers.
#' When the fitted object contains several statistically distinct
#' \code{row_type} values, that column remains visible so a cross-fitted
#' selection is identifiable as \code{"recv"}, \code{"crossfit_homo_hyper"}, or
#' another declared type. The returned object still contains \code{run_id} and
#' \code{candidate_id}; set \code{details = TRUE} to print them when tracing
#' the selection back to \code{result$run_meta} and
#' \code{result$cpd_path$candidates}.
#'
#' @param x A model-selection result returned by a post-fit selector.
#' @param details Show the internal run and candidate identifiers.
#' @param max_rows Maximum number of selected model settings printed. The
#'   default is 10; use \code{Inf} to print every row.
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#'
#' @return \code{x}, invisibly. Changepoint locations are printed in full;
#'   the returned object retains them as a list-column for programmatic use.
#' @export
print.reliever_model_selection <- function(x, details = FALSE,
                                           max_rows = 10L, ...) {
  max_rows <- .reliever_validate_max_rows(max_rows)
  shown <- if (isTRUE(details)) {
    x
  } else {
    x[setdiff(names(x), c("run_id", "candidate_id"))]
  }
  attr(shown, "hyper_name") <- attr(x, "hyper_name", exact = TRUE)
  cat("Reliever model selection\n")
  .reliever_print_table(shown, max_rows = max_rows, ...)
  invisible(x)
}

# Public post-fit selection ----------------------------------------------------

#' Select K separately within fitted loss paths
#'
#' A run is one stored loss path. This function selects K separately in every
#' requested run and returns one row per run. Use \code{run_type} or
#' \code{run_ids} to choose the paths. Use \code{select_across_runs()} when one
#' model setting and K should be selected jointly across comparable paths.
#' If a WBS-family fit supplied \code{wbs_stop_crit}, only the segmentations
#' mapped from those thresholds and the no-change baseline are compared.
#' Omit \code{wbs_stop_crit} when selection should traverse the full K path.
#'
#' @param result A \code{reliever_result} returned by \code{reliever()},
#'   \code{reliever_generic()}, or a focused model-specific wrapper.
#' @param run_ids Optional exact integer identifiers from
#'   \code{result$run_meta}. Use these when individual paths must be controlled
#'   directly. If neither \code{run_ids} nor \code{run_type} is supplied, K is
#'   selected separately for every run.
#' @param run_type Optional run-type labels, such as \code{"recv"} or
#'   \code{"crossfit_homo_hyper"}. These are matched against the
#'   \code{row_type} metadata field. This is the preferred way to select a
#'   statistically defined group of stored paths. It filters the runs before
#'   candidate values of K are compared; it is not itself a selection
#'   criterion. Supply only one of \code{run_ids} and \code{run_type}.
#' @param cpn_crit Required rule used to compare candidates with different K.
#'   Use \code{"loss"} to choose the smallest stored loss with no K penalty.
#'   Use \code{"aic"}, \code{"hqc"}, or
#'   \code{"sic"} to minimize
#'   \eqn{loss + \gamma K}, where \eqn{\gamma} is \eqn{2},
#'   \eqn{2\log\log(n)}, or \eqn{\log(n)}. These additive rules may also be
#'   applied directly to RSS. The alternatives \code{"rss_aic"},
#'   \code{"rss_hqc"}, and \code{"rss_sic"} instead minimize
#'   \eqn{n\log(RSS/n)/2 + \gamma K}; this log-RSS form may be preferable when
#'   residual variance is unknown or potentially heterogeneous. A non-negative
#'   number supplies \eqn{\gamma} directly in \eqn{loss + \gamma K}. See
#'   \code{\link{reliever}()} for the same rules in the fitting API.
#' @param n Sample size used in HQC/SIC penalties and the RSS transformation.
#'   The default reads the original observation count from
#'   \code{result$settings}; a compressed candidate grid does not change it.
#'   Supply \code{n} only for a manually constructed result that lacks this
#'   information.
#'
#' @return One row per requested run, containing \code{K_hat},
#'   \code{cpd_hat}, the minimized \code{score}, and available model metadata.
#'   Internal \code{run_id} and \code{candidate_id} columns link back to the
#'   fitted paths and are shown by \code{print(x, details = TRUE)}. A requested
#'   run with no finite candidate score is retained, with \code{NA} selection
#'   fields.
#' @seealso \code{\link{summary.reliever_result}()},
#'   \code{\link{select_across_runs}()},
#'   \code{\link{select_holdout}()}
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
#' fit <- reliever(X = x, cpn_max = 7, dm = 30, cov_rate = 0.8, method = "SN")
#' selected <- select_by_run(result = fit, cpn_crit = "rss_sic")
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
select_by_run <- function(result, run_ids = NULL, run_type = NULL, cpn_crit,
                          n = NULL) {
  if (missing(cpn_crit)) {
    stop("cpn_crit must be supplied.", call. = FALSE)
  }
  selection <- .model_select_candidates(
    result, run_ids = run_ids, run_type = run_type,
    cpn_crit = cpn_crit, n = n
  )
  selected <- lapply(selection$run_ids, function(run_id) {
    candidates <- selection$candidates[
      selection$candidates$run_id == run_id, , drop = FALSE
    ]
    .model_select_one_or_missing(candidates)
  })
  selected <- do.call(rbind, selected)
  .model_select_add_metadata(
    .model_select_format(selected, selection$penalty_info), result
  )
}

#' Select one model across comparable fitted loss paths
#'
#' Select one candidate jointly across comparable stored loss paths.
#' Homogeneous-hyperparameter ReCV (\code{recv_homo_hyper}) uses:
#' \preformatted{
#' select_across_runs(
#'   result = fit, run_type = "crossfit_homo_hyper", cpn_crit = "loss"
#' )
#' }
#' This jointly selects one fixed hyperparameter, K, and the changepoints.
#' Compared paths must use the same loss definition and scale.
#'
#' @inheritParams select_by_run
#' @param run_ids Optional exact integer identifiers from
#'   \code{result$run_meta}. Use these for direct path-level control. Supply
#'   either this argument or \code{run_type} so that unrelated losses, such as
#'   in-sample and cross-fitted losses, are not compared accidentally.
#' @param run_type Optional row types to compare. Supply either this argument
#'   or \code{run_ids}.
#'   \itemize{
#'     \item \code{"recv"} identifies interval-adaptive ReCV paths.
#'     \item \code{"crossfit_homo_hyper"}: fixed-hyperparameter cross-fitted
#'     paths.
#'     \item Other types remain available but trigger a comparability warning.
#'     A type containing only one run also suggests \code{select_by_run()} when
#'     only K is being selected.
#'   }
#'   For a custom \code{reg_fun}, advanced users should prefer explicit
#'   \code{run_ids} for deliberately comparable paths.
#' @param cpn_crit Required rule used to compare K and runs. It accepts the
#'   rules documented in \code{\link{select_by_run}()}; \code{"loss"} is
#'   unpenalized. Supply this argument explicitly so the
#'   statistical choice is visible.
#'
#' @section Guidance warnings:
#' Selection continues after each guidance warning:
#' \itemize{
#'   \item Unusual run types signal
#'   \code{"reliever_across_run_comparability_warning"}.
#'   \item A type resolving to one run signals
#'   \code{"reliever_single_run_selection_warning"}.
#'   \item Different declared loss kinds signal
#'   \code{"reliever_loss_kind_comparability_warning"}.
#' }
#' Explicit \code{run_ids} bypass the run-type guidance warnings.
#'
#' @return One row containing the jointly selected \code{K_hat},
#'   \code{cpd_hat}, minimized \code{score}, and available model metadata.
#'   Internal path identifiers are shown by
#'   \code{print(x, details = TRUE)}; see \code{\link{select_by_run}()} for the
#'   common column definitions.
#' @seealso \code{\link{summary.reliever_result}()},
#'   \code{\link{select_by_run}()},
#'   \code{\link{select_holdout}()}
#' @export
#'
#' @examples
#' \donttest{
#' # Fitting the cross-fitted lasso paths takes more than five seconds.
#' set.seed(2026)
#' n <- 450
#' p <- 20
#' tau <- c(150, 300)
#' b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5))
#' delta <- cbind(-2 * b0, 1.8 * b0)
#' dat <- dgp_linear_regression(n, p, tau, b0, delta, sig = 1)$data
#' fit <- reliever(
#'   X = dat[, -1, drop = FALSE], y = dat[, 1],
#'   cpd_family = "lasso_crossfit",
#'   cpn_max = 5, dm = 15, cov_rate = 0.8, method = "SN",
#'   nfolds = 2,
#'   loss_output_types = c("recv", "crossfit_homo_hyper")
#' )
#' # crossfit_homo_hyper means cross-fitted loss with one interval-independent
#' # hyperparameter value.
#' selected <- select_across_runs(
#'   result = fit, run_type = "crossfit_homo_hyper", cpn_crit = "loss"
#' )
#' selected
#' stopifnot(identical(selected$K_hat, 2L))
#' }
select_across_runs <- function(result, run_ids = NULL, run_type = NULL,
                               cpn_crit,
                               n = NULL) {
  if (is.null(run_ids) && is.null(run_type)) {
    stop("run_ids or run_type must be supplied for across-run selection.",
         call. = FALSE)
  }
  if (missing(cpn_crit)) {
    stop("cpn_crit must be supplied.", call. = FALSE)
  }
  resolved_run_ids <- .model_select_run_ids(
    result, run_ids = run_ids, run_type = run_type
  )
  if (!is.null(run_type)) {
    run_type_label <- paste0(
      "\"", unique(as.character(run_type)), "\"",
      collapse = ", "
    )
    recommended_types <- c("recv", "crossfit_homo_hyper")
    if (any(!unique(as.character(run_type)) %in% recommended_types)) {
      warning(warningCondition(
        paste0(
          "run_type = ", run_type_label,
          " may select incomparable losses; selection will continue. Prefer ",
          "crossfit or cv.reliever() for built-in models; custom reg_fun users ",
          "should prefer run_ids."
        ),
        call = NULL,
        class = "reliever_across_run_comparability_warning"
      ))
    }
    if (length(resolved_run_ids) == 1L) {
      warning(warningCondition(
        paste0(
          "run_type = ", run_type_label,
          " selects one run; selection will continue. Use select_by_run() ",
          "when selecting only K."
        ),
        call = NULL,
        class = "reliever_single_run_selection_warning"
      ))
    }
  }
  selected_meta <- .model_select_metadata(result)
  selected_meta <- selected_meta[
    match(resolved_run_ids, selected_meta$run_id), , drop = FALSE
  ]
  .reliever_warn_incomparable_loss_kinds(selected_meta)
  selection <- .model_select_candidates(
    result, run_ids = resolved_run_ids,
    cpn_crit = cpn_crit, n = n
  )
  selected <- .model_select_one(selection$candidates)
  if (is.null(selected)) {
    stop("No finite candidate is available for selection.", call. = FALSE)
  }
  .model_select_add_metadata(
    .model_select_format(selected, selection$penalty_info), result
  )
}
