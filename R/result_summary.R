.reliever_summary_columns <- c("rule", "K_hat", "cpd_hat")

.reliever_empty_summary <- function() {
  out <- data.frame(
    rule = character(),
    K_hat = integer(),
    stringsAsFactors = FALSE
  )
  out$cpd_hat <- I(vector("list", 0L))
  out <- out[.reliever_summary_columns]
  class(out) <- c("reliever_summary", "data.frame")
  out
}

.reliever_selection_summary <- function(selected) {
  required <- c("rule", "K_hat", "cpd_hat")
  if (!is.data.frame(selected) || !all(required %in% names(selected))) {
    stop(
      "selected must contain rule, K_hat, and cpd_hat.",
      call. = FALSE
    )
  }
  columns <- c(
    "rule",
    if ("row_type" %in% names(selected)) "row_type",
    if ("hyper_value" %in% names(selected)) "hyper_value",
    "K_hat",
    "cpd_hat"
  )
  out <- as.data.frame(selected[columns])
  rownames(out) <- NULL
  class(out) <- c("reliever_summary", "data.frame")
  hyper_name <- attr(selected, "hyper_name", exact = TRUE)
  if (!is.null(hyper_name) && "hyper_value" %in% names(out)) {
    attr(out, "hyper_name") <- hyper_name
  }
  out
}

.reliever_drop_null <- function(x) {
  x[!vapply(x, is.null, logical(1L))]
}

.reliever_result_settings <- function(setup, cpn_crit, cpn_penalty_info,
                                      pen_val, wbs_seed, owner_key,
                                      ratio = NULL, dc_grid_size = NULL) {
  if (is.character(cpn_crit)) {
    cpn_crit <- cpn_penalty_info$label
  }
  search_n <- if (setup$n != setup$n_original) setup$n
  search_dm <- if (setup$dm != setup$dm_original) setup$dm
  is_wbs <- setup$method %in% c("WBS", "WBS_recursive")
  is_wbs_family <- setup$method %in%
    c("WBS", "WBS_recursive", "SeedBS", "BS")

  .reliever_drop_null(list(
    method = setup$method,
    n = setup$n_original,
    dm = setup$dm_original,
    cpn_max = setup$cpn_max,
    search_n = search_n,
    search_dm = search_dm,
    dc_grid_size = dc_grid_size,
    dc_grid = setup$dc_grid,
    cov_rate = setup$cov_rate,
    pen_val = if (setup$method %in% c("PELT", "OP")) as.numeric(pen_val),
    prune_value = if (setup$method == "PELT") setup$prune_value,
    M = if (is_wbs) setup$M,
    wbs_seed = if (is_wbs) wbs_seed,
    wbs_stop_crit = if (is_wbs_family) setup$wbs_stop_crit,
    cache_backend = setup$cache_backend,
    owner_key = if (setup$cache_backend == "by_loss_block") isTRUE(owner_key),
    ratio = ratio,
    cpn_crit = cpn_crit,
    cpn_penalty = cpn_penalty_info$value,
    cpn_criterion = cpn_penalty_info$criterion
  ))
}

.reliever_column_has_value <- function(x) {
  values <- if (is.list(x)) x else as.list(x)
  any(vapply(values, function(value) {
    if (length(value) == 0L) {
      return(FALSE)
    }
    if (is.character(value)) {
      return(any(!is.na(value) & nzchar(value)))
    }
    !is.atomic(value) || any(!is.na(value))
  }, logical(1L)))
}

.reliever_simplify_meta_column <- function(x) {
  if (!is.list(x) || !all(vapply(x, length, integer(1L)) == 1L) ||
      !all(vapply(x, is.atomic, logical(1L)))) {
    return(x)
  }
  unlist(x, recursive = FALSE, use.names = FALSE)
}

.reliever_hyper_name <- function(meta, run_ids = NULL) {
  if (!is.data.frame(meta) || !"hyper_name" %in% names(meta)) {
    return(NULL)
  }
  if (!is.null(run_ids) && "run_id" %in% names(meta)) {
    meta <- meta[meta$run_id %in% run_ids, , drop = FALSE]
  }
  values <- unique(as.character(meta$hyper_name))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 1L) values else NULL
}

.reliever_default_run_ids <- function(meta, unmarked = c("all", "none")) {
  unmarked <- match.arg(unmarked)
  if (is.null(meta) || !"default_selection" %in% names(meta)) {
    return(NULL)
  }
  default_selection <- meta$default_selection
  if (!is.logical(default_selection) || anyNA(default_selection)) {
    stop("run_meta$default_selection must contain TRUE or FALSE.",
         call. = FALSE)
  }
  run_ids <- as.integer(meta$run_id[default_selection])
  if (length(run_ids) > 0L || unmarked == "none") {
    return(run_ids)
  }
  NULL
}

.reliever_hyper_label <- function(hyper_name) {
  if (is.null(hyper_name)) {
    return("Hyperparameter")
  }
  label <- gsub("_", " ", hyper_name, fixed = TRUE)
  paste0(toupper(substr(label, 1L, 1L)), substring(label, 2L))
}

.reliever_validate_max_rows <- function(max_rows) {
  if (!is.numeric(max_rows) || length(max_rows) != 1L || is.na(max_rows) ||
      max_rows <= 0 || (!is.infinite(max_rows) && max_rows != floor(max_rows))) {
    stop("max_rows must be a positive integer or Inf.", call. = FALSE)
  }
  max_rows
}

.reliever_print_table <- function(x, max_rows = 10L, ...) {
  shown <- x
  hyper_name <- attr(x, "hyper_name", exact = TRUE)
  if (!is.null(hyper_name) && "hyper_value" %in% names(shown) &&
      !hyper_name %in% setdiff(names(shown), "hyper_value")) {
    names(shown)[names(shown) == "hyper_value"] <- hyper_name
  }
  n_total <- nrow(shown)
  if (is.finite(max_rows) && n_total > max_rows) {
    shown <- shown[seq_len(as.integer(max_rows)), , drop = FALSE]
  }
  if ("cpd_hat" %in% names(shown)) {
    shown$cpd_hat <- vapply(shown$cpd_hat, function(cpd) {
      if (length(cpd) == 0L) "<none>" else paste(cpd, collapse = ", ")
    }, character(1L))
  }
  print.data.frame(shown, row.names = FALSE, ...)
  if (nrow(shown) < n_total) {
    cat(sprintf(
      "Showing first %d of %d rows; increase max_rows to display more.\n",
      nrow(shown), n_total
    ))
  }
}

#' Print a compact Reliever summary
#'
#' Print every changepoint location without changing the data-frame interface
#' returned by \code{summary()} and stored in the corresponding result
#' object's \code{summary} component. When the rows select K separately for
#' several hyperparameter values, the print method also explains that no
#' hyperparameter has been selected across those rows. An empty summary
#' explains that fitting and model selection are separate steps and points to
#' the recommended workflows in \code{\link{reliever}()}.
#'
#' @param x A compact summary from a Reliever or outer-CV result.
#' @param max_rows Maximum number of rows printed. Use \code{Inf} to print all
#'   rows.
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#'
#' @return \code{x}, invisibly.
#' @export
print.reliever_summary <- function(x, max_rows = 10L, ...) {
  max_rows <- .reliever_validate_max_rows(max_rows)
  if (nrow(x) == 0L) {
    .reliever_print_unselected_workflow()
    return(invisible(x))
  }
  .reliever_print_table(x, max_rows = max_rows, ...)
  selection_note <- attr(x, "selection_note", exact = TRUE)
  if (!is.null(selection_note)) {
    cat(selection_note, "\n", sep = "")
  }
  invisible(x)
}

.reliever_print_unselected_workflow <- function() {
  cat(
    "No model or changepoint number has been selected.\n",
    "By default, reliever() fits candidate changepoint models; it does not ",
    "select a model or K.\n",
    "For a compatible single loss path or ordinary lasso path, use ",
    "cv.reliever() to select K and, when applicable, the model setting; use ",
    "cv.reliever_generic() for a custom reg_fun.\n",
    "For a cross-fitted loss family, select the ReCV path with ",
    "select_by_run(..., run_type = \"recv\", cpn_crit = \"loss\").\n",
    "See ?reliever for the recommended workflows.\n",
    sep = ""
  )
}

.reliever_build_summary <- function(result) {
  meta <- result$run_meta
  cpn_crit <- result$settings$cpn_crit
  default_run_ids <- .reliever_default_run_ids(meta, unmarked = "none")
  if (!is.null(default_run_ids) && length(default_run_ids) == 0L) {
    return(.reliever_empty_summary())
  }
  if (is.null(cpn_crit)) {
    cpn_crit <- "none"
  }
  if (.reliever_cpn_is_none(cpn_crit)) {
    selector_map <- result$cpd_path$selector_map
    if (is.null(selector_map)) {
      return(.reliever_empty_summary())
    }
    if (!is.null(default_run_ids)) {
      selector_map <- selector_map[
        selector_map$run_id %in% default_run_ids, , drop = FALSE
      ]
    }
    selector_map <- selector_map[is.finite(selector_map$select_value), ,
                                 drop = FALSE]
    if (nrow(selector_map) == 0L) {
      return(.reliever_empty_summary())
    }
    candidates <- result$cpd_path$candidates
    candidate_id <- match(
      paste(selector_map$run_id, selector_map$candidate_id),
      paste(candidates$run_id, candidates$candidate_id)
    )
    if (anyNA(candidate_id)) {
      stop("A path selector refers to an unavailable candidate.",
           call. = FALSE)
    }
    out <- data.frame(selector_map$select_value, stringsAsFactors = FALSE)
    names(out) <- result$cpd_path$select_by
    out$K_hat <- as.integer(candidates$K[candidate_id])
    out$cpd_hat <- I(lapply(candidates$cpd[candidate_id], as.integer))
    meta_id <- if (!is.null(meta)) {
      match(selector_map$run_id, meta$run_id)
    } else {
      rep.int(NA_integer_, nrow(selector_map))
    }
    if (!is.null(meta) && "row_type" %in% names(meta)) {
      row_types <- as.character(meta$row_type)
      valid_types <- unique(row_types[!is.na(row_types) & nzchar(row_types)])
      if (length(valid_types) > 1L) {
        out$row_type <- row_types[meta_id]
      }
    }
    if (!is.null(meta) && "hyper_value" %in% names(meta)) {
      hyper_value <- meta$hyper_value[meta_id]
      if (.reliever_column_has_value(hyper_value)) {
        out$hyper_value <- .reliever_simplify_meta_column(hyper_value)
      }
    }
    out <- out[c(
      intersect("row_type", names(out)),
      intersect("hyper_value", names(out)),
      result$cpd_path$select_by, "K_hat", "cpd_hat"
    )]
    rownames(out) <- NULL
    return(out)
  }

  selected <- select_by_run(
    result, run_ids = default_run_ids, cpn_crit = cpn_crit
  )
  if ("hyper_value" %in% names(selected) &&
      !.reliever_column_has_value(selected$hyper_value)) {
    selected$hyper_value <- NULL
  }
  .reliever_selection_summary(selected)
}

.reliever_finalize_result <- function(result) {
  loss_spec <- result$loss_spec
  diagnostics <- result$diagnostics
  if (is.null(diagnostics)) {
    diagnostics <- list()
  }
  diagnostics <- .reliever_drop_null(diagnostics)
  diagnostics <- diagnostics[!vapply(diagnostics, function(x) {
    length(x) == 0L
  }, logical(1L))]

  settings <- .reliever_drop_null(result$settings)
  if (!is.null(loss_spec) && length(loss_spec$args) > 0L) {
    settings$loss_args <- loss_spec$args
  }

  timing_values <- result[c(
    "n_model_fit", "model_fit_time", "cpd_time", "total_time"
  )]
  timing <- if (length(timing_values) == 4L &&
                !any(vapply(timing_values, is.null, logical(1L)))) {
    data.frame(
      run_id = result$run_meta$run_id,
      n_model_fit = as.numeric(result$n_model_fit),
      model_fit_time = as.numeric(result$model_fit_time),
      cpd_time = as.numeric(result$cpd_time),
      total_time = as.numeric(result$total_time),
      stringsAsFactors = FALSE
    )
  }

  summary <- .reliever_build_summary(result)
  class(summary) <- c("reliever_summary", "data.frame")
  hyper_name <- .reliever_hyper_name(result$run_meta)
  if (!is.null(hyper_name) && "hyper_value" %in% names(summary)) {
    attr(summary, "hyper_name") <- hyper_name
  }
  if ("hyper_value" %in% names(summary)) {
    hyper_values <- if (is.list(summary$hyper_value)) {
      summary$hyper_value
    } else {
      as.list(summary$hyper_value)
    }
    hyper_values <- hyper_values[vapply(hyper_values, function(value) {
      length(value) > 0L && !(length(value) == 1L && is.na(value))
    }, logical(1L))]
    if (length(unique(hyper_values)) > 1L) {
      hyper_term <- tolower(.reliever_hyper_label(hyper_name))
      attr(summary, "selection_note") <- paste0(
        "K was selected separately within each ", hyper_term,
        " value. This summary does not choose one ", hyper_term,
        "; use fixed-hyperparameter cross-fitted losses from a crossfit fit, ",
        "outer CV, or independent hold-out loss for a joint choice."
      )
    }
  }

  out <- .reliever_drop_null(list(
    summary = summary,
    cpd_path = result$cpd_path,
    run_meta = result$run_meta,
    settings = settings,
    timing = timing,
    diagnostics = if (length(diagnostics) == 0L) NULL else diagnostics,
    cache_profile = result$cache_profile
  ))
  class(out) <- c("reliever_result", "list")
  if (!is.null(loss_spec)) {
    attr(out, ".reliever_loss_spec") <- loss_spec
  }
  out
}

.reliever_hyper_value_label <- function(value, index) {
  if (is.data.frame(value) && nrow(value) == 1L) {
    value <- as.list(value[1L, , drop = FALSE])
  } else if (is.matrix(value) && nrow(value) == 1L) {
    value <- as.list(value[1L, , drop = TRUE])
  }
  if (is.list(value) && length(value) == 1L &&
      is.null(names(value))) {
    value <- value[[1L]]
  }
  if (is.atomic(value) && length(value) > 0L && !anyNA(value)) {
    return(paste(format(value, digits = 5L, trim = TRUE), collapse = ", "))
  }
  if (is.list(value) && length(value) > 0L) {
    simple <- vapply(value, function(element) {
      is.atomic(element) && length(element) == 1L && !is.na(element)
    }, logical(1L))
    if (all(simple)) {
      labels <- vapply(value, function(element) {
        format(element, digits = 5L, trim = TRUE)
      }, character(1L))
      value_names <- names(value)
      if (!is.null(value_names) && all(nzchar(value_names))) {
        labels <- paste0(value_names, "=", labels)
      }
      return(paste(labels, collapse = ", "))
    }
  }
  paste("setting", index)
}

.reliever_plot_hyper_axis <- function(values) {
  values <- if (is.list(values)) values else as.list(values)
  numeric_value <- vapply(values, function(value) {
    if (!is.numeric(value) || length(value) != 1L ||
        is.na(value) || !is.finite(value)) {
      return(NA_real_)
    }
    as.numeric(value)
  }, numeric(1L))
  if (!anyNA(numeric_value)) {
    return(list(
      position = numeric_value,
      display = numeric_value,
      tick_at = NULL,
      tick_label = NULL,
      is_numeric = TRUE
    ))
  }
  labels <- vapply(seq_along(values), function(index) {
    .reliever_hyper_value_label(values[[index]], index)
  }, character(1L))
  list(
    position = seq_along(values),
    display = labels,
    tick_at = seq_along(values),
    tick_label = labels,
    is_numeric = FALSE
  )
}

.reliever_search_axis_label <- function(select_by) {
  switch(
    select_by,
    wbs_stop_crit = "WBS stopping criterion",
    pen_val = "Penalty",
    select_by
  )
}

.reliever_plot_selected_keys <- function(candidates) {
  selected <- lapply(split(candidates, candidates$run_id), .model_select_one)
  vapply(selected, function(row) {
    if (is.null(row)) return(NA_character_)
    paste(row$run_id, row$candidate_id)
  }, character(1L))
}

.reliever_plot_run_labels <- function(result, run_ids) {
  labels <- paste("run", run_ids)
  meta <- result$run_meta
  if (is.null(meta)) {
    return(stats::setNames(labels, run_ids))
  }
  meta <- meta[match(run_ids, meta$run_id), , drop = FALSE]
  if ("hyper_value" %in% names(meta)) {
    values <- if (is.list(meta$hyper_value)) {
      meta$hyper_value
    } else {
      as.list(meta$hyper_value)
    }
    has_value <- vapply(values, function(value) {
      length(value) == 1L && !is.na(value)
    }, logical(1L))
    labels[has_value] <- vapply(values[has_value], function(value) {
      paste0("hyper=", format(value, digits = 5L, trim = TRUE))
    }, character(1L))
  }
  if ("row_type" %in% names(meta) && length(unique(meta$row_type)) > 1L) {
    labels <- paste0(meta$row_type, ": ", labels)
  }
  stats::setNames(labels, run_ids)
}

.reliever_draw_path <- function(data, x_value, y_value, group, selected,
                                 se = NULL, labels = NULL, col = NULL,
                                 lty = 1, pch = 19, xlab, ylab, main, ...) {
  tick_at <- attr(x_value, "tick_at", exact = TRUE)
  tick_label <- attr(x_value, "tick_label", exact = TRUE)
  keep <- is.finite(x_value) & is.finite(y_value)
  if (!any(keep)) {
    stop("No finite values are available to plot.", call. = FALSE)
  }
  data <- data[keep, , drop = FALSE]
  x_value <- x_value[keep]
  y_value <- y_value[keep]
  group <- group[keep]
  selected <- selected[keep]
  if (!is.null(se)) {
    se <- se[keep]
  }

  groups <- unique(group)
  n_groups <- length(groups)
  if (is.null(col)) {
    col <- grDevices::hcl.colors(max(3L, n_groups), "Dark 3")[
      seq_len(n_groups)
    ]
  } else {
    col <- rep(col, length.out = n_groups)
  }
  lty <- rep(lty, length.out = n_groups)
  pch <- rep(pch, length.out = n_groups)

  dots <- list(...)
  dots$type <- NULL
  if (!is.null(tick_at)) {
    dots$xaxt <- "n"
  }
  if (is.null(dots$ylim) && !is.null(se)) {
    finite_se <- is.finite(se)
    dots$ylim <- range(c(
      y_value[finite_se] - se[finite_se],
      y_value[finite_se] + se[finite_se],
      y_value
    ), finite = TRUE)
  }
  do.call(graphics::plot.default, c(list(
    x = x_value,
    y = y_value,
    type = "n",
    xlab = xlab,
    ylab = ylab,
    main = main
  ), dots))
  if (!is.null(tick_at)) {
    graphics::axis(1, at = tick_at, labels = tick_label)
  }

  for (i in seq_along(groups)) {
    id <- which(group == groups[i])
    id <- id[order(x_value[id], y_value[id])]
    if (length(id) > 1L && !anyDuplicated(x_value[id])) {
      graphics::lines(x_value[id], y_value[id], col = col[i], lty = lty[i])
    }
    if (!is.null(se)) {
      finite_se <- id[is.finite(se[id])]
      if (length(finite_se)) {
        graphics::arrows(
          x_value[finite_se], y_value[finite_se] - se[finite_se],
          x_value[finite_se], y_value[finite_se] + se[finite_se],
          angle = 90, code = 3, length = 0.03, col = col[i]
        )
      }
    }
    graphics::points(x_value[id], y_value[id], col = col[i], pch = pch[i])
  }

  if (any(selected)) {
    graphics::points(
      x_value[selected], y_value[selected], pch = 1, cex = 1.8, lwd = 2
    )
  }
  if (n_groups > 1L && n_groups <= 12L) {
    if (is.null(labels)) {
      labels <- paste("run", groups)
    } else if (!is.null(names(labels))) {
      labels <- labels[as.character(groups)]
    }
    graphics::legend(
      "topright", legend = labels, col = col, lty = lty, pch = pch,
      bty = "n"
    )
  }
  invisible(data)
}

#' Plot Reliever model-selection profiles
#'
#' Plot candidate scores against the number of changepoints, inspect a
#' solution path across numeric or categorical hyperparameter values, or plot
#' a WBS stopping criterion or PELT/OP penalty on its original scale. The
#' default uses the selection criterion stored in the fit. Set
#' \code{cpn_crit = "loss"} to show and minimize the unpenalized fitted loss
#' instead.
#'
#' With \code{x_axis = "K"}, one curve is drawn for every requested run. A
#' selected K is circled only when the fit stores a selection criterion or
#' \code{cpn_crit} is supplied. With
#' \code{x_axis = "hyperparameter"}, \code{K} fixes the changepoint number and
#' produces a fixed-K loss or criterion curve. If \code{K = NULL}, a stored or
#' supplied criterion is required and first selects K separately within each
#' run. Hyperparameter plots require one non-missing \code{hyper_value} per run
#' in \code{x$run_meta}. Numeric values use their numeric scale; categorical
#' values follow their declared run order.
#' If the corresponding metadata also contains one common
#' \code{hyper_name}, that name is used for the horizontal-axis label and plot
#' title.
#'
#' Raw fitted losses from an ordinary solution path are not generally
#' comparable across hyperparameters. Such a profile is therefore shown
#' without circling a preferred hyperparameter or interpreting an endpoint
#' minimum. For a ReCV result, omitting both \code{run_ids} and \code{run_type}
#' uses the default \code{recv} run. Select fixed-hyperparameter cross-fitted
#' runs explicitly with \code{run_type = "crossfit_homo_hyper"}. When a
#' selection criterion is stored or supplied, those comparable losses have
#' their minimum circled, and a fixed-K minimum at either grid endpoint
#' triggers a warning because the candidate range may be too narrow.
#' Outer-CV profiles have the same grid-check behavior; see
#' \code{plot.cv_reliever_result()}.
#'
#' With \code{x_axis = "search_value"}, the horizontal axis is the supplied
#' \code{wbs_stop_crit} or \code{pen_val}. This is useful when several search
#' values produce the same K. If the fit has \code{cpn_crit = "none"}, the
#' loss profile is shown without circling a supposedly selected K. Supply an
#' explicit \code{cpn_crit}, including \code{"loss"}, when a selected point is
#' wanted.
#'
#' Two-step searches do not compute a comparable whole-segmentation loss, so
#' their candidate table cannot be drawn as a loss or IC profile. Inspect
#' \code{x$cpd_path$candidates} instead.
#'
#' @param x A result returned by \code{reliever()} or a related wrapper.
#' @param x_axis Horizontal-axis choice:
#'   \itemize{
#'     \item \code{"K"} for changepoint number;
#'     \item \code{"hyperparameter"} for values from \code{x$run_meta}; or
#'     \item \code{"search_value"} for WBS, PELT, or OP controls.
#'   }
#' @param K Optional non-negative changepoint number. It is used only for a
#'   fixed-K hyperparameter plot.
#' @param run_ids Optional positive integer run identifiers. By default, loss
#'   outputs marked as preferred are used when available; for cross-fitted
#'   families this is the \code{recv} run. Otherwise all runs are shown.
#'   Supply only one of \code{run_ids} and \code{run_type}.
#' @param run_type Optional stored row types, such as \code{"recv"} or
#'   \code{"crossfit_homo_hyper"}. Supply only one of \code{run_type} and
#'   \code{run_ids}.
#' @param cpn_crit Vertical-axis criterion. The default \code{NULL} uses the
#'   criterion stored in the fit. Use \code{"loss"} for raw fitted loss; the
#'   post-fit selector documentation lists the other criteria.
#' @param col,lty,pch Graphical parameters for curves and points.
#' @param xlab,ylab,main Optional axis labels and title.
#' @param ... Additional arguments passed to \code{graphics::plot.default()}.
#'
#' @return The data used for plotting, invisibly. It includes a logical
#'   \code{selected} column identifying the circled points. This column is all
#'   \code{FALSE} when the fit has no K-selection rule and none is supplied, or
#'   when an ordinary fitted-loss solution path is plotted across
#'   hyperparameters that are not directly comparable.
#' @seealso \code{\link{summary.reliever_result}()},
#'   \code{\link{select_by_run}()},
#'   \code{\link{plot.cv_reliever_result}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' x <- rbind(
#'   matrix(rnorm(300 * 3), 300, 3),
#'   matrix(rnorm(300 * 3, mean = 3), 300, 3),
#'   matrix(rnorm(300 * 3, mean = -3), 300, 3)
#' )
#' fit <- reliever(X = x, cpn_max = 5, dm = 30, cov_rate = 0.8)
#' plot(x = fit)
plot.reliever_result <- function(x, x_axis = c("K", "hyperparameter",
                                               "search_value"),
                                 K = NULL, run_ids = NULL, cpn_crit = NULL,
                                 col = NULL, lty = 1, pch = 19,
                                 xlab = NULL, ylab = NULL, main = NULL, ...,
                                 run_type = NULL) {
  x_axis <- match.arg(x_axis)
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

  if (is.null(run_ids) && is.null(run_type)) {
    run_ids <- .reliever_default_run_ids(x$run_meta)
  }
  if (is.null(cpn_crit)) cpn_crit <- x$settings$cpn_crit
  has_selection <- !is.null(cpn_crit) &&
    !.reliever_cpn_is_none(cpn_crit)
  score_crit <- if (has_selection) cpn_crit else "loss"
  if (all(is.na(x$cpd_path$candidates$loss))) {
    stop(
      "This result does not contain whole-segmentation losses, so a loss or criterion path cannot be plotted. For twostep() results, inspect x$cpd_path$candidates.",
      call. = FALSE
    )
  }
  selection <- .model_select_candidates(
    x, run_ids = run_ids, run_type = run_type, cpn_crit = score_crit
  )
  plot_data <- selection$candidates
  run_ids <- selection$run_ids

  if (x_axis == "search_value") {
    selector_map <- x$cpd_path$selector_map
    select_by <- x$cpd_path$select_by
    if (is.null(selector_map) || is.null(select_by)) {
      stop(
        "This result is indexed by K and has no WBS/PELT/OP search values.",
        call. = FALSE
      )
    }
    selector_map <- selector_map[
      selector_map$run_id %in% run_ids &
        is.finite(selector_map$select_value), , drop = FALSE
    ]
    candidate_id <- match(
      paste(selector_map$run_id, selector_map$candidate_id),
      paste(plot_data$run_id, plot_data$candidate_id)
    )
    if (nrow(selector_map) == 0L || anyNA(candidate_id)) {
      stop("No finite search values are available to plot.", call. = FALSE)
    }
    plot_data <- plot_data[candidate_id, , drop = FALSE]
    plot_data[[select_by]] <- selector_map$select_value
    if (has_selection) {
      selected_key <- .reliever_plot_selected_keys(selection$candidates)
      plot_data$selected <-
        paste(plot_data$run_id, plot_data$candidate_id) %in% selected_key
    } else {
      plot_data$selected <- FALSE
    }
    x_value <- plot_data[[select_by]]
    group <- plot_data$run_id
    axis_label <- .reliever_search_axis_label(select_by)
    xlab <- if (is.null(xlab)) axis_label else xlab
    main <- if (is.null(main)) "Search-value profile" else main
  } else if (x_axis == "hyperparameter") {
    if (is.null(K)) {
      if (!has_selection) {
        stop(
          "K must be supplied for a hyperparameter plot when the fit has ",
          "no K-selection rule; alternatively supply cpn_crit.",
          call. = FALSE
        )
      }
      selected_rows <- lapply(split(plot_data, plot_data$run_id),
                              .model_select_one)
      selected_rows <- selected_rows[!vapply(selected_rows, is.null,
                                              logical(1L))]
      plot_data <- do.call(rbind, selected_rows)
    } else {
      plot_data <- plot_data[plot_data$K == K, , drop = FALSE]
      if (nrow(plot_data) == 0L) {
        stop("No requested run contains the selected K.", call. = FALSE)
      }
      if (anyDuplicated(plot_data$run_id)) {
        stop(
          "K does not identify one candidate per run; use x_axis = \"K\".",
          call. = FALSE
        )
      }
    }
    meta <- .model_select_metadata(x)
    meta_id <- match(plot_data$run_id, meta$run_id)
    if (!"hyper_value" %in% names(meta) || anyNA(meta_id)) {
      stop("The selected runs do not contain hyper_value metadata.",
           call. = FALSE)
    }
    hyper_value <- meta$hyper_value[meta_id]
    has_hyper_value <- vapply(
      as.list(hyper_value), .reliever_column_has_value, logical(1L)
    )
    if (!all(has_hyper_value)) {
      stop(
        "The selected runs do not contain non-missing hyper_value metadata. ",
        "For a cross-fitted hyperparameter plot, supply ",
        "run_type = \"crossfit_homo_hyper\".",
        call. = FALSE
      )
    }
    comparable_hyper_loss <- "row_type" %in% names(meta) &&
      isTRUE(all(meta$row_type[meta_id] == "crossfit_homo_hyper"))
    hyper_name <- .reliever_hyper_name(meta, plot_data$run_id)
    hyper_label <- .reliever_hyper_label(hyper_name)
    hyper_axis <- .reliever_plot_hyper_axis(hyper_value)
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
    best <- .model_select_one(plot_data)
    if (is.null(best)) {
      stop("No finite values are available to plot.", call. = FALSE)
    }
    plot_data$selected <- if (has_selection && comparable_hyper_loss) {
      paste(plot_data$run_id, plot_data$candidate_id) ==
        paste(best$run_id, best$candidate_id)
    } else {
      FALSE
    }
    selected_id <- which(plot_data$selected)
    if (hyper_axis$is_numeric && comparable_hyper_loss &&
        has_selection && !is.null(K) &&
        nrow(plot_data) > 1L && length(selected_id) == 1L &&
        selected_id %in% c(1L, nrow(plot_data))) {
      boundary <- if (selected_id == 1L) "lower" else "upper"
      warning(
        "The fixed-K ", tolower(hyper_label),
        " profile is minimized at the ", boundary,
        " endpoint of the ", tolower(hyper_label),
        " grid; consider a wider candidate range.",
        call. = FALSE
      )
    }
    x_value <- hyper_position
    if (!hyper_axis$is_numeric) {
      attr(x_value, "tick_at") <- hyper_position
      attr(x_value, "tick_label") <- plot_data$hyper_value
    }
    group <- rep.int(1L, nrow(plot_data))
    xlab <- if (is.null(xlab)) hyper_label else xlab
    main <- if (is.null(main)) {
      if (comparable_hyper_loss) {
        paste("Cross-fitted", tolower(hyper_label), "profile")
      } else {
        paste("Fitted-loss", tolower(hyper_label), "profile")
      }
    } else {
      main
    }
  } else {
    if (has_selection) {
      selected_key <- .reliever_plot_selected_keys(plot_data)
      plot_data$selected <-
        paste(plot_data$run_id, plot_data$candidate_id) %in% selected_key
    } else {
      plot_data$selected <- FALSE
    }
    x_value <- plot_data$K
    group <- plot_data$run_id
    xlab <- if (is.null(xlab)) "Number of changepoints (K)" else xlab
    main <- if (is.null(main)) "Changepoint selection path" else main
  }

  ylab <- if (is.null(ylab)) {
    if (selection$penalty_info$criterion == "loss" &&
        selection$penalty_info$value == 0) {
      "Loss"
    } else {
      paste0("Selection score (", selection$penalty_info$label, ")")
    }
  } else {
    ylab
  }
  labels <- .reliever_plot_run_labels(x, unique(group))
  .reliever_draw_path(
    data = plot_data,
    x_value = x_value,
    y_value = plot_data$score,
    group = group,
    selected = plot_data$selected,
    labels = labels,
    col = col,
    lty = lty,
    pch = pch,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
}

#' Plot data with selected changepoints
#'
#' Draw one or more observed series against their row index and add vertical
#' lines at selected changepoints. Compact Reliever results do not retain the
#' fitting input, so the data are supplied explicitly. A non-native fit made
#' with \code{detail = TRUE} or \code{"cache"} may retain it inside its
#' reusable cache profile.
#'
#' By default, the function uses the only row of \code{summary(result)}. If the
#' summary has several rows, such as one per solution-path hyperparameter,
#' supply either a row number or a one-row selection table. A row returned by
#' \code{select_by_run()},
#' \code{select_across_runs()}, or \code{select_holdout()} can therefore be
#' plotted directly.
#'
#' Changepoint locations identify the final observation in the preceding
#' segment, so their vertical lines are drawn halfway between that observation
#' and the next one.
#'
#' @param result A \code{reliever_result} or \code{cv_reliever_result}.
#' @param data The numeric vector, matrix, or data frame used for fitting. Its
#'   number of rows must match \code{result$settings$n}.
#' @param selection Which changepoint selection to plot. The default
#'   \code{NULL} uses \code{summary(result)} when it has exactly one row. A
#'   positive integer chooses one summary row. Alternatively, supply a
#'   one-row data frame containing a \code{cpd_hat} list-column.
#' @param columns Numeric or character columns of \code{data} to draw. The
#'   default draws the first column.
#' @param type,col,lty Graphical parameters for the observed series.
#' @param cpd_col,cpd_lty,cpd_lwd Graphical parameters for the changepoint
#'   lines.
#' @param legend Add a legend when more than one series is drawn.
#' @param xlab,ylab,main Optional axis labels and title.
#' @param ... Additional arguments passed to \code{graphics::matplot()}.
#'
#' @return The plotted data, invisibly, with the selected changepoints in the
#'   \code{"cpd_hat"} attribute.
#' @seealso \code{\link{plot.reliever_result}()},
#'   \code{\link{select_by_run}()}
#' @export
#'
#' @examples
#' set.seed(2026)
#' n_seg <- 300
#' x <- c(
#'   rnorm(n_seg, 0, 0.5),
#'   rnorm(n_seg, 4, 0.5),
#'   rnorm(n_seg, -4, 0.5)
#' )
#' fit <- reliever(X = x, cpn_max = 7, dm = 30, cov_rate = 0.8)
#' selected <- select_by_run(result = fit, cpn_crit = "rss_sic")
#' stopifnot(identical(selected$K_hat, 2L))
#' plot_reliever_data(result = fit, data = x, selection = selected)
plot_reliever_data <- function(result, data, selection = NULL, columns = 1L,
                               type = "l", col = NULL, lty = 1,
                               cpd_col = 2, cpd_lty = 2, cpd_lwd = 1,
                               legend = TRUE, xlab = "Observation",
                               ylab = NULL,
                               main = "Data and selected changepoints", ...) {
  if (!inherits(result, c("reliever_result", "cv_reliever_result"))) {
    stop("result must be a reliever_result or cv_reliever_result.",
         call. = FALSE)
  }
  if (is.data.frame(data)) {
    is_numeric <- vapply(data, is.numeric, logical(1L))
    if (!all(is_numeric)) {
      stop("Every data-frame column must be numeric.",
           call. = FALSE)
    }
    data <- as.matrix(data)
  } else if (is.numeric(data) && is.null(dim(data))) {
    data <- matrix(data, ncol = 1L, dimnames = list(NULL, "value"))
  } else if (!is.matrix(data) || !is.numeric(data)) {
    stop("data must be a numeric vector, matrix, or data frame.",
         call. = FALSE)
  }
  if (ncol(data) == 0L || nrow(data) == 0L || any(!is.finite(data))) {
    stop("data must contain finite numeric observations.", call. = FALSE)
  }
  expected_n <- result$settings$n
  if (!is.null(expected_n) && nrow(data) != expected_n) {
    stop(
      "nrow(data) must match result$settings$n (", expected_n, ").",
      call. = FALSE
    )
  }

  result_summary <- summary(result)
  if (is.null(selection)) {
    if (nrow(result_summary) == 0L) {
      stop(
        "result has no selected summary row; supply a selection returned by ",
        "select_by_run(), select_across_runs(), or select_holdout().",
        call. = FALSE
      )
    }
    if (nrow(result_summary) != 1L) {
      stop(
        "result has ", nrow(result_summary),
        " summary rows; supply selection as a row number or a one-row ",
        "selection table.",
        call. = FALSE
      )
    }
    selection <- result_summary
  } else if (is.numeric(selection) && length(selection) == 1L &&
             !is.na(selection) && selection == floor(selection)) {
    if (selection < 1L || selection > nrow(result_summary)) {
      stop("selection row number is outside summary(result).",
           call. = FALSE)
    }
    selection <- result_summary[as.integer(selection), , drop = FALSE]
  } else if (!is.data.frame(selection) || nrow(selection) != 1L ||
             !"cpd_hat" %in% names(selection)) {
    stop(
      "selection must be a summary row number or a one-row data frame ",
      "containing cpd_hat.",
      call. = FALSE
    )
  }
  cpd_hat <- selection$cpd_hat[[1L]]
  if (!is.numeric(cpd_hat) || anyNA(cpd_hat) ||
      any(cpd_hat != floor(cpd_hat)) ||
      any(cpd_hat < 1L | cpd_hat >= nrow(data))) {
    stop("selection$cpd_hat must contain integer locations in [1, n - 1].",
         call. = FALSE)
  }
  cpd_hat <- as.integer(cpd_hat)

  if (is.character(columns)) {
    if (is.null(colnames(data)) || anyNA(match(columns, colnames(data)))) {
      stop("columns contains names that are not present in data.",
           call. = FALSE)
    }
    columns <- match(columns, colnames(data))
  }
  if (!is.numeric(columns) || length(columns) == 0L || anyNA(columns) ||
      any(columns != floor(columns)) ||
      any(columns < 1L | columns > ncol(data))) {
    stop("columns must identify one or more columns of data.", call. = FALSE)
  }
  columns <- unique(as.integer(columns))
  values <- data[, columns, drop = FALSE]
  series_names <- colnames(values)
  if (is.null(series_names)) {
    series_names <- paste("series", columns)
    colnames(values) <- series_names
  }
  if (is.null(col)) {
    col <- seq_len(ncol(values))
  }
  if (is.null(ylab)) {
    ylab <- if (ncol(values) == 1L) series_names else "Value"
  }

  index <- seq_len(nrow(values))
  graphics::matplot(
    index, values, type = type, col = col, lty = lty,
    xlab = xlab, ylab = ylab, main = main, ...
  )
  if (length(cpd_hat) > 0L) {
    graphics::abline(
      v = cpd_hat + 0.5, col = cpd_col, lty = cpd_lty, lwd = cpd_lwd
    )
  }
  if (isTRUE(legend) && ncol(values) > 1L) {
    graphics::legend(
      "topright", legend = series_names, col = col, lty = lty, bty = "n"
    )
  }

  plotted <- data.frame(index = index, values, check.names = FALSE)
  attr(plotted, "cpd_hat") <- cpd_hat
  invisible(plotted)
}

#' Print a Reliever result
#'
#' Display any explicitly configured selection in compact form. If
#' the fit was created with \code{cpn_crit = "none"}, supplied WBS stopping
#' thresholds or PELT/OP penalties are shown directly. If neither a selection
#' rule nor such search values were supplied, or no loss output was marked as
#' eligible for the explicit fit-time criterion, the method points to the
#' complete candidate path in \code{x$cpd_path}. When several hyperparameter
#' values are printed, K has
#' been selected separately within each value; the hyperparameters themselves
#' have not been compared.
#'
#' @param x A result returned by \code{reliever()} or a related wrapper.
#' @param max_rows Maximum number of result rows printed. The default is 10;
#'   use \code{Inf} to print every model setting.
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#'
#' @return \code{x}, invisibly. Changepoint locations are printed in full;
#'   \code{x$summary} retains them as a list-column for programmatic use.
#' @seealso \code{\link{summary.reliever_result}()},
#'   \code{\link{plot.reliever_result}()},
#'   \code{\link{reliever}()}, \code{\link{reliever_generic}()}
#' @export
print.reliever_result <- function(x, max_rows = 10L, ...) {
  max_rows <- .reliever_validate_max_rows(max_rows)
  cat("Reliever changepoint result\n")
  context <- character()
  if (!is.null(x$settings$cpd_family)) {
    context <- c(context, paste("Family:", x$settings$cpd_family))
  }
  if (!is.null(x$settings$method)) {
    context <- c(context, paste("Method:", x$settings$method))
  }
  if (length(context)) {
    cat(paste(context, collapse = " | "), "\n", sep = "")
  }
  if (nrow(x$summary) == 0L) {
    .reliever_print_unselected_workflow()
    cat("Complete candidate paths are stored in x$cpd_path.\n")
  } else {
    print(x$summary, max_rows = max_rows, ...)
  }
  invisible(x)
}

#' Summarize a Reliever result
#'
#' Return the explicit selection stored in \code{object$summary}.
#' \itemize{
#'   \item With the ordinary \code{cpn_crit = "none"} default, K-indexed paths
#'   have an empty summary.
#'   \item An information-criterion fit reports its rule, estimated
#'   changepoints, and model setting. It selects K within each run, without
#'   comparing hyperparameters across runs.
#'   \item Supplied WBS thresholds or PELT/OP penalties retain their associated
#'   segmentations even when no K-selection criterion is requested.
#' }
#' Use \code{select_across_runs()} for joint selection across comparable paths.
#'
#' @param object A result returned by \code{reliever()} or a related wrapper.
#' @param ... Unused.
#'
#' @return The same compact data frame as \code{object$summary}, with one row
#'   per explicitly selected model setting or supplied search value.
#'   \code{K_hat} is the estimated number of changepoints and \code{cpd_hat} is
#'   a list-column containing their locations. A solution-path fit can return
#'   several rows, for example one per lambda, when its criterion selects K
#'   within each model setting rather than comparing settings. Printing that
#'   summary reports this distinction without adding another data-frame column.
#'   If the fitted object contains several statistically distinct loss-row
#'   types, \code{row_type} identifies the selected source, such as
#'   \code{"recv"}.
#'   A stored \code{hyper_value} is printed under the model-specific
#'   \code{hyper_name}, such as \code{lambda}, when that metadata is available.
#' @seealso \code{\link{reliever}()},
#'   \code{\link{plot.reliever_result}()}, and
#'   \code{\link{select_by_run}()}.
#' @export
summary.reliever_result <- function(object, ...) {
  object$summary
}
