# Hyperparameter path helpers --------------------------------------------------

.crossfit_hyper_count <- function(hyper_set) {
  if (is.data.frame(hyper_set) || is.matrix(hyper_set)) {
    nrow(hyper_set)
  } else {
    length(hyper_set)
  }
}

.crossfit_hyper_slice <- function(hyper_set, idx) {
  if (is.data.frame(hyper_set) || is.matrix(hyper_set)) {
    hyper_set[idx, , drop = FALSE]
  } else {
    hyper_set[idx]
  }
}

.crossfit_hyper_as_list <- function(hyper) {
  if (is.data.frame(hyper)) {
    out <- as.list(hyper[1L, , drop = FALSE])
  } else if (is.matrix(hyper)) {
    out <- as.list(hyper[1L, , drop = TRUE])
  } else if (is.list(hyper) && length(hyper) == 1L &&
             is.null(names(hyper)) && is.list(hyper[[1L]])) {
    out <- hyper[[1L]]
  } else if (is.list(hyper)) {
    out <- hyper
  } else {
    out <- list(value = hyper)
  }
  if (is.null(names(out))) {
    return(out)
  }
  out[!startsWith(names(out), ".")]
}

.crossfit_external_train_n <- function(data, caller) {
  n_train <- attr(data, ".reliever_external_train_n", exact = TRUE)
  if (is.null(n_train)) {
    return(NULL)
  }
  n_total <- .reliever_nobs(data)
  if (!is.numeric(n_train) || length(n_train) != 1L || is.na(n_train) ||
      !is.finite(n_train) || n_train != floor(n_train) ||
      n_train < 1L || n_train >= n_total) {
    stop(caller, " received an invalid external training-row marker.",
         call. = FALSE)
  }
  as.integer(n_train)
}

.crossfit_validate_loss_output_types <- function(loss_output_types) {
  allowed <- c("recv", "incv", "crossfit_homo_hyper")
  if (!is.character(loss_output_types) || length(loss_output_types) == 0L ||
      anyNA(loss_output_types) || any(!nzchar(loss_output_types))) {
    stop("loss_output_types must contain non-empty character values.",
         call. = FALSE)
  }
  unsupported <- setdiff(unique(loss_output_types), allowed)
  if (length(unsupported) > 0L) {
    stop(
      "Unsupported loss_output_types: ", paste(unsupported, collapse = ", "),
      ". Available values are: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!"recv" %in% loss_output_types) {
    stop("loss_output_types must include \"recv\".", call. = FALSE)
  }
  allowed[allowed %in% loss_output_types]
}

.crossfit_loss_output_meta <- function(hyper_set, hyper_name = NULL,
                                       loss_output_types = "recv") {
  loss_output_types <- .crossfit_validate_loss_output_types(
    loss_output_types
  )
  n_hyper <- .crossfit_hyper_count(hyper_set)
  hyper_values <- I(lapply(seq_len(n_hyper), function(i) {
    .crossfit_hyper_slice(hyper_set, i)
  }))
  row_type <- c(
    "recv",
    if ("incv" %in% loss_output_types) "incv",
    if ("crossfit_homo_hyper" %in% loss_output_types) {
      rep("crossfit_homo_hyper", n_hyper)
    }
  )
  hyper_id <- c(
    NA_integer_,
    if ("incv" %in% loss_output_types) NA_integer_,
    if ("crossfit_homo_hyper" %in% loss_output_types) seq_len(n_hyper)
  )
  meta <- data.frame(
    loss_output_id = seq_along(row_type),
    row_type = row_type,
    hyper_id = hyper_id
  )
  meta$hyper_value <- I(c(
    list(NULL),
    if ("incv" %in% loss_output_types) list(NULL),
    if ("crossfit_homo_hyper" %in% loss_output_types) hyper_values
  ))
  meta$default_selection <- meta$row_type == "recv"
  meta$loss_kind <- ifelse(
    meta$row_type %in% c("recv", "crossfit_homo_hyper"),
    "crossfit_loss",
    NA_character_
  )
  if (!is.null(hyper_name)) {
    if (!is.character(hyper_name) || length(hyper_name) != 1L ||
        is.na(hyper_name) || !nzchar(hyper_name)) {
      stop("hyper_name must be NULL or one non-empty character string.",
           call. = FALSE)
    }
    meta$hyper_name <- ifelse(
      is.na(meta$hyper_id), NA_character_, hyper_name
    )
  }
  rownames(meta) <- NULL
  meta
}

.crossfit_loss_outputs <- function(crossfit_homo_hyper,
                                   in_sample_hyper_loss = NULL,
                                   selection_loss = crossfit_homo_hyper,
                                   loss_output_types = "recv") {
  loss_output_types <- .crossfit_validate_loss_output_types(
    loss_output_types
  )
  n_hyper <- ncol(crossfit_homo_hyper)
  best_hyper_id <- which.min(colSums(selection_loss))[1L]
  loss_columns <- list(recv = crossfit_homo_hyper[, best_hyper_id])
  if ("incv" %in% loss_output_types) {
    if (is.null(in_sample_hyper_loss)) {
      stop("In-sample hyperparameter losses are required for incv output.",
           call. = FALSE)
    }
    loss_columns$incv <- in_sample_hyper_loss[, best_hyper_id]
  }
  if ("crossfit_homo_hyper" %in% loss_output_types) {
    loss_columns <- c(
      loss_columns,
      lapply(seq_len(n_hyper), function(j) crossfit_homo_hyper[, j])
    )
  }
  loss <- do.call(cbind, loss_columns)
  colnames(loss) <- c(
    "recv",
    if ("incv" %in% loss_output_types) "incv",
    if ("crossfit_homo_hyper" %in% loss_output_types) {
      rep("crossfit_homo_hyper", n_hyper)
    }
  )
  list(loss = loss, model = NULL)
}

.crossfit_foldid <- function(n_total, full_indices, l, r, l_end, r_end,
                       nfolds, fold_type, op_size,
                       fold_stable_const) {
  n_full <- r_end - l_end + 1L
  n_core <- r - l + 1L
  if (fold_type == "op") {
    rep(rep(seq_len(nfolds), each = op_size), length.out = n_full)
  } else if (fold_type == "blk") {
    idx_l_rel <- l - l_end + 1L
    idx_r_rel <- r - l_end + 1L
    foldid_full <- integer(n_full)
    foldid_core <- as.integer(cut(seq_len(n_core), breaks = nfolds,
                                  labels = FALSE))
    foldid_full[idx_l_rel:idx_r_rel] <- foldid_core
    if (idx_l_rel > 1L) {
      foldid_full[seq_len(idx_l_rel - 1L)] <- 1L
    }
    if (idx_r_rel < n_full) {
      foldid_full[(idx_r_rel + 1L):n_full] <- nfolds
    }
    foldid_full
  } else if (fold_type == "stable_blk") {
    if (fold_stable_const > 1) {
      exponent <- ceiling(log(n_total / n_core, base = fold_stable_const))
      n_stable <- floor(n_total / (fold_stable_const ^ exponent))
    } else {
      n_stable <- n_core
    }
    block_size <- floor(n_stable / nfolds)
    if (block_size < 1L) {
      stop(
        "stable_blk produced only ", n_stable,
        " stable positions for nfolds = ", nfolds,
        ". Reduce nfolds or fold_stable_const, or increase dm so fitted ",
        "intervals are wider.",
        call. = FALSE
      )
    }
    global_mid_point <- (1 + n_total) / 2
    block_indices <- floor((full_indices - global_mid_point) / block_size)
    (block_indices %% nfolds) + 1L
  } else {
    stop("fold_type must be one of \"op\", \"blk\", or \"stable_blk\".",
         call. = FALSE)
  }
}

.crossfit_loss_matrix <- function(loss, n_eval, n_hyper, caller) {
  loss <- .reliever_loss_matrix(loss)
  if (nrow(loss) != n_eval || ncol(loss) != n_hyper) {
    stop(
      caller, " must return a loss matrix with length(eval_id) rows and ",
      "length(hyper_set) columns.",
      call. = FALSE
    )
  }
  loss
}

.crossfit_logspace_add <- function(x, y) {
  largest <- pmax(x, y)
  largest + log(exp(x - largest) + exp(y - largest))
}

.crossfit_classifier_loss <- function(data, train_id, eval_id, hyper_set,
                                  l, r, l_end, r_end, n_total,
                                  prob_fun, eps, ...) {
  n_hyper <- .crossfit_hyper_count(hyper_set)
  y <- as.integer(seq_len(n_total) >= l & seq_len(n_total) <= r)
  target_prior <- (r - l + 1L) / n_total
  if (target_prior == 1) {
    return(matrix(0, length(eval_id), n_hyper))
  }
  train_prior <- mean(y[train_id])
  if (!is.finite(train_prior) || train_prior <= 0 || train_prior >= 1 ||
      target_prior <= 0 || target_prior > 1) {
    stop(
      "Classifier interval loss requires both inside- and outside-interval ",
      "training observations.",
      call. = FALSE
    )
  }
  loss_mat <- matrix(0, length(eval_id), n_hyper)
  for (j in seq_len(n_hyper)) {
    prob <- prob_fun(
      data = data,
      y = y,
      train_id = train_id,
      eval_id = eval_id,
      hyper = .crossfit_hyper_as_list(.crossfit_hyper_slice(hyper_set, j)),
      hyper_id = j,
      ...
    )
    prob <- as.numeric(prob)
    if (length(prob) != length(eval_id) || any(!is.finite(prob))) {
      stop("prob_fun must return one finite probability per eval_id.",
           call. = FALSE)
    }
    prob <- pmin(pmax(prob, eps), 1 - eps)
    log_density_odds <- stats::qlogis(prob) +
      log1p(-train_prior) - log(train_prior)
    log_mixture_ratio <- .crossfit_logspace_add(
      log(target_prior) + log_density_odds,
      log1p(-target_prior)
    )
    loss_mat[, j] <- log_mixture_ratio - log_density_odds
  }
  loss_mat
}

.reg_fun_crossfit_template_impl <- function(data, l, r, l_end, r_end,
                                      loss_fun, hyper_set, hyper_name,
                                      loss_output_types,
                                      crossfit_mode,
                                      nfolds, fold_type, op_size,
                                      buffer_lag, fold_stable_const,
                                      save_model, is_virtual_run, ...) {
  nfolds <- .reliever_validate_positive_integer(nfolds, "nfolds")
  if (nfolds < 2L) {
    stop("nfolds must be at least 2.", call. = FALSE)
  }
  op_size <- .reliever_validate_positive_integer(op_size, "op_size")
  buffer_lag <- .reliever_validate_positive_integer(
    buffer_lag, "buffer_lag", allow_zero = TRUE
  )
  if (!is.numeric(fold_stable_const) || length(fold_stable_const) != 1L ||
      is.na(fold_stable_const) || !is.finite(fold_stable_const) ||
      fold_stable_const < 1) {
    stop("fold_stable_const must be a finite number >= 1.", call. = FALSE)
  }
  n_hyper <- .crossfit_hyper_count(hyper_set)
  if (n_hyper < 1L) {
    stop("hyper_set must contain at least one value.", call. = FALSE)
  }
  loss_output_types <- .crossfit_validate_loss_output_types(
    loss_output_types
  )
  fold_type <- match.arg(fold_type, c("op", "blk", "stable_blk"))
  if (identical(crossfit_mode, "clf") && !identical(fold_type, "op")) {
    stop(
      "Classifier cross-fitting uses full-data folds and only supports ",
      "fold_type = \"op\".",
      call. = FALSE
    )
  }
  if (is_virtual_run) {
    loss_output_meta <- .crossfit_loss_output_meta(
      hyper_set, hyper_name, loss_output_types
    )
    return(list(
      n_loss_outputs = nrow(loss_output_meta),
      loss_output_meta = loss_output_meta
    ))
  }
  external_train_n <- .crossfit_external_train_n(
    data, "The generic crossfit template"
  )
  if (!is.null(external_train_n)) {
    template_name <- if (identical(crossfit_mode, "clf")) {
      "classifier crossfit template"
    } else {
      "generic crossfit template"
    }
    stop(
      "External hold-out evaluation is not supported by the ", template_name,
      ". No evaluation rows were used for fitting. Use a model-specific ",
      "external scorer; the built-in lasso-crossfit and KDE-NLL-crossfit ",
      "losses provide one.",
      call. = FALSE
    )
  }
  if (!is.function(loss_fun)) {
    stop("loss_fun must be a function.", call. = FALSE)
  }

  n_total <- .reliever_nobs(data)
  n_full <- r_end - l_end + 1L
  n_core <- r - l + 1L
  if (n_core < nfolds) {
    stop(
      "Fitted interval length ", n_core,
      " is smaller than nfolds = ", nfolds,
      ". Reduce nfolds or increase dm so every fitted interval can contain ",
      "all folds.",
      call. = FALSE
    )
  }
  full_indices <- l_end:r_end
  core_indices <- l:r
  row_idx <- function(idx) idx - l_end + 1L

  global_foldid <- NULL
  if (identical(crossfit_mode, "clf")) {
    global_indices <- seq_len(n_total)
    global_foldid <- .crossfit_foldid(
      n_total = n_total,
      full_indices = global_indices,
      l = 1L, r = n_total, l_end = 1L, r_end = n_total,
      nfolds = nfolds,
      fold_type = "op",
      op_size = op_size,
      fold_stable_const = fold_stable_const
    )
    foldid_full <- global_foldid[full_indices]
  } else {
    foldid_full <- .crossfit_foldid(
      n_total = n_total,
      full_indices = full_indices,
      l = l, r = r, l_end = l_end, r_end = r_end,
      nfolds = nfolds,
      fold_type = fold_type,
      op_size = op_size,
      fold_stable_const = fold_stable_const
    )
  }
  core_foldid <- foldid_full[row_idx(core_indices)]
  if (length(unique(core_foldid)) < 2L) {
    stop(
      "The fitted interval falls in only one fold, leaving no crossfit ",
      "training rows. Reduce op_size or use a wider interval.",
      call. = FALSE
    )
  }

  crossfit_homo_hyper <- matrix(0, n_full, n_hyper)
  for (f in seq_len(nfolds)) {
    id_te_full <- full_indices[foldid_full == f]
    id_tr_full <- full_indices[foldid_full != f]
    id_tr_proxy <- if (identical(crossfit_mode, "clf")) {
      which(global_foldid != f)
    } else {
      id_tr_full[id_tr_full >= l & id_tr_full <= r]
    }

    if (buffer_lag > 0) {
      if (n_core / nfolds <= 2 * buffer_lag) {
        buffer_lag_eff <- floor(n_core / (2 * nfolds))
      } else {
        buffer_lag_eff <- buffer_lag
      }
      forbidden_indices <- as.vector(
        outer(id_te_full, seq.int(-buffer_lag_eff, buffer_lag_eff), "+")
      )
      id_tr_proxy <- setdiff(id_tr_proxy, forbidden_indices)
    }
    if (length(id_tr_proxy) == 0L) {
      stop(
        "training set is empty. Check buffer_lag, nfolds, op_size, and ",
        "interval length.",
        call. = FALSE
      )
    }
    if (length(id_te_full) == 0L) {
      next
    }

    loss_fold <- loss_fun(
      data = data,
      train_id = id_tr_proxy,
      eval_id = id_te_full,
      hyper_set = hyper_set,
      l = l,
      r = r,
      l_end = l_end,
      r_end = r_end,
      n_total = n_total,
      ...
    )
    loss_fold <- .crossfit_loss_matrix(loss_fold, length(id_te_full), n_hyper,
                                 "loss_fun")
    crossfit_homo_hyper[row_idx(id_te_full), ] <- loss_fold
  }

  in_sample_hyper_loss <- NULL
  if ("incv" %in% loss_output_types) {
    full_train_id <- if (crossfit_mode == "clf") {
      seq_len(n_total)
    } else {
      core_indices
    }
    in_sample_hyper_loss <- loss_fun(
      data = data,
      train_id = full_train_id,
      eval_id = full_indices,
      hyper_set = hyper_set,
      l = l,
      r = r,
      l_end = l_end,
      r_end = r_end,
      n_total = n_total,
      ...
    )
    in_sample_hyper_loss <- .crossfit_loss_matrix(
      in_sample_hyper_loss, n_full, n_hyper, "loss_fun"
    )
  }

  core_rows <- row_idx(core_indices)
  .crossfit_loss_outputs(
    crossfit_homo_hyper,
    in_sample_hyper_loss,
    selection_loss = crossfit_homo_hyper[core_rows, , drop = FALSE],
    loss_output_types = loss_output_types
  )
}

# Main cross-fitting reg_fun templates ----------------------------------------

#' Generic cross-fitted interval-loss template
#'
#' Turn a model-specific observation-loss function into a complete crossfit
#' output set. Users implement only \code{loss_fun}: fit the model on
#' \code{train_id} and return the loss of every row in \code{eval_id}, once for
#' each candidate hyperparameter. The template constructs folds, combines
#' their losses, and returns the requested recycled-CV (ReCV) outputs. By
#' default only the interval-adaptive \code{recv} loss is returned.
#'
#' @param data A matrix-like object accepted by \code{loss_fun}, with
#'   observations in rows.
#' @param l,r Integer scalars giving the first and last rows used to fit the
#'   interval model.
#' @param l_end,r_end Integer scalars giving the first and last rows for which
#'   losses must be returned. This interval may extend beyond \code{l:r}
#'   because Reliever reuses one fitted model for nearby candidate intervals.
#' @param loss_fun A function called with \code{data}, \code{train_id},
#'   \code{eval_id}, \code{hyper_set}, the four interval endpoints, and
#'   \code{n_total}. It must return a numeric matrix with one row per
#'   \code{eval_id} and one column per hyperparameter. Every column must use
#'   the same observation-level scoring rule, units, normalization, and fold
#'   weighting. With one hyperparameter, a numeric vector is also accepted.
#' @param hyper_set A vector, list, matrix, or data frame of candidate
#'   hyperparameters. A vector or list uses one element per candidate; a matrix
#'   or data frame uses one row per candidate.
#' @param hyper_name \code{NULL} or a non-empty character scalar naming the
#'   tuning parameter, such as \code{"lambda"}. It is carried in the
#'   loss-output metadata so summaries and tuning-path plots can use a
#'   model-specific label. Use \code{NULL} when \code{hyper_set} is not
#'   described by one common name.
#' @param loss_output_types A character vector specifying the crossfit outputs.
#'   It is treated as an unordered set: input order does not affect either the
#'   calculation or output-column order. The default \code{"recv"} returns
#'   only the interval-adaptive cross-fitted loss. Use
#'   \code{c("recv", "incv")} to additionally return the matching in-sample
#'   loss chosen by the same interval-specific hyperparameter. Include
#'   \code{"crossfit_homo_hyper"} to additionally return one cross-fitted
#'   output per fixed hyperparameter for subsequent
#'   \code{select_across_runs()} selection. \code{"recv"} is required; an
#'   input that omits it produces an error. Returned columns always follow the
#'   order \code{recv}, \code{incv} when requested, then
#'   \code{crossfit_homo_hyper} in \code{hyper_set} order.
#' @param nfolds An integer scalar giving the number of folds constructed within
#'   each fitted interval. The interval length need not be divisible by
#'   \code{nfolds}.
#'   It must be an integer of at least 2 and cannot exceed the number of rows in
#'   the fitted interval.
#' @param fold_type A character scalar selecting the fold construction.
#'   \code{"op"} assigns interlaced,
#'   order-preserving fold labels; \code{"blk"} divides the current fitting
#'   interval into contiguous blocks; \code{"stable_blk"} uses contiguous
#'   blocks anchored to the full time axis so that nearby intervals have more
#'   stable fold labels. The default is \code{"op"}.
#' @param op_size A positive integer scalar giving, for
#'   \code{fold_type = "op"}, the number of consecutive
#'   observations given the same fold label before cycling to the next fold. It
#'   must be a positive integer and small enough that the fitted interval spans
#'   at least two fold labels.
#' @param buffer_lag A non-negative integer scalar giving the number of time
#'   indices excluded on each side of held-out
#'   observations when fitting a fold model. Use a positive value to reduce
#'   leakage under local temporal dependence; 0 applies no buffer. Within each
#'   fitted interval, the effective lag is capped at
#'   \code{floor(n_core / (2 * nfolds))} so that every fold remains trainable;
#'   it can therefore be 0 on short intervals.
#' @param fold_stable_const A numeric scalar of at least 1 giving, for
#'   \code{fold_type = "stable_blk"}, the geometric
#'   scale factor used to choose a globally anchored block width. The block
#'   width is based on the largest value in
#'   \eqn{n, n/c, n/c^2, \ldots} that does not exceed the fitted interval
#'   length, where \eqn{c} is \code{fold_stable_const}, and is then divided
#'   among \code{nfolds}. Thus values larger than 1 keep nearby interval lengths
#'   on the same global fold grid; the default 1 scales the width directly with
#'   the current interval.
#' @param save_model A logical scalar accepted for compatibility with the interval-loss
#'   protocol. This template returns losses only and does not retain fitted
#'   models.
#' @param is_virtual_run A logical scalar used as a metadata-query flag by
#'   \code{reliever_generic()}. When
#'   \code{TRUE}, return only the number and metadata of loss outputs without
#'   fitting any model.
#' @param ... Additional arguments passed to \code{loss_fun}.
#'
#' @return A named list following the custom \code{reg_fun} calling convention
#'   documented in the \emph{Writing a custom reg_fun} section of
#'   \code{\link{reliever_generic}()}. With
#'   \code{is_virtual_run = TRUE}, the list contains the integer scalar
#'   \code{n_loss_outputs} and the data frame \code{loss_output_meta}.
#'   Otherwise, it contains the numeric observation-by-output matrix
#'   \code{loss} and \code{model = NULL}. According to
#'   \code{loss_output_types}, the matrix can contain:
#'   \itemize{
#'     \item \code{recv}: within each fitted interval, choose the
#'     hyperparameter with the smallest cross-fitted loss and return that
#'     cross-fitted loss.
#'     \item \code{incv}: use the same interval-specific hyperparameter chosen
#'     by \code{recv}, but return its in-sample loss.
#'     \item \code{crossfit_homo_hyper}: one cross-fitted-loss output for every
#'     hyperparameter, keeping that value fixed rather than choosing it
#'     separately within each interval.
#'   }
#'   Selection follows two workflows:
#'   \itemize{
#'     \item Interval-adaptive ReCV uses \code{select_by_run()} with run type
#'     \code{"recv"} and criterion \code{"loss"}.
#'     \item Homogeneous-hyperparameter ReCV chooses K and one globally fixed
#'     hyperparameter from stored fixed-hyperparameter rows.
#'   }
#'   Apply \code{select_across_runs()} for the homogeneous workflow.
#'
#'   Its run type is \code{"crossfit_homo_hyper"}; request those rows first.
#' @seealso \code{\link{reliever_generic}()},
#'   \code{\link{reg_fun_clf_crossfit_template}()}, and
#'   \code{\link{select_by_run}()}.
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
#' loss_fun <- function(data, train_id, eval_id, hyper_set, ...) {
#'   center <- colMeans(data[train_id, , drop = FALSE])
#'   vapply(hyper_set, function(shrinkage) {
#'     fitted_center <- shrinkage * center
#'     rowMeans(
#'       sweep(data[eval_id, , drop = FALSE], 2, fitted_center, "-")^2
#'     )
#'   }, numeric(length(eval_id)))
#' }
#' out <- reg_fun_crossfit_template(
#'   data = x, l = 241, r = 660, l_end = 211, r_end = 690,
#'   loss_fun = loss_fun, hyper_set = c(0.8, 1),
#'   nfolds = 2
#' )
#' dim(out$loss)
reg_fun_crossfit_template <- function(data, l, r, l_end = l, r_end = r,
                                loss_fun, hyper_set, hyper_name = NULL,
                                loss_output_types = "recv",
                                nfolds = 5,
                                fold_type = c("op", "blk", "stable_blk"),
                                op_size = 1,
                                buffer_lag = 0,
                                fold_stable_const = 1,
                                save_model = FALSE,
                                is_virtual_run = FALSE,
                                ...) {
  .reg_fun_crossfit_template_impl(
    data = data,
    l = l,
    r = r,
    l_end = l_end,
    r_end = r_end,
    loss_fun = loss_fun,
    hyper_set = hyper_set,
    hyper_name = hyper_name,
    loss_output_types = loss_output_types,
    crossfit_mode = "normal",
    nfolds = nfolds,
    fold_type = match.arg(fold_type),
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = save_model,
    is_virtual_run = is_virtual_run,
    ...
  )
}

#' Cross-fitted probability-classifier interval-loss template
#'
#' Specialization of \code{reg_fun_crossfit_template()} for distribution-change
#' detection by binary classification. The template labels observations inside
#' the fitted interval as class 1 and observations outside it as class 0. User
#' code supplies only \code{prob_fun}, which fits the classifier and predicts a
#' class-1 probability for each requested evaluation row.
#'
#' @inheritParams reg_fun_crossfit_template
#' @param data A matrix-like object accepted by \code{prob_fun}, with
#'   observations in rows.
#' @param prob_fun A function called with \code{data}, binary labels \code{y},
#'   \code{train_id}, \code{eval_id}, one hyperparameter value \code{hyper}, and
#'   its index \code{hyper_id}. It must return one class-1 probability per
#'   evaluation row. Probabilities must be calibrated to the empirical class
#'   proportion in \code{y[train_id]}. If fitting uses class weights or
#'   class-balanced sampling, recalibrate the probabilities to that empirical
#'   training prior before returning them.
#' @param loss A character scalar selecting the segment loss computed from
#'   predicted probabilities. The
#'   supported value \code{"density_ratio_nll"} converts class probabilities
#'   to density-ratio negative log-likelihood. The conversion corrects for the
#'   inside-interval class proportion among the fold's training rows. If
#'   \eqn{p} is the predicted probability, \eqn{\pi_t} that training
#'   proportion, and \eqn{\pi=(r-l+1)/n}, it sets
#'   \eqn{q=\{p/(1-p)\}\{(1-\pi_t)/\pi_t\}} and returns
#'   \eqn{-\log[q/\{\pi q+(1-\pi)\}]}.
#' @param eps Probability-clipping constant in
#'   \code{[.Machine$double.eps, 0.5)}. Clipping precedes logarithms and keeps
#'   both endpoints distinct from zero and one in double precision.
#' @param fold_type A character scalar selecting fold construction. Classifier
#'   cross-fitting fits against
#'   the full data and supports only the globally anchored, interlaced
#'   \code{"op"} construction.
#' @param nfolds An integer scalar of at least 2 giving the number of global
#'   \code{"op"} folds. The fitted interval must contain at least
#'   \code{nfolds} observations and span at least two fold labels.
#' @param fold_stable_const Unused for classifier cross-fitting, which supports
#'   only \code{"op"} folds. Retained for argument compatibility.
#' @param ... Additional arguments passed to \code{prob_fun}.
#'
#' @return A named list. With \code{is_virtual_run = TRUE}, it contains an
#'   integer scalar \code{n_loss_outputs} and a data frame
#'   \code{loss_output_meta}. Otherwise, it contains the requested numeric
#'   observation-by-output \code{loss} matrix and \code{model = NULL}. See
#'   \code{\link{reg_fun_crossfit_template}()} for the full return convention
#'   and column definitions.
#'
#' @seealso \code{\link{reg_fun_crossfit_template}()},
#'   \code{\link{reg_fun_ranger_crossfit}()}, and
#'   \code{\link{reg_fun_mlp_crossfit}()}.
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
#' prob_fun <- function(data, y, train_id, eval_id, hyper, ...) {
#'   rep(pmin(pmax(mean(y[train_id]) + hyper$offset, 0.01), 0.99),
#'       length(eval_id))
#' }
#' out <- reg_fun_clf_crossfit_template(
#'   data = x, l = 241, r = 660, l_end = 211, r_end = 690,
#'   prob_fun = prob_fun, hyper_set = data.frame(offset = c(0, 0.05)),
#'   nfolds = 2
#' )
#' dim(out$loss)
reg_fun_clf_crossfit_template <- function(data, l, r, l_end = l, r_end = r,
                                    prob_fun, hyper_set,
                                    hyper_name = NULL,
                                    loss_output_types = "recv",
                                    nfolds = 5,
                                    fold_type = "op",
                                    op_size = 1,
                                    buffer_lag = 0,
                                    fold_stable_const = 1,
                                    save_model = FALSE,
                                    is_virtual_run = FALSE,
                                    loss = "density_ratio_nll",
                                    eps = 1e-12,
                                    ...) {
  if (!is.function(prob_fun)) {
    stop("prob_fun must be a function.", call. = FALSE)
  }
  loss <- match.arg(loss, "density_ratio_nll")
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) ||
      !is.finite(eps) || eps < .Machine$double.eps || eps >= 0.5) {
    stop(
      "eps must be a finite number in [.Machine$double.eps, 0.5).",
      call. = FALSE
    )
  }

  out <- .reg_fun_crossfit_template_impl(
    data = data,
    l = l,
    r = r,
    l_end = l_end,
    r_end = r_end,
    loss_fun = .crossfit_classifier_loss,
    hyper_set = hyper_set,
    hyper_name = hyper_name,
    loss_output_types = loss_output_types,
    crossfit_mode = "clf",
    nfolds = nfolds,
    fold_type = fold_type,
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = save_model,
    is_virtual_run = is_virtual_run,
    prob_fun = prob_fun,
    eps = eps,
    ...
  )
  if (is_virtual_run) {
    missing_kind <- is.na(out$loss_output_meta$loss_kind)
    out$loss_output_meta$loss_kind[missing_kind] <-
      "negative_log_likelihood"
  }
  out
}

# glmnet/lasso crossfit reg_fun -----------------------------------------------

.loss_lasso_crossfit <- function(data, train_id, eval_id, hyper_set,
                           family = "gaussian", thresh = 1e-7,
                           intercept = FALSE, standardize = FALSE, ...) {
  model <- glmnet::glmnet(
    data[train_id, -1, drop = FALSE],
    data[train_id, 1],
    family = family,
    lambda = hyper_set / sqrt(length(train_id)),
    intercept = intercept,
    standardize = standardize,
    control = list(thresh = thresh)
  )
  .reliever_glmnet_loss(
    model,
    data[eval_id, -1, drop = FALSE],
    data[eval_id, 1],
    family
  )
}

.lasso_crossfit_external <- function(data, l, r, l_end, r_end,
                                     lam_set, family, thresh,
                                     nfolds, fold_type, op_size,
                                     buffer_lag, fold_stable_const,
                                     loss_output_types) {
  n_train <- .crossfit_external_train_n(
    data, "External lasso-crossfit evaluation"
  )
  if (is.null(n_train)) {
    stop("External lasso-crossfit evaluation requires a training-row marker.",
         call. = FALSE)
  }
  if (!is.matrix(data) || ncol(data) < 3L ||
      l < 1L || r > n_train || l > r ||
      l_end <= n_train || r_end > nrow(data) || l_end > r_end) {
    stop("External lasso-crossfit evaluation has invalid data or endpoints.",
         call. = FALSE)
  }

  train_data <- data[seq_len(n_train), , drop = FALSE]
  training_fit <- reg_fun_crossfit_template(
    data = train_data,
    l = l,
    r = r,
    l_end = l,
    r_end = r,
    loss_fun = .loss_lasso_crossfit,
    hyper_set = lam_set,
    hyper_name = "lambda",
    loss_output_types = c("recv", "crossfit_homo_hyper"),
    nfolds = nfolds,
    fold_type = fold_type,
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = FALSE,
    is_virtual_run = FALSE,
    family = family,
    thresh = thresh
  )

  external_fixed <- .loss_lasso_crossfit(
    data = data,
    train_id = l:r,
    eval_id = l_end:r_end,
    hyper_set = lam_set,
    family = family,
    thresh = thresh
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

#' Cross-fitted lasso interval loss
#'
#' Fit a glmnet lasso path within each interval fold. The interval-adaptive
#' \code{recv} column chooses a normalized lambda setting in each fitted
#' interval. When requested, the
#' \code{crossfit_homo_hyper} columns support joint selection of one normalized
#' setting and K. Most users should call \code{reliever()} with
#' \code{cpd_family = "lasso_crossfit"}.
#' The built-in model uses \code{intercept = FALSE} and
#' \code{standardize = FALSE}; predictors should already be on their intended
#' scale.
#'
#' @param data Numeric matrix with the response in the first column and at least
#'   two predictor columns, as required by \code{glmnet}.
#' @inheritParams reg_fun_crossfit_template
#' @param lam_set A numeric vector of finite positive lambda values on a
#'   sample-size-independent
#'   scale. A fold containing \eqn{m} training observations passes
#'   \code{lam_set / sqrt(m)} to glmnet.
#' @param family A character scalar selecting the response family:
#'   \code{"gaussian"}, \code{"binomial"}, or
#'   \code{"poisson"}. Gaussian returns squared prediction loss; binomial and
#'   Poisson return negative log-likelihood loss. Binomial responses must be
#'   0 or 1; Poisson responses may contain any finite non-negative values.
#' @param thresh Positive convergence threshold passed to \code{glmnet}.
#'
#' @return A named list. With \code{is_virtual_run = TRUE}, it contains the
#'   integer scalar \code{n_loss_outputs} and data frame
#'   \code{loss_output_meta}. Otherwise, it contains the requested numeric
#'   observation-by-output \code{loss} matrix and \code{model = NULL}. See
#'   \code{\link{reg_fun_crossfit_template}()} for the full return convention
#'   and column definitions.
#'
#' @seealso \code{\link{reliever_lasso_crossfit}()},
#'   \code{\link{reg_fun_lasso_solpath}()},
#'   \code{\link{reg_fun_crossfit_template}()}
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
#' out <- reg_fun_lasso_crossfit(
#'   data = data, l = 241, r = 660, l_end = 211, r_end = 690,
#'   lam_set = lam_set, nfolds = 2
#' )
#' colSums(out$loss)
reg_fun_lasso_crossfit <- function(data, l, r, l_end = l, r_end = r,
                             nfolds = 5,
                             fold_type = c("op", "blk", "stable_blk"),
                             op_size = 1,
                             buffer_lag = 0,
                             fold_stable_const = 1,
                             loss_output_types = "recv",
                             save_model = FALSE, is_virtual_run = FALSE,
                             lam_set = NULL,
                             family = "gaussian",
                             thresh = 1e-7) {
  external_train_n <- attr(
    data, ".reliever_external_train_n", exact = TRUE
  )
  data <- .reliever_lasso_data_matrix(data)
  family <- .reliever_validate_glmnet_family(family)
  involved <- sort(unique(c(seq.int(l, r), seq.int(l_end, r_end))))
  .reliever_validate_lasso_response(data[involved, 1L], family)
  if (!is.null(external_train_n)) {
    attr(data, ".reliever_external_train_n") <- external_train_n
  }
  if (is.null(lam_set)) {
    stop(
      paste0(
        "lam_set is required by reg_fun_lasso_crossfit(); use ",
        "reliever(..., cpd_family = \"lasso_crossfit\") for an automatic path."
      ),
      call. = FALSE
    )
  }
  lam_set <- .reliever_validate_lam_set(lam_set)
  fold_type <- match.arg(fold_type)
  if (is_virtual_run) {
    meta <- .crossfit_loss_output_meta(
      lam_set, "lambda", loss_output_types
    )
    is_crossfit <- meta$row_type %in% c("recv", "crossfit_homo_hyper")
    meta$loss_kind <- ifelse(
      is_crossfit,
      "crossfit_loss",
      if (identical(family, "gaussian")) {
        "rss"
      } else {
        "negative_log_likelihood"
      }
    )
    return(list(
      n_loss_outputs = nrow(meta),
      loss_output_meta = meta
    ))
  }
  if (!is_virtual_run && !is.null(external_train_n)) {
    return(.lasso_crossfit_external(
      data = data,
      l = l,
      r = r,
      l_end = l_end,
      r_end = r_end,
      lam_set = lam_set,
      family = family,
      thresh = thresh,
      nfolds = nfolds,
      fold_type = fold_type,
      op_size = op_size,
      buffer_lag = buffer_lag,
      fold_stable_const = fold_stable_const,
      loss_output_types = loss_output_types
    ))
  }
  reg_fun_crossfit_template(
    data = data,
    l = l,
    r = r,
    l_end = l_end,
    r_end = r_end,
    loss_fun = .loss_lasso_crossfit,
    hyper_set = lam_set,
    hyper_name = "lambda",
    loss_output_types = loss_output_types,
    nfolds = nfolds,
    fold_type = fold_type,
    op_size = op_size,
    buffer_lag = buffer_lag,
    fold_stable_const = fold_stable_const,
    save_model = save_model,
    is_virtual_run = is_virtual_run,
    family = family,
    thresh = thresh
  )
}
