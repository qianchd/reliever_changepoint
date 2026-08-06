.reliever_validate_loss_output_meta <- function(loss_output_meta,
                                                n_loss_outputs) {
  if (is.null(loss_output_meta)) {
    return(NULL)
  }
  loss_output_meta <- as.data.frame(loss_output_meta)
  if (nrow(loss_output_meta) != n_loss_outputs) {
    stop("loss_output_meta must have one row per loss output.", call. = FALSE)
  }
  if (!"loss_output_id" %in% names(loss_output_meta)) {
    loss_output_meta$loss_output_id <- seq_len(n_loss_outputs)
  }
  loss_output_meta$loss_output_id <- .reliever_validate_positive_integer_vector(
    loss_output_meta$loss_output_id, "loss_output_meta$loss_output_id"
  )
  if (anyDuplicated(loss_output_meta$loss_output_id) ||
      !setequal(loss_output_meta$loss_output_id, seq_len(n_loss_outputs))) {
    stop(
      "loss_output_meta$loss_output_id must contain each loss-output index exactly once.",
      call. = FALSE
    )
  }
  if ("default_selection" %in% names(loss_output_meta) &&
      (!is.logical(loss_output_meta$default_selection) ||
       anyNA(loss_output_meta$default_selection))) {
    stop("loss_output_meta$default_selection must contain TRUE or FALSE.",
         call. = FALSE)
  }
  if ("loss_kind" %in% names(loss_output_meta)) {
    loss_kind <- loss_output_meta$loss_kind
    if (!is.character(loss_kind) ||
        any(!is.na(loss_kind) & !nzchar(loss_kind))) {
      stop(
        "loss_output_meta$loss_kind must contain non-empty strings or NA.",
        call. = FALSE
      )
    }
  }
  rownames(loss_output_meta) <- NULL
  loss_output_meta
}

.reliever_run_meta <- function(loss_output_meta, run_cpd_ids) {
  if (is.null(loss_output_meta)) {
    return(data.frame(
      run_id = seq_along(run_cpd_ids),
      loss_output_id = as.integer(run_cpd_ids)
    ))
  }
  id <- match(run_cpd_ids, loss_output_meta$loss_output_id)
  out <- loss_output_meta[id, , drop = FALSE]
  out$run_id <- seq_along(run_cpd_ids)
  readable_first <- c("run_id", "row_type", "hyper_name", "hyper_value")
  out <- out[, c(
    intersect(readable_first, names(out)),
    setdiff(names(out), readable_first)
  ), drop = FALSE]
  rownames(out) <- NULL
  out
}

.reliever_run_meta_structure <- function(run_meta,
                                         compare_hyper_value = TRUE) {
  run_meta <- as.data.frame(run_meta)
  structural_columns <- c(
    "run_id", "loss_output_id", "hyper_id",
    if (compare_hyper_value) "hyper_value"
  )
  out <- run_meta[
    , intersect(structural_columns, names(run_meta)), drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

.reliever_validate_run_meta_structure <- function(actual, expected,
                                                  error_message,
                                                  compare_hyper_value = TRUE) {
  if (!identical(
    .reliever_run_meta_structure(actual, compare_hyper_value),
    .reliever_run_meta_structure(expected, compare_hyper_value)
  )) {
    stop(error_message, call. = FALSE)
  }

  if ("loss_kind" %in% names(actual) &&
      "loss_kind" %in% names(expected)) {
    actual_kind <- as.character(actual$loss_kind)
    expected_kind <- as.character(expected$loss_kind)
    comparable <- !is.na(actual_kind) & nzchar(actual_kind) &
      !is.na(expected_kind) & nzchar(expected_kind)
    if (any(comparable & actual_kind != expected_kind)) {
      warning(warningCondition(
        paste0(
          "Corresponding runs declare different loss_kind values; losses ",
          "may not be comparable. Evaluation will continue."
        ),
        call = NULL,
        class = "reliever_loss_kind_comparability_warning"
      ))
    }
  }
  invisible(actual)
}

.reliever_reg_fun_virtual_info <- function(reg_fun, data, para_list) {
  if (!is.function(reg_fun)) {
    stop("reg_fun must be a function.", call. = FALSE)
  }
  formal_names <- names(formals(reg_fun))
  if (!"..." %in% formal_names) {
    required <- c(
      "data", "l", "r", "l_end", "r_end", "save_model",
      "is_virtual_run"
    )
    missing_formals <- setdiff(required, formal_names)
    if (length(missing_formals) > 0L) {
      stop(
        "reg_fun must accept ", paste(required, collapse = ", "),
        " (or ...). Missing: ", paste(missing_formals, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  virtual_args <- c(
    list(data = data, l = 1L, r = 1L, is_virtual_run = TRUE),
    para_list
  )
  virtual_res <- do.call(reg_fun, virtual_args)
  if (is.list(virtual_res) && !is.null(virtual_res$n_loss_outputs)) {
    n_loss_outputs <- virtual_res$n_loss_outputs
    loss_output_meta <- virtual_res$loss_output_meta
  } else {
    n_loss_outputs <- virtual_res
    loss_output_meta <- NULL
  }
  if (!is.numeric(n_loss_outputs) || length(n_loss_outputs) != 1 ||
      is.na(n_loss_outputs) || n_loss_outputs != floor(n_loss_outputs) ||
      n_loss_outputs < 1) {
    stop("reg_fun(..., is_virtual_run = TRUE) must return the number of loss outputs or a list containing n_loss_outputs.",
         call. = FALSE)
  }
  n_loss_outputs <- as.integer(n_loss_outputs)
  list(
    n_loss_outputs = n_loss_outputs,
    loss_output_meta = .reliever_validate_loss_output_meta(
      loss_output_meta, n_loss_outputs
    )
  )
}

.reliever_search_dm_from_grid <- function(dm, dc_grid) {
  if (is.null(dc_grid)) {
    return(dm)
  }

  grid_boundaries <- c(0L, dc_grid)
  search_n <- length(dc_grid)
  for (k in seq_len(search_n)) {
    span <- grid_boundaries[(k + 1L):(search_n + 1L)] -
      grid_boundaries[seq_len(search_n + 1L - k)]
    if (all(span >= dm)) {
      return(as.integer(k))
    }
  }
  stop("dc_grid cannot support the requested minimum segment length dm.",
       call. = FALSE)
}

.reliever_validate_wbs_stop_crit <- function(wbs_stop_crit, method) {
  if (is.null(wbs_stop_crit)) {
    return(NULL)
  }
  if (!method %in% c("WBS", "WBS_recursive", "SeedBS", "BS")) {
    stop("wbs_stop_crit can only be used with method = \"WBS\", \"WBS_recursive\", \"SeedBS\", or \"BS\".",
         call. = FALSE)
  }
  if (!is.numeric(wbs_stop_crit) || length(wbs_stop_crit) == 0L ||
      anyNA(wbs_stop_crit) || any(!is.finite(wbs_stop_crit))) {
    stop("wbs_stop_crit must be a finite numeric vector.", call. = FALSE)
  }
  as.numeric(wbs_stop_crit)
}

.reliever_prepare_search_setup <- function(data, cpn_max, dm, pen_val,
                                           cov_rate, method, M, prune_value,
                                           dc_grid, cache_backend,
                                           wbs_stop_crit = NULL) {
  method <- match.arg(
    method, c("SN", "WBS", "WBS_recursive", "SeedBS", "BS", "PELT", "OP")
  )
  cache_backend <- match.arg(cache_backend, c("by_loss_block", "by_cost_mat"))
  wbs_stop_crit <- .reliever_validate_wbs_stop_crit(wbs_stop_crit, method)
  cpn_max <- .reliever_validate_positive_integer(
    cpn_max, "cpn_max", allow_zero = TRUE
  )
  dm <- .reliever_validate_positive_integer(dm, "dm")
  M <- .reliever_validate_positive_integer(M, "M", allow_zero = TRUE)
  if (!is.numeric(cov_rate) || length(cov_rate) != 1 || is.na(cov_rate) ||
      cov_rate <= 0 || cov_rate > 1) {
    stop("cov_rate must be a single number in (0, 1].", call. = FALSE)
  }
  if (method == "BS") {
    M <- 0
  }
  if (method %in% c("PELT", "OP")) {
    if (!is.numeric(pen_val) || is.complex(pen_val) ||
        !is.null(dim(pen_val)) || length(pen_val) == 0L ||
        anyNA(pen_val) || any(!is.finite(pen_val)) || any(pen_val < 0)) {
      stop("pen_val must be a non-empty vector of finite, non-negative numbers for PELT and OP.",
           call. = FALSE)
    }
    pen_val <- as.numeric(pen_val)
  }
  if (method == "PELT" &&
      (!is.numeric(prune_value) || is.complex(prune_value) ||
       length(prune_value) != 1L ||
       is.na(prune_value) || !is.finite(prune_value))) {
    stop("prune_value must be a finite numeric scalar for PELT.",
         call. = FALSE)
  }
  if (method == "OP") {
    prune_value <- -Inf
  }

  n_original <- .reliever_nobs(data)
  n <- n_original
  dm_original <- dm
  if (!is.null(dc_grid)) {
    dc_grid_int <- suppressWarnings(as.integer(dc_grid))
    if (!is.numeric(dc_grid) || length(dc_grid) == 0L ||
        anyNA(dc_grid) || anyNA(dc_grid_int) ||
        any(dc_grid != dc_grid_int)) {
      stop("dc_grid must contain integer sample indices.", call. = FALSE)
    }
    dc_grid <- sort(dc_grid_int)
    if (anyDuplicated(dc_grid) || any(dc_grid < 1L) ||
        tail(dc_grid, 1L) != n) {
      stop("dc_grid must contain unique integer indices in [1, n] and end with n.",
           call. = FALSE)
    }
    n <- length(dc_grid)
    dm <- .reliever_search_dm_from_grid(dm_original, dc_grid)
  }
  if (dm > n) {
    stop("n must be at least dm for changepoint detection.", call. = FALSE)
  }

  if (cov_rate >= (n - 1) / n) {
    cov_rate <- 1
  }
  is_full <- cov_rate == 1
  if (is_full && identical(cache_backend, "by_loss_block")) {
    warning(warningCondition(
      paste0(
        "cache_backend = \"by_loss_block\" is not used for full-search runs; ",
        "using cache_backend = \"by_cost_mat\" instead. Set ",
        "cache_backend = \"by_cost_mat\" explicitly to suppress this warning."
      ),
      call = NULL,
      class = "reliever_cache_backend_warning"
    ))
    cache_backend <- "by_cost_mat"
  }
  int_set <- if (!is_full) create_relief_itv(n, cov_rate, dm) else matrix(NA, 0, 0)
  if ((cpn_max + 1) * dm > n) {
    cpn_max <- floor(n / dm) - 1
  }
  if (method %in% c("PELT", "OP")) {
    cpn_max <- length(pen_val)
  }

  list(
    method = method,
    cache_backend = cache_backend,
    cpn_max = cpn_max,
    dm = dm,
    dm_original = dm_original,
    M = M,
    prune_value = prune_value,
    n_original = n_original,
    n = n,
    dc_grid = dc_grid,
    cov_rate = cov_rate,
    is_full = is_full,
    int_set = int_set,
    wbs_stop_crit = wbs_stop_crit
  )
}

.reliever_prepare <- function(data, cpn_max, dm, pen_val, cov_rate, reg_fun,
                              method, M, prune_value, run_cpd_ids, dc_grid,
                              cache_backend, wbs_stop_crit, para_list) {
  setup <- .reliever_prepare_search_setup(
    data = data,
    cpn_max = cpn_max,
    dm = dm,
    pen_val = pen_val,
    cov_rate = cov_rate,
    method = method,
    M = M,
    prune_value = prune_value,
    dc_grid = dc_grid,
    cache_backend = cache_backend,
    wbs_stop_crit = wbs_stop_crit
  )

  virtual_info <- .reliever_reg_fun_virtual_info(reg_fun, data, para_list)
  n_loss_outputs <- virtual_info$n_loss_outputs

  if (is.null(run_cpd_ids)) {
    run_cpd_ids <- seq_len(n_loss_outputs)
  }
  run_cpd_ids <- unique(.reliever_validate_positive_integer_vector(
    run_cpd_ids, "run_cpd_ids"
  ))
  if (max(run_cpd_ids) > n_loss_outputs) {
    stop(
      "run_cpd_ids contains an index larger than the number of loss outputs.",
      call. = FALSE
    )
  }

  c(
    setup,
    list(
      n_loss_outputs = n_loss_outputs,
      run_cpd_ids = run_cpd_ids,
      run_meta = .reliever_run_meta(
        virtual_info$loss_output_meta, run_cpd_ids
      )
    )
  )
}
