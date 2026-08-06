.reliever_cache_profile_backend <- function(cache_profile) {
  profile_backend <- cache_profile$backend
  if (is.null(profile_backend) && !is.null(cache_profile$summary)) {
    profile_backend <- cache_profile$summary$backend
  }
  profile_backend
}

.reliever_detail_mode <- function(detail) {
  if (isTRUE(detail)) {
    return("cache")
  }
  if (identical(detail, FALSE) || is.null(detail)) {
    return("none")
  }
  if (is.character(detail) && length(detail) == 1L && !is.na(detail)) {
    mode <- match.arg(tolower(detail), c("none", "cache"))
    return(mode)
  }
  stop("detail must be TRUE, FALSE, \"cache\", or \"none\".",
       call. = FALSE)
}

.reliever_validate_cache_profile_backend <- function(cache_profile, cache_backend) {
  if (is.null(cache_profile)) {
    return(invisible(NULL))
  }
  if (!is.list(cache_profile)) {
    stop("cache_profile must be a reliever cache_profile list.", call. = FALSE)
  }
  profile_backend <- .reliever_cache_profile_backend(cache_profile)
  if (!identical(profile_backend, cache_backend)) {
    stop(
      sprintf(
        "cache_profile backend is %s, but this run uses cache_backend = %s.",
        dQuote(profile_backend),
        dQuote(cache_backend)
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.reliever_validate_arma_double_matrix <- function(x, name) {
  if (is.null(x)) {
    stop(sprintf("%s must be a double matrix, not NULL.", name),
         call. = FALSE)
  }
  if (!is.matrix(x) || !identical(typeof(x), "double")) {
    stop(
      sprintf(
        "%s must be a double matrix so it can be passed to arma::mat without copying. Use numeric/double matrices such as matrix(Inf, ...), matrix(0, ...), or matrix(NA_real_, ...); integer/logical matrices such as matrix(1:10, ...), matrix(TRUE, ...), or matrix(NA, ...) are not allowed.",
        name
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.reliever_backend_cache_from_profile <- function(cache_profile, setup,
                                                 loss_context) {
  profile_backend <- if (is.null(cache_profile)) {
    NULL
  } else {
    .reliever_cache_profile_backend(cache_profile)
  }
  if (setup$is_full && identical(setup$cache_backend, "by_cost_mat") &&
      identical(profile_backend, "by_loss_block")) {
    warning(
      "Ignoring by_loss_block cache_profile because full-search runs use ",
      "cache_backend = \"by_cost_mat\". Omit cache_profile or supply a ",
      "profile created with cache_backend = \"by_cost_mat\".",
      call. = FALSE
    )
    cache_profile <- NULL
  }
  .reliever_validate_cache_profile_backend(cache_profile, setup$cache_backend)
  if (is.null(cache_profile)) {
    return(list())
  }
  if (is.null(cache_profile$loss_context) ||
      !identical(cache_profile$loss_context, loss_context)) {
    stop(
      paste0(
        "cache_profile was created for different data, reg_fun, dc_grid, ",
        "or loss-function arguments. Recompute it for the current loss context."
      ),
      call. = FALSE
    )
  }
  switch(
    setup$cache_backend,
    by_cost_mat = list(
      cost_mat = .reliever_cost_mat_from_cache_profile(
        cache_profile = cache_profile,
        n_loss_outputs = setup$n_loss_outputs,
        n = setup$n
      )
    ),
    by_loss_block = list(
      loss_block_cache = .reliever_loss_block_cache_from_profile(
        cache_profile = cache_profile,
        setup = setup
      )
    )
  )
}

.reliever_cost_mat_from_cache_profile <- function(cache_profile,
                                                  n_loss_outputs, n) {
  if (is.null(cache_profile)) {
    return(NULL)
  }
  .reliever_validate_cache_profile_backend(cache_profile, "by_cost_mat")
  cost_mat <- cache_profile$objects$cost_mat
  if (is.null(cost_mat)) {
    stop(
      "cache_profile must contain a double matrix at cache_profile$objects$cost_mat. Rerun with detail = TRUE before reusing it.",
      call. = FALSE
    )
  }
  .reliever_validate_arma_double_matrix(
    cost_mat, "cache_profile$objects$cost_mat"
  )
  expected_dim <- c(n_loss_outputs, 1 + n * (n + 1) / 2)
  if (length(dim(cost_mat)) != 2L || any(dim(cost_mat) != expected_dim)) {
    stop(
      sprintf(
        "cache_profile$objects$cost_mat has dimension %s, but this reliever run expects %s.",
        paste(dim(cost_mat), collapse = " x "),
        paste(expected_dim, collapse = " x ")
      ),
      call. = FALSE
    )
  }
  cost_mat
}

.reliever_loss_block_cache_from_profile <- function(cache_profile, setup) {
  .reliever_validate_cache_profile_backend(cache_profile, "by_loss_block")
  cache_state <- cache_profile$objects$loss_block_cache
  if (is.null(cache_state)) {
    stop(
      "cache_profile must contain cache_profile$objects$loss_block_cache. Rerun with detail = TRUE before reusing it.",
      call. = FALSE
    )
  }
  if (!identical(as.integer(cache_state$n), as.integer(setup$n))) {
    stop("cache_profile$objects$loss_block_cache was built for a different n.",
         call. = FALSE)
  }
  cached_int_set <- cache_profile$objects$int_set
  if (!is.null(cached_int_set) &&
      (!identical(cached_int_set$miss_cover_len, setup$int_set$miss_cover_len) ||
       !identical(cached_int_set$int_len, setup$int_set$int_len) ||
       !identical(cached_int_set$layer_point, setup$int_set$layer_point) ||
       !identical(cached_int_set$int_eps, setup$int_set$int_eps))) {
    stop("by_loss_block cache_profile reuse requires the same Reliever interval set.",
         call. = FALSE)
  }
  .reliever_validate_loss_block_cache_for_arma(cache_state)
  cache_state
}

.reliever_validate_loss_block_cache_for_arma <- function(cache_state) {
  if (is.null(cache_state$blocks)) {
    return(invisible(NULL))
  }
  if (!is.list(cache_state$blocks)) {
    stop("cache_profile$objects$loss_block_cache$blocks must be a list.",
         call. = FALSE)
  }
  for (i in seq_along(cache_state$blocks)) {
    block <- cache_state$blocks[[i]]
    if (!is.list(block) || is.null(block$loss)) {
      stop(
        sprintf(
          "cache_profile$objects$loss_block_cache$blocks[[%d]] must contain a double matrix named loss.",
          i
        ),
        call. = FALSE
      )
    }
    .reliever_validate_arma_double_matrix(
      block$loss,
      sprintf(
        "cache_profile$objects$loss_block_cache$blocks[[%d]]$loss",
        i
      )
    )
  }
  invisible(NULL)
}

.reliever_by_loss_block_cache_profile <- function(cache_profile, int_set,
                                                  loss_context = NULL) {
  out <- list(
    backend = "by_loss_block",
    objects = list(
      int_set = int_set,
      loss_block_cache = cache_profile$cache_state
    )
  )
  if (!is.null(loss_context)) {
    out$loss_context <- loss_context
  }
  out
}

.reliever_by_cost_mat_cache_profile <- function(cost_mat, int_set,
                                                loss_context = NULL) {
  out <- list(
    backend = "by_cost_mat",
    objects = list(
      cost_mat = cost_mat,
      int_set = int_set
    )
  )
  if (!is.null(loss_context)) {
    out$loss_context <- loss_context
  }
  out
}
