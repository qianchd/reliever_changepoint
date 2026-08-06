# Search interval construction ------------------------------------------------

.reliever_search_intervals <- function(method, n, dm, M, wbs_seed) {
  if (method %in% c("WBS", "WBS_recursive", "SeedBS", "BS") &&
      n < 2L * dm) {
    return(matrix(integer(), 0L, 2L))
  }
  if (method %in% c("WBS", "WBS_recursive", "BS")) {
    create_wbs_itv(n, 2 * dm, M, wbs_seed)
  } else if (method == "SeedBS") {
    create_seed_itv(n, 2 * dm)
  } else {
    matrix(integer(), 0L, 2L)
  }
}

# Backend result formatting ----------------------------------------------------

.reliever_run_candidate_indices <- function(backend_res, run_id) {
  start <- as.integer(backend_res$run_start[run_id])
  len <- as.integer(backend_res$run_len[run_id])
  if (len == 0L) {
    return(integer())
  }
  start + seq_len(len)
}

.reliever_decode_cpd <- function(backend_res, candidate_ids, dc_grid) {
  flat <- as.integer(backend_res$cpd_flat)
  start <- as.integer(backend_res$cpd_start)
  len <- as.integer(backend_res$cpd_len)
  lapply(candidate_ids, function(id) {
    if (len[id] == 0L) {
      return(integer())
    }
    idx <- start[id] + seq_len(len[id])
    cpd <- flat[idx]
    if (!is.null(dc_grid)) {
      finite_idx <- is.finite(cpd)
      cpd[finite_idx] <- dc_grid[cpd[finite_idx]]
    }
    cpd
  })
}

.reliever_run_path_score <- function(backend_res, run_id) {
  start <- as.integer(backend_res$path_score_start[run_id])
  len <- as.integer(backend_res$path_score_len[run_id])
  if (len == 0L) {
    return(numeric())
  }
  as.numeric(backend_res$path_score_flat[start + seq_len(len)])
}

.reliever_wbs_stop_indices <- function(K, path_score, wbs_stop_crit) {
  gain <- as.numeric(path_score)
  max_k <- length(K) - 1L
  c(1L, vapply(wbs_stop_crit, function(crit) {
    stop_at <- which(gain <= crit)[1L]
    k <- if (is.na(stop_at)) max_k else stop_at - 1L
    max(0L, min(max_k, k)) + 1L
  }, integer(1L)))
}

.reliever_run_selector_map <- function(backend_res, run_id, setup, pen_val) {
  candidate_ids <- .reliever_run_candidate_indices(backend_res, run_id)
  n_candidate <- length(candidate_ids)
  if (setup$method %in% c("PELT", "OP")) {
    objective <- as.numeric(backend_res$objective[candidate_ids])
    if (length(objective) < n_candidate) {
      objective <- rep(NA_real_, n_candidate)
    }
    return(data.frame(
      run_id = run_id,
      select_value = c(Inf, pen_val)[seq_len(n_candidate)],
      candidate_id = seq_len(n_candidate),
      objective = objective[seq_len(n_candidate)]
    ))
  }
  if (setup$method %in% c("WBS", "WBS_recursive", "SeedBS", "BS") &&
      !is.null(setup$wbs_stop_crit)) {
    K <- as.integer(backend_res$cps_num[candidate_ids])
    return(data.frame(
      run_id = run_id,
      select_value = c(Inf, setup$wbs_stop_crit),
      candidate_id = .reliever_wbs_stop_indices(
        K, .reliever_run_path_score(backend_res, run_id),
        setup$wbs_stop_crit
      )
    ))
  }
  NULL
}

.reliever_deduplicate_penalty_candidates <- function(candidates,
                                                      selector_map) {
  runs <- split(candidates, candidates$run_id)
  for (i in seq_along(runs)) {
    run <- runs[[i]]
    cpd_key <- vapply(run$cpd, paste, character(1L), collapse = ",")
    key <- paste(run$K, cpd_key, sep = ":")
    old_to_new <- match(key, unique(key))
    map_id <- selector_map$run_id == run$run_id[1L]
    selector_map$candidate_id[map_id] <- old_to_new[
      selector_map$candidate_id[map_id]
    ]
    run <- run[!duplicated(key), , drop = FALSE]
    run$candidate_id <- seq_len(nrow(run))
    runs[[i]] <- run
  }
  candidates <- do.call(rbind, runs)
  rownames(candidates) <- NULL
  list(candidates = candidates, selector_map = selector_map)
}

.reliever_format_backend_cpd_result <- function(backend_res, setup, pen_val) {
  n_run <- length(backend_res$run_len)
  select_by <- if (setup$method %in% c("PELT", "OP")) {
    "pen_val"
  } else if (setup$method %in% c("WBS", "WBS_recursive", "SeedBS", "BS") &&
             !is.null(setup$wbs_stop_crit)) {
    "wbs_stop_crit"
  } else {
    "K"
  }
  candidates <- do.call(rbind, lapply(seq_len(n_run), function(i) {
    candidate_ids <- .reliever_run_candidate_indices(backend_res, i)
    data.frame(
      run_id = i,
      candidate_id = seq_along(candidate_ids),
      K = as.integer(backend_res$cps_num[candidate_ids]),
      loss = as.numeric(backend_res$loss[candidate_ids]),
      cpd = I(.reliever_decode_cpd(
        backend_res, candidate_ids, setup$dc_grid
      )),
      stringsAsFactors = FALSE
    )
  }))
  rownames(candidates) <- NULL
  selector_map <- lapply(seq_len(n_run), function(i) {
    .reliever_run_selector_map(backend_res, i, setup, pen_val)
  })
  selector_map <- Filter(Negate(is.null), selector_map)
  selector_map <- if (length(selector_map) == 0L) {
    NULL
  } else {
    out <- do.call(rbind, selector_map)
    rownames(out) <- NULL
    out
  }
  if (setup$method %in% c("PELT", "OP")) {
    deduplicated <- .reliever_deduplicate_penalty_candidates(
      candidates, selector_map
    )
    candidates <- deduplicated$candidates
    selector_map <- deduplicated$selector_map
  }
  candidates <- candidates[c("run_id", "K", "cpd", "loss", "candidate_id")]
  cpd_path <- list(
    select_by = select_by,
    candidates = candidates,
    selector_map = selector_map
  )
  path_score <- lapply(seq_len(n_run), function(i) {
    .reliever_run_path_score(backend_res, i)
  })
  if (!any(lengths(path_score) > 0L)) {
    path_score <- NULL
  }

  list(
    cpd_path = cpd_path,
    path_score = path_score,
    n_model_fit = as.numeric(backend_res$n_model_fit),
    model_fit_time = as.numeric(backend_res$model_fit_time),
    total_time = as.numeric(backend_res$total_time)
  )
}

.reliever_finalize_backend_result <- function(cpd_result, num_pruned,
                                               cache_profile) {
  list(
    cpd_path = cpd_result$cpd_path,
    n_model_fit = cpd_result$n_model_fit,
    model_fit_time = cpd_result$model_fit_time,
    cpd_time = cpd_result$total_time - cpd_result$model_fit_time,
    total_time = cpd_result$total_time,
    diagnostics = list(
      path_score = cpd_result$path_score,
      num_pruned = num_pruned
    ),
    cache_profile = cache_profile
  )
}

# Backend dispatch -------------------------------------------------------------

.reliever_run_by_loss_block_backend <- function(request) {
  setup <- request$setup
  run_cpd_ids <- as.integer(setup$run_cpd_ids)
  loss_block_cache <- request$backend_cache$loss_block_cache
  search_intervals <- .reliever_search_intervals(
    setup$method, setup$n, setup$dm, setup$M, request$wbs_seed
  )
  res_by_loss_block <- cpd_r_by_loss_block(
    setup$method, setup$n, setup$cpn_max, setup$dm, search_intervals,
    request$pen_val, setup$prune_value, run_cpd_ids,
    request$individual_loss_fun, setup$int_set$miss_cover_len,
    setup$int_set$int_len, setup$int_set$layer_point,
    setup$int_set$int_eps, loss_block_cache,
    request$return_cache_profile, request$owner_key
  )
  cpd_result <- .reliever_format_backend_cpd_result(
    backend_res = res_by_loss_block,
    setup = setup,
    pen_val = request$pen_val
  )
  cache_profile <- NULL
  if (request$return_cache_profile) {
    cache_profile <- .reliever_by_loss_block_cache_profile(
      res_by_loss_block$cache_profile,
      setup$int_set,
      request$loss_context
    )
  }

  num_pruned <- if (setup$method %in% c("PELT", "OP")) {
    res_by_loss_block$num_pruned
  } else {
    list()
  }

  .reliever_finalize_backend_result(
    cpd_result = cpd_result,
    num_pruned = num_pruned,
    cache_profile = cache_profile
  )
}

.reliever_run_by_cost_mat_backend <- function(request) {
  setup <- request$setup
  run_cpd_ids <- as.integer(setup$run_cpd_ids)
  cost_mat <- request$backend_cache$cost_mat

  if (is.null(cost_mat)) {
    cost_mat <- matrix(
      Inf, setup$n_loss_outputs, 1 + setup$n * (setup$n + 1) / 2
    )
  }

  if (setup$is_full) {
    miss_cover_len <- int_len <- layer_point <- integer()
    int_eps <- matrix(integer(), 0L, 2L)
  } else {
    miss_cover_len <- setup$int_set$miss_cover_len
    int_len <- setup$int_set$int_len
    layer_point <- setup$int_set$layer_point
    int_eps <- setup$int_set$int_eps
  }

  search_intervals <- .reliever_search_intervals(
    setup$method, setup$n, setup$dm, setup$M, request$wbs_seed
  )
  res_by_cost_mat <- cpd_r_by_cost_mat(
    setup$method, setup$n, setup$cpn_max, setup$dm, search_intervals,
    request$pen_val, setup$prune_value, run_cpd_ids, cost_mat,
    request$individual_loss_fun, miss_cover_len, int_len, layer_point, int_eps,
    setup$is_full
  )

  cpd_result <- .reliever_format_backend_cpd_result(
    backend_res = res_by_cost_mat,
    setup = setup,
    pen_val = request$pen_val
  )
  num_pruned <- if (setup$method %in% c("PELT", "OP")) {
    res_by_cost_mat$num_pruned
  } else {
    list()
  }

  cache_profile <- NULL
  if (request$return_cache_profile) {
    cache_profile <- .reliever_by_cost_mat_cache_profile(
      cost_mat,
      setup$int_set,
      request$loss_context
    )
  }

  .reliever_finalize_backend_result(
    cpd_result = cpd_result,
    num_pruned = num_pruned,
    cache_profile = cache_profile
  )
}
