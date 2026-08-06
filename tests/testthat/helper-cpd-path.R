cpd_path_candidates <- function(result) {
  if (is.null(result$cpd_path$candidates)) {
    stop("test result does not contain cpd_path$candidates.", call. = FALSE)
  }
  result$cpd_path$candidates
}

cpd_path_run <- function(result, run_id = 1L) {
  tab <- cpd_path_candidates(result)
  tab[tab$run_id == run_id, , drop = FALSE]
}

cpd_path_loss <- function(result, run_id = 1L) {
  cpd_path_run(result, run_id)$loss
}

cpd_path_K <- function(result, run_id = 1L) {
  cpd_path_run(result, run_id)$K
}

cpd_path_candidate <- function(result, run_id = 1L, K = NULL,
                               candidate_id = NULL) {
  tab <- cpd_path_run(result, run_id)
  if (!is.null(K)) {
    tab <- tab[tab$K == K, , drop = FALSE]
  }
  if (!is.null(candidate_id)) {
    tab <- tab[tab$candidate_id == candidate_id, , drop = FALSE]
  }
  if (nrow(tab) == 0L) {
    stop("no matching changepoint candidate.", call. = FALSE)
  }
  tab$cpd[[1L]]
}

cpd_path_matrix <- function(result, run_id = 1L) {
  tab <- cpd_path_run(result, run_id)
  tab <- tab[tab$K > 0L, , drop = FALSE]
  if (nrow(tab) == 0L) {
    return(matrix(NA_real_, 0L, 0L))
  }
  max_len <- max(lengths(tab$cpd))
  out <- matrix(NaN, nrow(tab), max_len)
  for (i in seq_len(nrow(tab))) {
    if (length(tab$cpd[[i]]) > 0L) {
      out[i, seq_along(tab$cpd[[i]])] <- tab$cpd[[i]]
    }
  }
  out
}

expect_cpd_error_lte <- function(estimate, truth, tolerance = 10L,
                                 info = NULL) {
  error <- relieverChangepoint::cp_error(
    as.integer(estimate),
    as.integer(truth)
  )
  label <- "changepoint-location error"
  if (!is.null(info)) {
    label <- paste0(label, " (", info, ")")
  }
  expect_lte(error, tolerance, label = label)
  invisible(error)
}

expect_same_cpd_path <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  x_tab <- cpd_path_candidates(x)
  y_tab <- cpd_path_candidates(y)
  expect_equal(x$cpd_path$select_by, y$cpd_path$select_by)
  expect_equal(x$cpd_path$selector_map, y$cpd_path$selector_map,
               tolerance = tolerance)
  key_cols <- c("run_id", "candidate_id", "K")
  expect_equal(x_tab[key_cols], y_tab[key_cols])
  expect_equal(x_tab$loss, y_tab$loss, tolerance = tolerance)
  expect_equal(x_tab$cpd, y_tab$cpd)
}

cpp_run_candidate_indices <- function(result, run_id = 1L) {
  start <- as.integer(result$run_start[run_id])
  len <- as.integer(result$run_len[run_id])
  if (len == 0L) {
    return(integer())
  }
  start + seq_len(len)
}

cpp_cpd_list <- function(result, run_id = 1L) {
  candidate_ids <- cpp_run_candidate_indices(result, run_id)
  flat <- as.integer(result$cpd_flat)
  start <- as.integer(result$cpd_start)
  len <- as.integer(result$cpd_len)
  lapply(candidate_ids, function(id) {
    if (len[id] == 0L) {
      return(integer())
    }
    flat[start[id] + seq_len(len[id])]
  })
}

cpp_cpd_matrix <- function(result, run_id = 1L) {
  cpd <- cpp_cpd_list(result, run_id)
  cpd <- cpd[lengths(cpd) > 0L]
  if (length(cpd) == 0L) {
    return(matrix(NA_real_, 0L, 0L))
  }
  max_len <- max(lengths(cpd))
  out <- matrix(NaN, length(cpd), max_len)
  for (i in seq_along(cpd)) {
    out[i, seq_along(cpd[[i]])] <- cpd[[i]]
  }
  out
}

cpp_loss <- function(result, run_id = 1L) {
  result$loss[cpp_run_candidate_indices(result, run_id)]
}

cpp_cps_num <- function(result, run_id = 1L) {
  result$cps_num[cpp_run_candidate_indices(result, run_id)]
}
