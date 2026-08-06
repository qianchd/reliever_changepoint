# Search algorithms -----------------------------------------------------------

test_that("reliever skeleton methods run with the null reg_fun", {
  with_test_timeout({
    data <- matrix(seq_len(24), ncol = 1)
    method_specs <- list(
      SN = list(),
      BS = list(M = 6L, wbs_seed = 1L),
      WBS = list(M = 6L, wbs_seed = 1L),
      WBS_recursive = list(M = 6L, wbs_seed = 1L),
      SeedBS = list(M = 6L),
      PELT = list(pen_val = c(0.1, 1)),
      OP = list(pen_val = c(0.1, 1))
    )

    for (method in names(method_specs)) {
      args <- c(
        list(
          data = data,
          cpn_max = 2L,
          dm = 5L,
          cov_rate = 0.8,
          reg_fun = reg_null,
          method = method,
          echo = FALSE
        ),
        method_specs[[method]]
      )
      res <- do.call(relieverChangepoint::reliever_generic, args)
      tab <- cpd_path_candidates(res)

      expected_rows <- if (method %in% c("PELT", "OP")) 1L else 3L
      expect_equal(nrow(tab), expected_rows)
      expect_equal(nrow(res$run_meta), 1L)
      expect_equal(res$run_meta$loss_output_id, 1L)
      expect_true(all(is.finite(tab$loss)))
      expect_true(all(is.finite(tab$K)))
      expect_false("cost_mat" %in% names(res))
      expect_false("int_set" %in% names(res))
    }
  })
})

test_that("SN supports zero requested changepoints", {
  with_test_timeout({
    data <- matrix(seq_len(20), ncol = 1L)

    res <- relieverChangepoint::reliever_generic(
      data,
      cpn_max = 0L,
      dm = 5L,
      cov_rate = 1,
      reg_fun = reg_null,
      method = "SN",
      cpn_crit = "sic",
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    expect_equal(nrow(cpd_path_candidates(res)), 1L)
    expect_equal(cpd_path_K(res), 0L)
    expect_equal(res$summary$K_hat, 0L)
    expect_equal(res$summary$cpd_hat[[1L]], integer())

    loss_block_res <- relieverChangepoint::reliever_generic(
      data,
      cpn_max = 0L,
      dm = 5L,
      cov_rate = 0.8,
      reg_fun = reg_null,
      method = "SN",
      cpn_crit = "sic",
      cache_backend = "by_loss_block",
      echo = FALSE
    )
    expect_equal(nrow(cpd_path_candidates(loss_block_res)), 1L)
    expect_equal(cpd_path_K(loss_block_res), 0L)
    expect_equal(loss_block_res$summary$K_hat, 0L)
    expect_equal(loss_block_res$summary$cpd_hat[[1L]], integer())

    mean_res <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 0L,
      dm = 5L,
      cov_rate = 1,
      method = "SN",
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    expect_equal(nrow(cpd_path_candidates(mean_res)), 1L)
    expect_equal(cpd_path_K(mean_res), 0L)
    expect_equal(nrow(mean_res$summary), 0L)

    mean_loss_block_res <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 0L,
      dm = 5L,
      cov_rate = 0.8,
      method = "SN",
      cache_backend = "by_loss_block",
      echo = FALSE
    )
    expect_equal(nrow(cpd_path_candidates(mean_loss_block_res)), 1L)
    expect_equal(cpd_path_K(mean_loss_block_res), 0L)
    expect_equal(nrow(mean_loss_block_res$summary), 0L)
  })
})

test_that("WBS stops cleanly when no finite split gain is available", {
  with_test_timeout({
    n <- 90L
    cost_mat <- matrix(Inf, nrow = 1L, ncol = 1L + n * (n + 1L) / 2L)
    inf_loss_fun <- function(l, r, l_end, r_end) {
      matrix(Inf, nrow = r_end - l_end + 1L, ncol = 1L)
    }
    lr_m <- relieverChangepoint::create_wbs_itv(n, dm = 30L, M = 5L, wbs_seed = 1L)

    res <- relieverChangepoint:::cpd_r_by_cost_mat(
      "WBS", n, 2L, 30L, lr_m, numeric(), 0, as.integer(1L),
      cost_mat, inf_loss_fun,
      integer(), integer(), integer(), matrix(integer(), 0L, 2L),
      TRUE
    )

    expect_equal(cpp_cpd_matrix(res), matrix(NA_real_, 0L, 0L))
    expect_equal(cpp_cps_num(res), 0L)
  })
})

# Cache profiles and reuse ----------------------------------------------------

test_that("detail returns reusable cache objects for by_cost_mat and by_loss_block backends", {
  with_test_timeout({
    data <- matrix(seq_len(24), ncol = 1)
    common_args <- list(
      data = data,
      cpn_max = 1L,
      dm = 5L,
      cov_rate = 0.8,
      reg_fun = reg_null,
      method = "SN",
      echo = FALSE
    )

    by_cost_mat <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(detail = TRUE, cache_backend = "by_cost_mat"))
    )
    expect_false(is.null(by_cost_mat$cache_profile))
    expect_equal(by_cost_mat$cache_profile$backend, "by_cost_mat")
    expect_false("cost_mat" %in% names(by_cost_mat))
    expect_false("int_set" %in% names(by_cost_mat))
    expect_true(is.matrix(by_cost_mat$cache_profile$objects$cost_mat))
    expect_false(is.null(by_cost_mat$cache_profile$objects$int_set))

    by_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(detail = TRUE, cache_backend = "by_loss_block"))
    )
    expect_false(is.null(by_loss_block$cache_profile))
    expect_equal(by_loss_block$cache_profile$backend, "by_loss_block")
    expect_false("cost_mat" %in% names(by_loss_block))
    expect_false("int_set" %in% names(by_loss_block))
    expect_null(by_loss_block$cache_profile$objects$cost_mat)
    expect_false(is.null(by_loss_block$cache_profile$objects$int_set))
    expect_false(is.null(by_loss_block$cache_profile$objects$loss_block_cache))
    expect_true(length(by_loss_block$cache_profile$objects$loss_block_cache$blocks) > 0L)
    expect_true("owner_key" %in% names(by_loss_block$cache_profile$objects$loss_block_cache))
    expect_type(by_loss_block$cache_profile$objects$loss_block_cache$owner_key,
                "integer")

    compact_by_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_loss_block"))
    )
    expect_null(compact_by_loss_block$cache_profile)
    expect_false("cost_mat" %in% names(compact_by_loss_block))
    expect_false("int_set" %in% names(compact_by_loss_block))

    no_owner_key <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(
        detail = TRUE,
        cache_backend = "by_loss_block",
        owner_key = FALSE
      ))
    )
    expect_same_cpd_path(no_owner_key, by_loss_block)
    expect_false(
      "owner_key" %in%
        names(no_owner_key$cache_profile$objects$loss_block_cache)
    )

    cost_mat_no_owner_key <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(detail = TRUE, owner_key = FALSE))
    )
    expect_same_cpd_path(cost_mat_no_owner_key, by_cost_mat)
    expect_false("owner_key" %in% names(cost_mat_no_owner_key$cache_profile$objects))

  })
})

test_that("cov_rate at the full-search threshold is normalized to one", {
  with_test_timeout({
    data <- matrix(seq_len(12), ncol = 1)

    res <- relieverChangepoint::reliever_generic(
      data = data,
      cpn_max = 1L,
      dm = 4L,
      cov_rate = 0.95,
      reg_fun = reg_null,
      method = "SN",
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    expect_equal(res$settings$cov_rate, 1)
    expect_true(all(is.finite(cpd_path_loss(res))))

    mean_res <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 1L,
      dm = 4L,
      cov_rate = 0.95,
      method = "SN",
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    expect_equal(mean_res$settings$cov_rate, 1)
    expect_true(all(is.finite(cpd_path_loss(mean_res))))

    switched <- NULL
    expect_warning(
      switched <- relieverChangepoint::reliever_generic(
        data = data,
        cpn_max = 1L,
        dm = 4L,
        cov_rate = 0.95,
        reg_fun = reg_null,
        method = "SN",
        cache_backend = "by_loss_block",
        echo = FALSE
      ),
      "using cache_backend = \"by_cost_mat\""
    )
    expect_equal(switched$settings$cache_backend, "by_cost_mat")
    expect_equal(switched$settings$cov_rate, 1)
    expect_true(all(is.finite(cpd_path_loss(switched))))
  })
})

test_that("by_cost_mat reliever can reuse a cache profile", {
  with_test_timeout({
    signal_reg <- function(data, l, r, l_end = l, r_end = r,
                           save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      loss <- matrix((y - center)^2, ncol = 1)
      list(loss = loss, model = NULL)
    }

    data <- matrix(c(rep(0, 10), rep(5, 10)), ncol = 1)
    args <- list(
      data = data,
      cpn_max = 1L,
      dm = 4L,
      cov_rate = 0.8,
      reg_fun = signal_reg,
      method = "SN",
      detail = TRUE,
      cache_backend = "by_cost_mat",
      echo = FALSE
    )

    fresh <- do.call(relieverChangepoint::reliever_generic, args)
    reused <- do.call(
      relieverChangepoint::reliever_generic,
      c(args, list(cache_profile = fresh$cache_profile))
    )

    expect_same_cpd_path(reused, fresh)
    expect_equal(sum(reused$timing$n_model_fit), 0)
  })
})

test_that("cache_profile reuse rejects matrices that cannot zero-copy into arma::mat", {
  with_test_timeout({
    data <- matrix(seq_len(24), ncol = 1)
    common_args <- list(
      data = data,
      cpn_max = 1L,
      dm = 5L,
      cov_rate = 0.8,
      reg_fun = reg_null,
      method = "SN",
      detail = TRUE,
      echo = FALSE
    )

    by_cost_mat <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_cost_mat"))
    )
    cost_dim <- dim(by_cost_mat$cache_profile$objects$cost_mat)
    bad_cost_profiles <- list(
      logical_na = matrix(NA, cost_dim[1L], cost_dim[2L]),
      integer = matrix(1L, cost_dim[1L], cost_dim[2L]),
      logical = matrix(TRUE, cost_dim[1L], cost_dim[2L]),
      null = NULL
    )
    for (bad_cost_mat in bad_cost_profiles) {
      bad_profile <- by_cost_mat$cache_profile
      bad_profile$objects$cost_mat <- bad_cost_mat
      expect_error(
        do.call(
          relieverChangepoint::reliever_generic,
          c(common_args, list(
            cache_backend = "by_cost_mat",
            cache_profile = bad_profile
          ))
        ),
        "double matrix"
      )
    }

    by_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_loss_block"))
    )
    loss_dim <- dim(
      by_loss_block$cache_profile$objects$loss_block_cache$blocks[[1L]]$loss
    )
    bad_loss_blocks <- list(
      logical_na = matrix(NA, loss_dim[1L], loss_dim[2L]),
      integer = matrix(1L, loss_dim[1L], loss_dim[2L]),
      logical = matrix(TRUE, loss_dim[1L], loss_dim[2L]),
      null = NULL
    )
    for (bad_loss in bad_loss_blocks) {
      bad_profile <- by_loss_block$cache_profile
      bad_profile$objects$loss_block_cache$blocks[[1L]]$loss <- bad_loss
      expect_error(
        do.call(
          relieverChangepoint::reliever_generic,
          c(common_args, list(
            cache_backend = "by_loss_block",
            cache_profile = bad_profile
          ))
        ),
        "double matrix"
      )
    }
  })
})

test_that("by_cost_mat reuses fitted losses across loss outputs in one run", {
  with_test_timeout({
    two_row_reg <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(2L)
      }
      y <- seq_len(r_end - l_end + 1L)
      loss <- cbind(abs(y), y^2)
      list(loss = loss, model = NULL)
    }

    data <- matrix(seq_len(45), ncol = 1)
    res <- relieverChangepoint::reliever_generic(
      data = data,
      cpn_max = 2L,
      dm = 8L,
      cov_rate = 0.8,
      reg_fun = two_row_reg,
      method = "SN",
      run_cpd_ids = 1:2,
      cache_backend = "by_cost_mat",
      detail = TRUE,
      echo = FALSE
    )

    expect_equal(length(res$timing$n_model_fit), 2L)
    expect_gt(res$timing$n_model_fit[1L], 0)
    expect_equal(res$timing$n_model_fit[2L], 0)
    expect_gt(sum(is.finite(res$cache_profile$objects$cost_mat[1L, ])), 0)
    expect_equal(
      sum(is.finite(res$cache_profile$objects$cost_mat[2L, ])),
      sum(is.finite(res$cache_profile$objects$cost_mat[1L, ]))
    )
  })
})

test_that("by_cost_mat cache profile can be partially reused and updated", {
  with_test_timeout({
    signal_reg <- function(data, l, r, l_end = l, r_end = r,
                           save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      loss <- matrix((y - center)^2, ncol = 1)
      list(loss = loss, model = NULL)
    }

    data <- matrix(c(rep(0, 12), rep(4, 12)), ncol = 1)
    common <- list(
      data = data,
      cpn_max = 1L,
      dm = 4L,
      cov_rate = 0.8,
      reg_fun = signal_reg,
      detail = TRUE,
      echo = FALSE
    )

    sn <- do.call(
      relieverChangepoint::reliever_generic,
      c(common, list(method = "SN", cache_backend = "by_cost_mat"))
    )
    filled_before <- sum(is.finite(sn$cache_profile$objects$cost_mat))
    reused_wbs <- do.call(
      relieverChangepoint::reliever_generic,
      c(
        common,
        list(
          method = "WBS",
          M = 8L,
          wbs_seed = 123L,
          cache_backend = "by_cost_mat",
          cache_profile = sn$cache_profile
        )
      )
    )
    fresh_wbs <- do.call(
      relieverChangepoint::reliever_generic,
      c(
        common,
        list(
          method = "WBS",
          M = 8L,
          wbs_seed = 123L,
          cache_backend = "by_cost_mat"
        )
      )
    )

    expect_same_cpd_path(reused_wbs, fresh_wbs)
    expect_lte(
      sum(reused_wbs$timing$n_model_fit),
      sum(fresh_wbs$timing$n_model_fit)
    )
    expect_gte(
      sum(is.finite(reused_wbs$cache_profile$objects$cost_mat)),
      filled_before
    )

    by_loss_block_profile <- do.call(
      relieverChangepoint::reliever_generic,
      c(common, list(method = "SN", cache_backend = "by_loss_block"))
    )
    expect_error(
      do.call(
        relieverChangepoint::reliever_generic,
        c(
          common,
          list(
            cache_backend = "by_cost_mat",
            cache_profile = by_loss_block_profile$cache_profile
          )
        )
      ),
      "by_cost_mat"
    )
  })
})

test_that("by_loss_block cache profile can be reused", {
  with_test_timeout({
    signal_reg <- function(data, l, r, l_end = l, r_end = r,
                           save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      loss <- matrix((y - center)^2, ncol = 1)
      list(loss = loss, model = NULL)
    }

    data <- matrix(c(rep(0, 12), rep(4, 12)), ncol = 1)
    args <- list(
      data = data,
      cpn_max = 1L,
      dm = 4L,
      cov_rate = 0.8,
      reg_fun = signal_reg,
      method = "SN",
      cache_backend = "by_loss_block",
      detail = TRUE,
      echo = FALSE
    )

    fresh <- do.call(relieverChangepoint::reliever_generic, args)
    reused <- do.call(
      relieverChangepoint::reliever_generic,
      c(args, list(cache_profile = fresh$cache_profile))
    )

    expect_same_cpd_path(reused, fresh)
    expect_equal(sum(reused$timing$n_model_fit), 0)
    expect_false(is.null(reused$cache_profile$objects$loss_block_cache))
    expect_true("owner_key" %in% names(reused$cache_profile$objects$loss_block_cache))
    expect_type(reused$cache_profile$objects$loss_block_cache$owner_key,
                "integer")
    expect_equal(
      length(reused$cache_profile$objects$loss_block_cache$blocks),
      length(fresh$cache_profile$objects$loss_block_cache$blocks)
    )

    mismatch_args <- args
    mismatch_args$cov_rate <- 0.7
    expect_error(
      do.call(
        relieverChangepoint::reliever_generic,
        c(mismatch_args, list(cache_profile = fresh$cache_profile))
      ),
      "same Reliever interval set"
    )
  })
})

test_that("cache profiles reject a different loss context", {
  with_test_timeout({
    scaled_mean <- function(data, l, r, l_end = l, r_end = r,
                            loss_scale = 1, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      center <- mean(data[l:r, 1L])
      list(loss = matrix(
        loss_scale * (data[l_end:r_end, 1L] - center)^2,
        ncol = 1L
      ))
    }

    data <- matrix(c(rep(0, 12), rep(3, 12)), ncol = 1L)
    for (backend in c("by_cost_mat", "by_loss_block")) {
      fit <- relieverChangepoint::reliever_generic(
        data = data, reg_fun = scaled_mean,
        cpn_max = 1L, dm = 4L, cov_rate = 0.8,
        detail = TRUE, cache_backend = backend,
        loss_scale = 1
      )

      expect_error(
        relieverChangepoint::reliever_generic(
          data = data, reg_fun = scaled_mean,
          cpn_max = 1L, dm = 4L, cov_rate = 0.8,
          detail = TRUE, cache_backend = backend,
          cache_profile = fit$cache_profile,
          loss_scale = 2
        ),
        "different data, reg_fun, dc_grid, or loss-function arguments",
        info = backend
      )
      expect_error(
        relieverChangepoint::reliever_generic(
          data = data + 1, reg_fun = scaled_mean,
          cpn_max = 1L, dm = 4L, cov_rate = 0.8,
          detail = TRUE, cache_backend = backend,
          cache_profile = fit$cache_profile,
          loss_scale = 1
        ),
        "different data, reg_fun, dc_grid, or loss-function arguments",
        info = backend
      )
    }
  })
})

# Search edge cases -----------------------------------------------------------

test_that("all search methods return the no-split path when n is below 2 * dm", {
  with_test_timeout({
    set.seed(2026)
    data <- matrix(stats::rnorm(20), 10L, 2L)
    methods <- c("SN", "WBS", "WBS_recursive", "SeedBS", "BS", "PELT", "OP")

    for (method in methods) {
      fit <- relieverChangepoint::reliever_mean(
        data, cpn_max = 3L, dm = 8L, cov_rate = 0.8, method = method,
        pen_val = c(0.5, 2), M = 8L, wbs_seed = 1L
      )
      expect_true(all(fit$cpd_path$candidates$K == 0L), info = method)
      expect_true(all(lengths(fit$cpd_path$candidates$cpd) == 0L), info = method)
      expect_true(all(is.finite(fit$cpd_path$candidates$loss)), info = method)
    }

    for (method in c("PELT", "OP")) {
      fit <- relieverChangepoint::reliever_mean(
        data, cpn_max = 3L, dm = 8L, cov_rate = 0.8, method = method,
        pen_val = c(0.5, 2), cache_backend = "by_cost_mat"
      )
      expect_equal(fit$cpd_path$candidates$K, 0L, info = method)
      expect_true(all(is.finite(fit$cpd_path$candidates$loss)), info = method)
    }
  })
})

test_that("PELT and OP validate penalty controls before entering C++", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1L)
    invalid_penalties <- list(NULL, NA_real_, Inf, -1, matrix(1, 1L, 1L))
    for (method in c("PELT", "OP")) {
      for (penalty in invalid_penalties) {
        expect_error(
          relieverChangepoint::reliever_mean(
            data, dm = 5L, method = method, pen_val = penalty
          ),
          "pen_val must be a non-empty vector of finite, non-negative numbers",
          info = method
        )
      }
    }
    expect_error(
      relieverChangepoint::reliever_mean(
        data, dm = 5L, method = "PELT", pen_val = 1, prune_value = NA_real_
      ),
      "prune_value must be a finite numeric scalar"
    )
  })
})

test_that("SeedBS changepoint paths do not depend on the global RNG state", {
  with_test_timeout({
    set.seed(2026)
    data <- matrix(c(
      stats::rnorm(40), stats::rnorm(40, 3), stats::rnorm(40, -2)
    ), ncol = 1L)

    set.seed(1)
    seed_one <- relieverChangepoint::reliever_mean(
      data, cpn_max = 3L, dm = 10L, cov_rate = 0.8, method = "SeedBS"
    )
    set.seed(2)
    seed_two <- relieverChangepoint::reliever_mean(
      data, cpn_max = 3L, dm = 10L, cov_rate = 0.8, method = "SeedBS"
    )
    expect_same_cpd_path(seed_one, seed_two)
  })
})

test_that("native and generic mean PELT expose initialized pruning diagnostics", {
  with_test_timeout({
    set.seed(2026)
    data <- matrix(stats::rnorm(60), 30L, 2L)
    args <- list(
      data = data, dm = 5L, cov_rate = 1, method = "PELT",
      pen_val = c(0.5, 2), prune_value = 0, detail = TRUE,
      cache_backend = "by_cost_mat"
    )
    native <- do.call(relieverChangepoint::reliever_mean, args)
    generic <- do.call(
      relieverChangepoint::reliever_generic,
      c(args, list(reg_fun = relieverChangepoint::reg_fun_mean))
    )

    expect_equal(native$diagnostics$num_pruned,
                 generic$diagnostics$num_pruned)
    expect_length(native$diagnostics$num_pruned, 1L)
    expect_true(is.matrix(native$diagnostics$num_pruned[[1L]]))
    expect_equal(native$diagnostics$num_pruned[[1L]][, seq_len(10L)],
                 matrix(0, 2L, 10L))
  })
})

test_that("PELT delayed pruning and OP match the exact full mean path", {
  with_test_timeout({
    set.seed(1)
    data <- stats::rnorm(30L)
    dm <- 3L
    pen_val <- c(0.01, 0.1, 1, 5)
    full_path <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = floor(length(data) / dm) - 1L,
      dm = dm,
      cov_rate = 1,
      method = "SN",
      cache_backend = "by_cost_mat"
    )
    exact_objective <- vapply(
      pen_val,
      function(penalty) {
        min(
          full_path$cpd_path$candidates$loss +
            penalty * full_path$cpd_path$candidates$K
        )
      },
      numeric(1L)
    )
    pelt <- relieverChangepoint::reliever_mean(
      data, dm = dm, cov_rate = 1, method = "PELT",
      pen_val = pen_val, cache_backend = "by_cost_mat"
    )
    op <- relieverChangepoint::reliever_mean(
      data, dm = dm, cov_rate = 1, method = "OP",
      pen_val = pen_val, cache_backend = "by_cost_mat"
    )

    expect_equal(
      pelt$cpd_path$selector_map$objective[-1L],
      exact_objective,
      tolerance = 1e-10
    )
    expect_equal(
      op$cpd_path$selector_map$objective[-1L],
      exact_objective,
      tolerance = 1e-10
    )

    expect_silent(
      relieverChangepoint::reliever_mean(
        data, dm = 1L, cov_rate = 1, method = "OP",
        pen_val = 0.1, cache_backend = "by_cost_mat"
      )
    )

    set.seed(1)
    scaled_data <- matrix(stats::rnorm(24L) * 1e6, ncol = 1L)
    scaled_path <- relieverChangepoint::reliever(
      scaled_data,
      cpd_family = "mean_crossfit",
      cpn_max = 7L,
      dm = 3L,
      cov_rate = 1,
      method = "SN",
      nfolds = 2L,
      cache_backend = "by_cost_mat"
    )
    scaled_penalty <- 1e9
    scaled_exact <- vapply(
      split(
        scaled_path$cpd_path$candidates,
        scaled_path$cpd_path$candidates$run_id
      ),
      function(path) {
        min(path$loss + scaled_penalty * path$K)
      },
      numeric(1L)
    )
    scaled_op <- relieverChangepoint::reliever(
      scaled_data,
      cpd_family = "mean_crossfit",
      dm = 3L,
      cov_rate = 1,
      method = "OP",
      pen_val = scaled_penalty,
      nfolds = 2L,
      cache_backend = "by_cost_mat"
    )
    expect_equal(
      scaled_op$cpd_path$selector_map$objective[
        is.finite(scaled_op$cpd_path$selector_map$select_value)
      ],
      unname(scaled_exact),
      tolerance = 1e-8
    )
  })
})

test_that("PELT runs silently across loss and search configurations", {
  with_test_timeout({
    set.seed(20260724)
    data <- matrix(stats::rnorm(30L * 2L), nrow = 30L)

    expect_silent(
      relieverChangepoint::reliever(
        data,
        cpd_family = "mean_crossfit",
        nfolds = 2L,
        dm = 5L,
        cov_rate = 1,
        method = "PELT",
        pen_val = 0.5,
        cache_backend = "by_cost_mat"
      )
    )
    expect_silent(
      relieverChangepoint::reliever_mean(
        data,
        dm = 5L,
        cov_rate = 0.8,
        method = "PELT",
        pen_val = 0.5
      )
    )
  })
})

# Cache diagnostics -----------------------------------------------------------

test_that("detail cache profiles can be reused across owner-key modes", {
  with_test_timeout({
    signal_reg <- function(data, l, r, l_end = l, r_end = r,
                           save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      loss <- matrix((y - center)^2, ncol = 1)
      list(loss = loss, model = NULL)
    }

    data <- matrix(c(rep(0, 14), rep(3, 14)), ncol = 1)
    base_args <- list(
      data = data,
      cpn_max = 1L,
      dm = 4L,
      cov_rate = 0.8,
      reg_fun = signal_reg,
      method = "SN",
      detail = TRUE,
      echo = FALSE
    )
    specs <- list(
      by_cost_mat = list(cache_backend = "by_cost_mat", owner_key = FALSE),
      by_loss_block_owner_false =
        list(cache_backend = "by_loss_block", owner_key = FALSE),
      by_loss_block_owner_true =
        list(cache_backend = "by_loss_block", owner_key = TRUE)
    )

    for (spec_name in names(specs)) {
      args <- c(base_args, specs[[spec_name]])
      fresh <- do.call(relieverChangepoint::reliever_generic, args)
      reused <- do.call(
        relieverChangepoint::reliever_generic,
        c(args, list(cache_profile = fresh$cache_profile))
      )

      expect_same_cpd_path(reused, fresh)
    expect_equal(sum(reused$timing$n_model_fit), 0,
                   info = spec_name)
      expect_false(is.null(reused$cache_profile),
                   info = spec_name)

      if (identical(specs[[spec_name]]$cache_backend, "by_cost_mat")) {
        expect_true(is.matrix(reused$cache_profile$objects$cost_mat),
                    info = spec_name)
        expect_false("owner_key" %in% names(reused$cache_profile$objects),
                     info = spec_name)
      } else {
        cache_names <- names(reused$cache_profile$objects$loss_block_cache)
        expect_equal(
          "owner_key" %in% cache_names,
          isTRUE(specs[[spec_name]]$owner_key),
          info = spec_name
        )
      }
    }
  })
})

test_that("by_cost_mat skeleton methods return reusable cache detail", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1)
    method_specs <- list(
      BS = list(M = 0L),
      WBS = list(M = 8L, wbs_seed = 10L),
      SeedBS = list(M = 8L),
      PELT = list(pen_val = c(0.1, 1)),
      OP = list(pen_val = c(0.1, 1))
    )

    for (method in names(method_specs)) {
      args <- c(
        list(
          data = data,
          cpn_max = 2L,
          dm = 5L,
          cov_rate = 0.85,
          reg_fun = reg_null,
          method = method,
          detail = TRUE,
          cache_backend = "by_cost_mat",
          echo = FALSE
        ),
        method_specs[[method]]
      )
      res <- do.call(relieverChangepoint::reliever_generic, args)

    expect_equal(res$settings$cache_backend, "by_cost_mat")
      expect_false(is.null(res$cache_profile))
      expect_false("cost_mat" %in% names(res))
      expect_false("int_set" %in% names(res))
      expect_true(is.matrix(res$cache_profile$objects$cost_mat))
      expect_false(is.null(res$cache_profile$objects$int_set))
      expect_true(all(is.finite(cpd_path_loss(res))))
    expect_true(all(is.finite(cpd_path_K(res))))
    }
  })
})

# Candidate paths and selectors ----------------------------------------------

test_that("reliever cpd_path records K, penalty, and WBS stopping selectors", {
  with_test_timeout({
    data <- matrix(c(
      rep(0, 30L),
      rep(3, 30L),
      rep(-2, 30L)
    ), ncol = 1L)

    sn <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 2L,
      dm = 15L,
      cov_rate = 0.6,
      method = "SN",
      cpn_crit = "loss",
      echo = FALSE
    )
    expect_equal(sn$cpd_path$select_by, "K")
    expect_true(all(c(
      "run_id", "candidate_id", "K", "loss", "cpd"
    ) %in% names(sn$cpd_path$candidates)))
    expect_null(sn$cpd_path$selector_map)
    expect_null(sn$diagnostics)
    expect_equal(sn$cpd_path$candidates$K, 0:2)
    expect_equal(
      lengths(sn$cpd_path$candidates$cpd), sn$cpd_path$candidates$K
    )

    pelt <- relieverChangepoint::reliever_generic(
      data = data,
      cpn_max = 2L,
      dm = 15L,
      cov_rate = 0.6,
      reg_fun = reg_null,
      method = "PELT",
      pen_val = c(0.5, 2),
      cpn_crit = "none",
      echo = FALSE
    )
    expect_equal(pelt$cpd_path$select_by, "pen_val")
    expect_equal(pelt$cpd_path$selector_map$select_value, c(Inf, 0.5, 2))
    expect_true(all(c("candidate_id", "objective") %in%
                      names(pelt$cpd_path$selector_map)))
    expect_equal(pelt$summary$pen_val, c(0.5, 2))
    expect_equal(
      pelt$summary$K_hat,
      pelt$cpd_path$candidates$K[
        pelt$cpd_path$selector_map$candidate_id[-1L]
      ]
    )
    expect_equal(
      length(unique(pelt$cpd_path$selector_map$candidate_id)),
      nrow(pelt$cpd_path$candidates)
    )
    expect_equal(pelt$settings$cpn_max, 2L)
    expect_equal(pelt$settings$pen_val, c(0.5, 2))
    expect_equal(pelt$settings$prune_value, 0)

    wbs_base <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 2L,
      dm = 15L,
      cov_rate = 0.6,
      method = "WBS",
      M = 30L,
      wbs_seed = 17L,
      cpn_crit = "loss",
      detail = TRUE,
      echo = FALSE
    )
    gains <- wbs_base$diagnostics$path_score[[1L]]
    expect_equal(length(gains), 2L)
    expect_true(all(is.finite(gains)))

    stop_crit <- c(max(gains) + 1, min(gains) - 1)
    wbs_stop <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 2L,
      dm = 15L,
      cov_rate = 0.6,
      method = "WBS",
      M = 30L,
      wbs_seed = 17L,
      wbs_stop_crit = stop_crit,
      cpn_crit = "none",
      echo = FALSE
    )
    expect_equal(wbs_stop$cpd_path$select_by, "wbs_stop_crit")
    expect_equal(wbs_stop$cpd_path$selector_map$select_value, c(Inf, stop_crit))
    expect_equal(wbs_stop$cpd_path$selector_map$candidate_id, c(1L, 1L, 3L))
    expect_equal(wbs_stop$cpd_path$candidates$K, 0:2)
    expect_equal(nrow(wbs_stop$cpd_path$candidates), 3L)
    expect_equal(wbs_stop$summary$wbs_stop_crit, stop_crit)
    expect_equal(wbs_stop$summary$K_hat, c(0L, 2L))
    expect_equal(wbs_stop$settings$cpn_max, 2L)
    expect_equal(wbs_stop$settings$M, 30L)
    expect_equal(wbs_stop$settings$wbs_seed, 17L)
    expect_equal(wbs_stop$settings$wbs_stop_crit, stop_crit)

    wbs_recursive_base <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 2L,
      dm = 15L,
      cov_rate = 0.6,
      method = "WBS_recursive",
      M = 30L,
      wbs_seed = 17L,
      cpn_crit = "loss",
      detail = TRUE,
      echo = FALSE
    )
    recursive_gains <- wbs_recursive_base$diagnostics$path_score[[1L]]
    expect_equal(length(recursive_gains), 2L)
    expect_true(all(is.finite(recursive_gains)))

    recursive_stop_crit <- c(max(recursive_gains) + 1, min(recursive_gains) - 1)
    wbs_recursive_stop <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 2L,
      dm = 15L,
      cov_rate = 0.6,
      method = "WBS_recursive",
      M = 30L,
      wbs_seed = 17L,
      wbs_stop_crit = recursive_stop_crit,
      cpn_crit = "none",
      echo = FALSE
    )
    expect_equal(wbs_recursive_stop$cpd_path$select_by, "wbs_stop_crit")
    expect_equal(
      wbs_recursive_stop$cpd_path$selector_map$select_value,
      c(Inf, recursive_stop_crit)
    )
    expect_equal(
      wbs_recursive_stop$cpd_path$selector_map$candidate_id,
      c(1L, 1L, 3L)
    )
    expect_equal(nrow(wbs_recursive_stop$cpd_path$candidates), 3L)
  })
})

test_that("wbs_stop_crit is restricted to WBS-family methods", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1L)
    expect_error(
      relieverChangepoint::reliever_mean(
        data,
        cpn_max = 1L,
        dm = 10L,
        cov_rate = 0.6,
        method = "SN",
        wbs_stop_crit = 1,
        echo = FALSE
      ),
      "wbs_stop_crit can only be used"
    )
  })
})

# Compressed search grids -----------------------------------------------------

test_that("dc_grid_size is the regular-grid fitting interface", {
  with_test_timeout({
    data <- matrix(seq_len(20), ncol = 1L)
    expected_grid <- c(6L, 12L, 18L, 20L)
    common_args <- list(
      X = data, cpn_max = 1L, dm = 2L, cov_rate = 0.6, method = "SN",
      cache_backend = "by_cost_mat"
    )

    regular <- do.call(
      relieverChangepoint::reliever,
      c(common_args, list(dc_grid_size = 6L))
    )
    explicit <- do.call(
      relieverChangepoint::reliever,
      c(common_args, list(dc_grid = expected_grid))
    )

    expect_equal(regular$settings$dc_grid_size, 6L)
    expect_equal(regular$settings$dc_grid, expected_grid)
    expect_same_cpd_path(regular, explicit)
    expect_error(
      do.call(
        relieverChangepoint::reliever,
        c(common_args, list(
          dc_grid_size = 6L,
          dc_grid = expected_grid
        ))
      ),
      "Supply only one"
    )
  })
})

test_that("dc_grid uses compressed grid internally and maps candidates back", {
  with_test_timeout({
    data <- matrix(seq_len(20), ncol = 1)
    dc_grid <- seq(2L, 20L, by = 2L)

    res <- relieverChangepoint::reliever_generic(
      data = data,
      reg_fun = reg_null,
      cpn_max = 1L,
      dm = 2L,
      cov_rate = 0.6,
      method = "SN",
      detail = TRUE,
      cache_backend = "by_cost_mat",
      echo = FALSE,
      dc_grid = dc_grid
    )

    expect_equal(
      ncol(res$cache_profile$objects$cost_mat),
      1L + length(dc_grid) * (length(dc_grid) + 1L) / 2L
    )
    expect_equal(res$cache_profile$objects$int_set$n, length(dc_grid))
    expect_true(all(unlist(res$cpd_path$candidates$cpd) %in% dc_grid))
    expect_equal(nrow(res$cpd_path$candidates), 2L)
    expect_equal(res$settings$dm, 2L)
    expect_equal(res$settings$search_dm, 1L)
    expect_equal(res$settings$dc_grid, dc_grid)
    expect_false(is.null(res$cache_profile))
  })
})

test_that("dc_grid keeps dm on the original sample scale", {
  with_test_timeout({
    data <- matrix(seq_len(90), ncol = 1)
    dc_grid <- seq(10L, 90L, by = 10L)

    res <- relieverChangepoint::reliever_generic(
      data = data,
      reg_fun = reg_null,
      cpn_max = 1L,
      dm = 30L,
      cov_rate = 0.6,
      method = "SN",
      detail = TRUE,
      cache_backend = "by_cost_mat",
      echo = FALSE,
      dc_grid = dc_grid
    )

    expect_equal(res$settings$n, 90L)
    expect_equal(res$settings$dm, 30L)
    expect_equal(res$settings$search_n, 9L)
    expect_equal(res$settings$search_dm, 3L)
    expect_equal(res$settings$dc_grid, dc_grid)
    expect_equal(res$cache_profile$objects$int_set$n, 9L)
    expect_equal(nrow(res$cpd_path$candidates), 2L)
  })
})

test_that("by_loss_block backend supports dc_grid candidate mapping", {
  with_test_timeout({
    data <- matrix(seq_len(24), ncol = 1)
    dc_grid <- seq(2L, 24L, by = 2L)
    method_specs <- list(
      SN = list(),
      WBS = list(M = 6L, wbs_seed = 13L)
    )

    for (method in names(method_specs)) {
      common_args <- c(
        list(
          data = data,
          reg_fun = reg_null,
          cpn_max = 1L,
          dm = 2L,
          cov_rate = 0.6,
          method = method,
          detail = TRUE,
          echo = FALSE,
          dc_grid = dc_grid
        ),
        method_specs[[method]]
      )

      by_cost_mat <- do.call(
        relieverChangepoint::reliever_generic,
        c(common_args, list(cache_backend = "by_cost_mat"))
      )
      by_loss_block <- do.call(
        relieverChangepoint::reliever_generic,
        c(common_args, list(cache_backend = "by_loss_block"))
      )

      expect_same_cpd_path(by_loss_block, by_cost_mat)
      expect_true(all(unlist(by_loss_block$cpd_path$candidates$cpd) %in% dc_grid))
      expect_equal(by_loss_block$cache_profile$objects$int_set$n, length(dc_grid))
    }
  })
})

# Post-search holdout evaluation ---------------------------------------------

test_that("generic NMCD holdout evaluation keeps training CDF cutpoints", {
  with_test_timeout({
    train <- c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8)
    eval <- train + c(rep(0.25, 6L), rep(20, 6L))
    common <- list(
      data = train,
      cpn_max = 1L,
      dm = 3L,
      cov_rate = 1,
      reg_fun = relieverChangepoint::reg_fun_nmcd,
      method = "SN",
      cache_backend = "by_cost_mat"
    )
    default_fit <- do.call(relieverChangepoint::reliever_generic, common)
    explicit_fit <- do.call(
      relieverChangepoint::reliever_generic,
      c(common, list(sort_X = sort(train)))
    )
    default_eval <- relieverChangepoint::evaluate_reliever_segments(
      default_fit, data = train, eval_data = eval
    )
    explicit_eval <- relieverChangepoint::evaluate_reliever_segments(
      explicit_fit, data = train, eval_data = eval
    )

    expect_equal(
      default_eval$candidates$eval_loss,
      explicit_eval$candidates$eval_loss
    )
  })
})

test_that("evaluate_reliever_segments accumulates holdout loss and saved models", {
  with_test_timeout({
    fit_calls <- 0L
    segment_reg <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(2L)
      }
      fit_calls <<- fit_calls + 1L
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      loss <- cbind(abs(y - center), (y - center)^2)
      model <- if (save_model) list(center = center, interval = c(l, r)) else NULL
      list(loss = loss, model = model)
    }

    data <- matrix(sin(seq_len(24) / 3), ncol = 1)
    eval_data <- matrix(data[, 1] + 0.25, ncol = 1)
    eval_index <- seq_len(nrow(data))
    res <- relieverChangepoint::reliever_generic(
      data = data,
      cpn_max = 2L,
      dm = 5L,
      cov_rate = 0.8,
      reg_fun = segment_reg,
      method = "SN",
      run_cpd_ids = 1:2,
      echo = FALSE
    )

    expect_false("models" %in% names(res))
    fit_calls <- 0L

    eval_res <- relieverChangepoint::evaluate_reliever_segments(
      result = res,
      data = data,
      eval_data = eval_data,
      eval_index = eval_index,
      save_model = TRUE
    )

    expect_equal(nrow(eval_res$candidates), 6L)
    expect_true(all(is.finite(eval_res$candidates$eval_loss)))
    expect_true(length(eval_res$models) > 0L)
    expect_true(all(grepl("^[0-9]+:[0-9]+$", names(eval_res$models))))
    expect_true(all(vapply(eval_res$models, function(x) {
      is.list(x)
    }, logical(1))))
    expected_segments <- unique(unlist(lapply(
      eval_res$candidates$cpd,
      function(cpd) {
        cpd <- as.integer(cpd)
        paste(c(1L, cpd + 1L), c(cpd, nrow(data)), sep = ":")
      }
    )))
    expect_equal(fit_calls, length(expected_segments))
    expect_setequal(names(eval_res$models), expected_segments)

    no_model <- relieverChangepoint::evaluate_reliever_segments(
      result = res,
      data = data,
      eval_data = eval_data,
      eval_index = eval_index
    )
    expect_null(no_model$models)
    expect_true(all(is.finite(no_model$candidates$eval_loss)))

    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result = res,
        data = data
      ),
      "eval_data is required"
    )
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result = res,
        data = data,
        eval_data = NULL,
        eval_index = eval_index
      ),
      "eval_data is required"
    )
    automatic_index <- relieverChangepoint::evaluate_reliever_segments(
      result = res,
      data = data,
      eval_data = eval_data
    )
    expect_equal(
      automatic_index$candidates$eval_loss,
      no_model$candidates$eval_loss
    )
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result = res,
        data = data,
        eval_data = eval_data,
        y = data[, 1L],
        eval_y = eval_data[, 1L]
      ),
      "only for a built-in response-based fit"
    )
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result = res,
        data = data,
        eval_data = eval_data[-1L, , drop = FALSE]
      ),
      "eval_index is required"
    )
  })
})

test_that("built-in lasso holdout input mirrors reliever X/y forms", {
  with_test_timeout({
    set.seed(20260719)
    n <- 48L
    p <- 4L
    x <- matrix(stats::rnorm(n * p), nrow = n)
    beta <- c(1, -0.7, 0.4, 0)
    y <- drop(x %*% beta + stats::rnorm(n, sd = 0.2))
    eval_x <- matrix(stats::rnorm(n * p), nrow = n)
    eval_y <- drop(eval_x %*% beta + stats::rnorm(n, sd = 0.2))
    lam_set <- c(0.5, 0.1)

    fit <- relieverChangepoint::reliever(
      X = x, y = y, cpd_family = "lasso",
      cpn_max = 1L, dm = 8L, cov_rate = 0.8,
      method = "SN",
      lam_set = lam_set
    )
    expect_identical(fit$settings$input_spec$type, "response_predictor")
    expect_identical(fit$settings$input_spec$n_predictors, p)
    expect_identical(fit$settings$input_spec$original_form, "separate_xy")

    separate <- relieverChangepoint::evaluate_reliever_segments(
      fit,
      data = x, y = y,
      eval_data = eval_x, eval_y = eval_y
    )
    response_first <- relieverChangepoint::evaluate_reliever_segments(
      fit,
      data = cbind(y, x),
      eval_data = cbind(eval_y, eval_x)
    )
    expect_equal(
      separate$candidates$eval_loss,
      response_first$candidates$eval_loss
    )

    selected_separate <- relieverChangepoint::select_holdout(
      fit,
      data = x, y = y,
      eval_data = eval_x, eval_y = eval_y
    )
    selected_response_first <- relieverChangepoint::select_holdout(
      fit,
      data = cbind(y, x),
      eval_data = cbind(eval_y, eval_x)
    )
    expect_equal(selected_separate$run_id, selected_response_first$run_id)
    expect_equal(selected_separate$K_hat, selected_response_first$K_hat)
    expect_equal(selected_separate$score, selected_response_first$score)

    expect_error(
      relieverChangepoint::select_holdout(
        fit, data = x, eval_data = eval_x
      ),
      "predictor columns but no response"
    )
    expect_error(
      relieverChangepoint::select_holdout(
        fit, data = x, y = y, eval_data = eval_x
      ),
      "Supply both y and eval_y"
    )
    expect_error(
      relieverChangepoint::select_holdout(
        fit,
        data = x[, -1L, drop = FALSE], y = y,
        eval_data = eval_x[, -1L, drop = FALSE], eval_y = eval_y
      ),
      "exactly 4 predictor columns"
    )

    focused <- relieverChangepoint::reliever_lasso(
      cbind(y, x), cpn_max = 1L, dm = 8L, cov_rate = 0.8,
      lam_set = lam_set, method = "SN"
    )
    expect_identical(
      focused$settings$input_spec$original_form, "response_first"
    )
    expect_s3_class(
      relieverChangepoint::select_holdout(
        focused,
        data = cbind(y, x),
        eval_data = cbind(eval_y, eval_x)
      ),
      "reliever_model_selection"
    )
  })
})

test_that("kernel holdout input mirrors raw and precomputed representations", {
  with_test_timeout({
    set.seed(20260719)
    n <- 24L
    x <- matrix(stats::rnorm(n * 2L), nrow = n)
    eval_x <- x + matrix(stats::rnorm(n * 2L, sd = 0.15), nrow = n)
    dist_sq <- relieverChangepoint:::.kernel_dist_sq(x)
    eval_dist_sq <- relieverChangepoint:::.kernel_dist_sq(eval_x, x)
    bandwidth_vec <- c(0.6, 1.2)

    raw_nll <- relieverChangepoint::reliever(
      x, cpd_family = "kde_nll",
      kernel = "student", kernel_args = list(df = 5),
      bandwidth_vec = bandwidth_vec,
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    distance_nll <- relieverChangepoint::reliever(
      dist_sq, cpd_family = "kde_nll",
      kernel = "student", kernel_args = list(df = 5),
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    raw_nll_eval <- relieverChangepoint::evaluate_reliever_segments(
      raw_nll, data = x, eval_data = eval_x
    )
    distance_nll_eval <- relieverChangepoint::evaluate_reliever_segments(
      distance_nll, data = dist_sq, eval_data = eval_dist_sq
    )
    expect_equal(
      raw_nll_eval$candidates$eval_loss,
      distance_nll_eval$candidates$eval_loss,
      tolerance = 1e-12
    )

    generic_nll <- relieverChangepoint::reliever_generic(
      dist_sq,
      reg_fun = relieverChangepoint::reg_fun_kde_nll_solpath,
      var_dim = 2L,
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    generic_nll_eval <- relieverChangepoint::evaluate_reliever_segments(
      generic_nll, data = dist_sq, eval_data = eval_dist_sq
    )
    expect_true(all(is.finite(generic_nll_eval$candidates$eval_loss)))

    raw_l2 <- relieverChangepoint::reliever(
      x, cpd_family = "kde_l2",
      kernel = "matern52", bandwidth = 1.1,
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    kernel_train <- relieverChangepoint:::.kernel_l2_matrix(
      x, kernel = "matern52", bandwidth = 1.1
    )
    kernel_eval <- relieverChangepoint:::.kernel_l2_matrix(
      eval_x, x, kernel = "matern52", bandwidth = 1.1
    )
    precomputed_l2 <- relieverChangepoint::reliever(
      kernel_train, cpd_family = "kde_l2",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    raw_l2_eval <- relieverChangepoint::evaluate_reliever_segments(
      raw_l2, data = x, eval_data = eval_x
    )
    precomputed_l2_eval <- relieverChangepoint::evaluate_reliever_segments(
      precomputed_l2, data = kernel_train, eval_data = kernel_eval
    )
    expect_equal(
      raw_l2_eval$candidates$eval_loss,
      precomputed_l2_eval$candidates$eval_loss,
      tolerance = 1e-12
    )
  }, seconds = 5)
})

test_that("KDE-NLL crossfit holdout selects bandwidth on training rows only", {
  with_test_timeout({
    set.seed(20260720)
    n <- 24L
    x <- matrix(stats::rnorm(n * 2L), nrow = n)
    eval_x <- x + matrix(
      stats::rnorm(n * 2L, sd = 0.15), nrow = n
    )
    dist_sq <- relieverChangepoint:::.kernel_dist_sq(x)
    eval_dist_sq <- relieverChangepoint:::.kernel_dist_sq(eval_x, x)
    bandwidth_vec <- c(0.6, 1.2)

    raw_fit <- relieverChangepoint::reliever(
      x, cpd_family = "kde_nll_crossfit",
      kernel = "laplace",
      bandwidth_vec = bandwidth_vec,
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      nfolds = 2L, fold_type = "blk",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    distance_fit <- relieverChangepoint::reliever(
      dist_sq, cpd_family = "kde_nll_crossfit",
      kernel = "laplace",
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      nfolds = 2L, fold_type = "blk",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    raw_eval <- relieverChangepoint::evaluate_reliever_segments(
      raw_fit, data = x, eval_data = eval_x
    )
    distance_eval <- relieverChangepoint::evaluate_reliever_segments(
      distance_fit, data = dist_sq, eval_data = eval_dist_sq
    )
    expect_equal(
      raw_eval$candidates$eval_loss,
      distance_eval$candidates$eval_loss,
      tolerance = 1e-12
    )

    generic_fit <- relieverChangepoint::reliever_generic(
      dist_sq,
      reg_fun = relieverChangepoint::reg_fun_kde_nll_crossfit,
      var_dim = 2L,
      nfolds = 2L, fold_type = "blk",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85
    )
    generic_eval <- relieverChangepoint::evaluate_reliever_segments(
      generic_fit, data = dist_sq, eval_data = eval_dist_sq
    )
    expect_true(all(is.finite(generic_eval$candidates$eval_loss)))

    training_fit <- relieverChangepoint::reg_fun_kde_nll_crossfit(
      sqrt(dist_sq), l = 1L, r = n,
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      kernel = "laplace", distance_power = 1L,
      nfolds = 2L, fold_type = "blk",
      loss_output_types = c("recv", "crossfit_homo_hyper")
    )
    fixed_cf <- training_fit$loss[
      , 1L + seq_along(bandwidth_vec), drop = FALSE
    ]
    best_bandwidth <- which.min(colSums(fixed_cf))[1L]
    external_fixed <- relieverChangepoint:::.kde_nll_loss(
      sqrt(eval_dist_sq), bandwidth_vec, var_dim = 2L,
      kernel = "laplace", kernel_args = list(),
      distance_power = 1L
    )
    expected_by_output <- c(
      sum(external_fixed[, best_bandwidth]),
      sum(external_fixed[, best_bandwidth]),
      colSums(external_fixed)
    )
    no_change <- raw_eval$candidates[
      raw_eval$candidates$K == 0L, , drop = FALSE
    ]
    output_id <- raw_fit$run_meta$loss_output_id[
      match(no_change$run_id, raw_fit$run_meta$run_id)
    ]
    expect_equal(
      no_change$eval_loss,
      expected_by_output[output_id],
      tolerance = 1e-12
    )

    raw_selected <- relieverChangepoint::select_holdout(
      raw_fit, data = x, eval_data = eval_x
    )
    distance_selected <- relieverChangepoint::select_holdout(
      distance_fit, data = dist_sq, eval_data = eval_dist_sq
    )
    expect_equal(raw_selected$run_id, distance_selected$run_id)
    expect_equal(raw_selected$K_hat, distance_selected$K_hat)
    expect_equal(raw_selected$cpd_hat, distance_selected$cpd_hat)
    expect_equal(raw_selected$score, distance_selected$score, tolerance = 1e-12)
  }, seconds = 5)
})

test_that("evaluate_reliever_segments reuses reg_fun_lasso_solpath on eval_data", {
  with_test_timeout({
    set.seed(20240711)
    n <- 36L
    x <- matrix(stats::rnorm(n * 4L), nrow = n)
    beta <- c(1.0, -0.7, 0, 0)
    y <- drop(x %*% beta + stats::rnorm(n, sd = 0.2))
    data <- cbind(y, x)
    eval_data <- data
    eval_data[, 1L] <- eval_data[, 1L] + stats::rnorm(n, sd = 0.05)
    lam_set <- c(0.5, 0.1)

    res <- relieverChangepoint::reliever_generic(
      data = data,
      cpn_max = 1L,
      dm = 6L,
      cov_rate = 0.8,
      reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
      method = "SN",
      run_cpd_ids = seq_along(lam_set),
      lam_set = lam_set,
      thresh = 1e-7,
      echo = FALSE
    )
    eval_res <- relieverChangepoint::evaluate_reliever_segments(
      result = res,
      data = data,
      eval_data = eval_data,
      eval_index = seq_len(n)
    )

    expect_equal(nrow(eval_res$candidates), nrow(res$cpd_path$candidates))
    expect_true(all(is.finite(eval_res$candidates$eval_loss)))
    expect_false(".reliever_loss_spec" %in% names(res))

    restored <- unserialize(serialize(res, NULL))
    restored_eval <- relieverChangepoint::evaluate_reliever_segments(
      result = restored,
      data = data,
      eval_data = eval_data,
      eval_index = seq_len(n)
    )
    expect_equal(
      restored_eval$candidates$eval_loss,
      eval_res$candidates$eval_loss
    )

    threshold_eval <- relieverChangepoint::evaluate_reliever_segments(
      result = res,
      data = data,
      eval_data = eval_data,
      eval_index = seq_len(n),
      thresh = 1e-6
    )
    expect_true(all(is.finite(threshold_eval$candidates$eval_loss)))

    tagged_lasso <- function(data, l, r, l_end = l, r_end = r,
                             save_model = FALSE,
                             is_virtual_run = FALSE, ...) {
      out <- relieverChangepoint::reg_fun_lasso_solpath(
        data = data, l = l, r = r, l_end = l_end, r_end = r_end,
        save_model = save_model, is_virtual_run = is_virtual_run, ...
      )
      if (is_virtual_run) {
        out$loss_output_meta$display_tag <- "external"
      }
      out
    }
    tagged_eval <- relieverChangepoint::evaluate_reliever_segments(
      result = res,
      data = data,
      eval_data = eval_data,
      eval_index = seq_len(n),
      reg_fun = tagged_lasso
    )
    expect_true(all(is.finite(tagged_eval$candidates$eval_loss)))

    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result = res,
        data = data,
        eval_data = eval_data,
        eval_index = seq_len(n),
        lam_set = c(0.6, 0.1)
      ),
      "same run structure and hyperparameter path"
    )
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result = res,
        data = data[-1L, , drop = FALSE],
        eval_data = eval_data[-1L, , drop = FALSE],
        eval_index = seq_len(n - 1L)
      ),
      "same number of observations"
    )
  })
})

test_that("native mean and twostep results retain their evaluation loss", {
  with_test_timeout({
    set.seed(20260716)
    data <- rbind(
      matrix(stats::rnorm(24L, mean = 0), 12L, 2L),
      matrix(stats::rnorm(24L, mean = 2), 12L, 2L),
      matrix(stats::rnorm(24L, mean = -2), 12L, 2L)
    )
    eval_data <- data + matrix(
      stats::rnorm(length(data), sd = 0.1), nrow = nrow(data)
    )

    mean_fit <- relieverChangepoint::reliever_mean(
      data, cpn_max = 2L, dm = 6L, cov_rate = 0.8, cpn_crit = "none"
    )
    mean_eval <- relieverChangepoint::evaluate_reliever_segments(
      mean_fit,
      data = data,
      eval_data = eval_data,
      eval_index = seq_len(nrow(data))
    )
    expect_true(all(is.finite(mean_eval$candidates$eval_loss)))

    twostep_fit <- relieverChangepoint::twostep(
      data, reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 2L, dm = 6L, method = "BS"
    )
    twostep_eval <- relieverChangepoint::evaluate_reliever_segments(
      twostep_fit,
      data = data,
      eval_data = eval_data,
      eval_index = seq_len(nrow(data))
    )
    expect_true(all(is.finite(twostep_eval$candidates$eval_loss)))
  })
})

test_that("select_holdout evaluates and selects by holdout loss", {
  with_test_timeout({
    segment_reg <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 2L,
          loss_output_meta = data.frame(
            loss_output_id = 1:2,
            row_type = "toy",
            hyper_id = 1:2,
            hyper_value = I(as.list(c(0.1, 0.2)))
          )
        ))
      }
      len <- r - l + 1L
      n_eval <- r_end - l_end + 1L
      loss <- cbind(
        rep((len - 6L)^2 + 1, n_eval),
        rep((len - 4L)^2, n_eval)
      )
      list(loss = loss, model = NULL)
    }
    cpd_rows <- list(
      integer(), 6L, c(4L, 8L),
      integer(), 6L, c(4L, 8L)
    )
    result <- list(
      cpd_path = list(
        select_by = "K",
        selector_map = NULL,
        candidates = data.frame(
          run_id = rep(1:2, each = 3L),
          candidate_id = rep(1:3, times = 2L),
          K = rep(0:2, times = 2L),
          loss = 0,
          cpd = I(cpd_rows)
        )
      ),
      run_meta = data.frame(
        run_id = 1:2,
        loss_output_id = 1:2,
        row_type = "toy",
        hyper_id = 1:2,
        hyper_value = I(as.list(c(0.1, 0.2)))
      )
    )
    data <- matrix(seq_len(12), ncol = 1)
    eval_data <- data
    eval_index <- seq_len(nrow(data))

    expect_error(
      relieverChangepoint::select_holdout(
        result, data = data,
        eval_data = eval_data, eval_index = eval_index
      ),
      "reg_fun is unavailable"
    )
    selected <- relieverChangepoint::select_holdout(
      result, data = data,
      eval_data = eval_data, eval_index = eval_index,
      reg_fun = segment_reg
    )
    expect_s3_class(selected, "reliever_model_selection")
    expect_equal(selected$rule, "holdout")
    expect_equal(selected$K_hat, 2L)
    expect_equal(selected$run_id, 2L)
    expect_equal(selected$cpd_hat[[1L]], c(4L, 8L))
    expect_equal(selected$score, 0)
    expect_equal(selected$hyper_value, 0.2)

    selected_k1 <- relieverChangepoint::select_holdout(
      result, data = data,
      eval_data = eval_data, eval_index = eval_index,
      reg_fun = segment_reg, K = 1L
    )
    expect_equal(selected_k1$K_hat, 1L)
    expect_equal(selected_k1$run_id, 1L)
    expect_equal(selected_k1$cpd_hat[[1L]], 6L)
    expect_equal(selected_k1$hyper_value, 0.1)
    expect_error(
      relieverChangepoint::select_holdout(
        result, data = data, reg_fun = segment_reg
      ),
      "eval_data is required"
    )
    expect_error(
      relieverChangepoint::select_holdout(
        result, data = data,
        eval_data = NULL,
        eval_index = eval_index,
        reg_fun = segment_reg
      ),
      "eval_data is required"
    )
    selected_automatic_index <- relieverChangepoint::select_holdout(
      result, data = data, eval_data = eval_data,
      reg_fun = segment_reg
    )
    expect_equal(selected_automatic_index$run_id, selected$run_id)
    expect_equal(selected_automatic_index$K_hat, selected$K_hat)
  })
})

test_that("evaluate_reliever_segments handles non-one-to-one eval samples", {
  with_test_timeout({
    segment_reg <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 2L,
          loss_output_meta = data.frame(
            loss_output_id = 1:2,
            row_type = c("absolute", "squared")
          )
        ))
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      loss <- cbind(abs(y - center), (y - center)^2)
      list(loss = loss, model = NULL)
    }
    segment_score <- function(data, eval_data, eval_index, left, right,
                              loss_output_id) {
      idx <- which(eval_index >= left & eval_index <= right)
      if (length(idx) == 0L) {
        return(0)
      }
      center <- mean(data[left:right, 1])
      y <- eval_data[idx, 1]
      if (loss_output_id == 1L) {
        sum(abs(y - center))
      } else {
        sum((y - center)^2)
      }
    }
    expected_column <- function(data, eval_data, eval_index, bounds,
                                loss_output_id) {
      sum(vapply(seq_len(nrow(bounds)), function(i) {
        segment_score(
          data, eval_data, eval_index,
          bounds[i, 1L], bounds[i, 2L], loss_output_id
        )
      }, numeric(1)))
    }

    data <- matrix(seq_len(12) / 10, ncol = 1)
    eval_index <- c(10L, 2L, 2L, 11L, 5L)
    eval_data <- matrix(c(2.0, -1.5, -1.1, 0.4, 3.0), ncol = 1)
    result <- list(
      cpd_path = list(
        select_by = "K",
        selector_map = NULL,
        candidates = data.frame(
          run_id = rep(1:2, each = 3L),
          candidate_id = rep(1:3, times = 2L),
          K = rep(0:2, times = 2L),
          loss = 0,
          cpd = I(list(
            integer(), 6L, c(4L, 8L),
            integer(), 6L, c(4L, 8L)
          ))
        )
      ),
      run_meta = data.frame(
        run_id = 1:2,
        row_type = c("absolute", "squared"),
        loss_output_id = 1:2
      )
    )

    got <- relieverChangepoint::evaluate_reliever_segments(
      result, data = data,
      eval_data = eval_data,
      eval_index = eval_index,
      reg_fun = segment_reg
    )$candidates$eval_loss
    expected <- numeric(6L)
    for (loss_output_id in 1:2) {
      row_offset <- (loss_output_id - 1L) * 3L
      expected[row_offset + 1L] <- expected_column(
        data, eval_data, eval_index,
        matrix(c(1L, 12L), ncol = 2L), loss_output_id
      )
      expected[row_offset + 2L] <- expected_column(
        data, eval_data, eval_index,
        matrix(c(1L, 6L, 7L, 12L), ncol = 2L, byrow = TRUE),
        loss_output_id
      )
      expected[row_offset + 3L] <- expected_column(
        data, eval_data, eval_index,
        matrix(c(1L, 4L, 5L, 8L, 9L, 12L), ncol = 2L, byrow = TRUE),
        loss_output_id
      )
    }

    expect_equal(got, expected, tolerance = 1e-12)
    run_two <- relieverChangepoint::evaluate_reliever_segments(
      result, data = data,
      eval_data = eval_data,
      eval_index = eval_index,
      reg_fun = segment_reg,
      run_ids = 2L
    )$candidates
    expect_equal(unique(run_two$run_id), 2L)
    expect_equal(run_two$eval_loss, expected[4:6], tolerance = 1e-12)
    duplicate_two <- relieverChangepoint::evaluate_reliever_segments(
      result, data = data,
      eval_data = eval_data,
      eval_index = eval_index,
      reg_fun = segment_reg,
      run_ids = c(2L, 2L)
    )$candidates
    expect_equal(duplicate_two, run_two)
    squared <- relieverChangepoint::evaluate_reliever_segments(
      result, data = data,
      eval_data = eval_data,
      eval_index = eval_index,
      reg_fun = segment_reg,
      run_type = "squared"
    )$candidates
    expect_equal(squared, run_two)
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result, data = data,
        eval_data = eval_data,
        eval_index = eval_index,
        reg_fun = segment_reg,
        run_ids = 2L,
        run_type = "squared"
      ),
      "only one of run_ids and run_type"
    )
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result, data = data,
        eval_data = eval_data,
        eval_index = c(1.5, 2, 3, 4, 5),
        reg_fun = segment_reg
      ),
      "integer original sample index"
    )
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result, data = data,
        eval_data = eval_data,
        eval_index = eval_index,
        reg_fun = segment_reg,
        run_ids = 1.2
      ),
      "run_ids"
    )
    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        result, data = data,
        eval_data = eval_data,
        eval_index = eval_index,
        reg_fun = segment_reg,
        data_stack_fun = function(data, eval_data) {
          rbind(data, eval_data, eval_data[1L, , drop = FALSE])
        }
      ),
      "exactly n \\+ n_eval rows"
    )
  })
})

# User-facing validation ------------------------------------------------------

test_that("reliever validates simple argument errors early", {
  with_test_timeout({
    data <- matrix(seq_len(12), ncol = 1)

    expect_error(
      relieverChangepoint::reliever_generic(data, cpn_max = 1L, dm = 0L, reg_fun = reg_null),
      "dm must be"
    )
    expect_error(
      relieverChangepoint::reliever_generic(data, cpn_max = 1.5, dm = 4L, reg_fun = reg_null),
      "cpn_max must be"
    )
    expect_error(
      relieverChangepoint::reliever_generic(data, cpn_max = 1L, dm = 4L, cov_rate = 1.1,
                    reg_fun = reg_null),
      "cov_rate"
    )
    expect_error(
      relieverChangepoint::reliever_generic(data, cpn_max = 1L, dm = 4L, method = "bad",
                    reg_fun = reg_null),
      "should be one of"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        data, cpn_max = 1L, dm = 4L, reg_fun = reg_null,
        run_cpd_ids = 1.2
      ),
      "run_cpd_ids"
    )
    deduplicated_outputs <- relieverChangepoint::reliever_generic(
      data, cpn_max = 1L, dm = 4L, reg_fun = reg_null,
      run_cpd_ids = c(1L, 1L)
    )
    expect_identical(deduplicated_outputs$run_meta$loss_output_id, 1L)
    expect_error(
      relieverChangepoint::reliever_generic(
        data, cpn_max = 1L, dm = 4L, reg_fun = reg_null,
        owner_key = NA
      ),
      "owner_key"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        data = data, reg_fun = reg_null,
        cpn_max = 1L, dm = 4L,
        dc_grid = c(2L, 2L, 12L)
      ),
      "dc_grid"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        data = data, reg_fun = reg_null,
        cpn_max = 1L, dm = 4L,
        dc_grid = c(2.5, 12)
      ),
      "dc_grid"
    )
  })
})

test_that("reliever reports reg_fun protocol violations clearly", {
  with_test_timeout({
    data <- matrix(seq_len(12), ncol = 1)
    bad_virtual <- function(data, l, r, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return("one")
      }
      list(loss = matrix(0, r - l + 1L, 1L))
    }
    noninteger_virtual <- function(data, l, r, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1.5)
      }
      list(loss = matrix(0, r - l + 1L, 1L))
    }
    missing_loss <- function(data, l, r, l_end = l, r_end = r,
                             is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      list()
    }
    wrong_outputs <- function(data, l, r, l_end = l, r_end = r,
                              is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(2L)
      }
      list(loss = matrix(0, r_end - l_end + 1L, 1L))
    }
    invalid_default <- function(data, l, r, l_end = l, r_end = r,
                                is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 2L,
          loss_output_meta = data.frame(
            default_selection = c(TRUE, NA)
          )
        ))
      }
      list(loss = matrix(0, r_end - l_end + 1L, 2L))
    }
    no_virtual_flag <- function(data, l, r, l_end = l, r_end = r,
                                save_model = FALSE) {
      list(loss = matrix(0, r_end - l_end + 1L, 1L))
    }

    expect_error(
      relieverChangepoint::reliever_generic(
        data, cpn_max = 1L, dm = 4L, reg_fun = NULL
      ),
      "reg_fun must be a function"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        data, cpn_max = 1L, dm = 4L, reg_fun = "reg_fun_mean"
      ),
      "reg_fun must be a function"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        data, cpn_max = 1L, dm = 4L, reg_fun = no_virtual_flag
      ),
      "is_virtual_run"
    )
    expect_error(
      relieverChangepoint::reliever_generic(data, cpn_max = 1L, dm = 4L, reg_fun = bad_virtual),
      "virtual_run"
    )
    expect_error(
      relieverChangepoint::reliever_generic(data, cpn_max = 1L, dm = 4L, reg_fun = noninteger_virtual),
      "virtual_run"
    )
    expect_error(
      relieverChangepoint::reliever_generic(data, cpn_max = 1L, dm = 4L, reg_fun = missing_loss),
      "containing `loss`"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        data,
        cpn_max = 1L,
        dm = 4L,
        reg_fun = wrong_outputs
      ),
      "virtual call reported 2"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        data, cpn_max = 1L, dm = 4L, reg_fun = invalid_default
      ),
      "default_selection"
    )
  })
})

test_that("a single-output reg_fun may return a numeric loss vector", {
  with_test_timeout({
    data <- matrix(seq_len(12), ncol = 1L)
    vector_loss <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE,
                            is_virtual_run = FALSE) {
      if (is_virtual_run) {
        return(1L)
      }
      center <- mean(data[l:r, 1L])
      list(loss = (data[l_end:r_end, 1L] - center)^2)
    }
    fit <- relieverChangepoint::reliever_generic(
      data = data, reg_fun = vector_loss,
      cpn_max = 1L, dm = 4L, cov_rate = 0.8, method = "SN"
    )
    expect_s3_class(fit, "reliever_result")
    expect_true(all(is.finite(fit$cpd_path$candidates$loss)))
  })
})

test_that("custom reg_fun metadata controls fit-time selection eligibility", {
  with_test_timeout({
    selective_reg <- function(data, l, r, l_end = l, r_end = r,
                              save_model = FALSE,
                              is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 2L,
          loss_output_meta = data.frame(
            row_type = c("candidate_a", "candidate_b"),
            hyper_value = c(1, 2),
            default_selection = c(FALSE, TRUE),
            loss_kind = "pinball_loss"
          )
        ))
      }
      n_eval <- r_end - l_end + 1L
      list(loss = matrix(1, nrow = n_eval, ncol = 2L))
    }
    fit <- relieverChangepoint::reliever_generic(
      matrix(seq_len(16), ncol = 1L),
      reg_fun = selective_reg,
      cpn_max = 1L, dm = 4L, cov_rate = 0.8,
      method = "SN", cpn_crit = "loss"
    )

    expect_equal(nrow(fit$summary), 1L)
    expect_equal(fit$summary$hyper_value, 2)
    expect_equal(nrow(fit$run_meta), 2L)
    expect_identical(fit$run_meta$loss_kind, rep("pinball_loss", 2L))
    expect_equal(
      nrow(relieverChangepoint::select_by_run(
        fit, cpn_crit = "loss"
      )),
      2L
    )

    suppressed_reg <- function(data, l, r, l_end = l, r_end = r,
                               save_model = FALSE,
                               is_virtual_run = FALSE, ...) {
      out <- selective_reg(
        data, l, r, l_end, r_end, save_model, is_virtual_run, ...
      )
      if (is_virtual_run) {
        out$loss_output_meta$default_selection <- FALSE
      }
      out
    }
    suppressed <- relieverChangepoint::reliever_generic(
      matrix(seq_len(16), ncol = 1L),
      reg_fun = suppressed_reg,
      cpn_max = 1L, dm = 4L, cov_rate = 0.8,
      method = "SN", cpn_crit = "loss"
    )
    expect_equal(nrow(suppressed$summary), 0L)
    expect_equal(nrow(suppressed$cpd_path$candidates), 4L)
  })
})
