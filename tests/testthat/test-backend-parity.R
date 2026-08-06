test_that("SN by_loss_block cache backend matches by_cost_mat reg_fun_lasso_solpath results", {
  with_test_timeout({
    set.seed(20240523)

    n <- 50L
    p <- 8L
    tau <- as.integer(c(0.34, 0.68) * n)
    sparsity <- 2L

    b0 <- numeric(p)
    b0[seq_len(sparsity)] <- test_scale_snr(
      stats::rnorm(sparsity), snr = 2
    )

    delta <- matrix(0, p, length(tau))
    delta[seq_len(sparsity), ] <- matrix(
      stats::rnorm(sparsity * length(tau)),
      nrow = sparsity
    )
    delta[seq_len(sparsity), ] <- apply(
      delta[seq_len(sparsity), , drop = FALSE],
      2,
      test_scale_snr,
      snr = 0.5
    )

    data <- relieverChangepoint::dgp_linear_regression(n, p, tau, b0, delta)$data
    lambda_set <- glmnet::glmnet(data[, -1], data[, 1], nlambda = 3)$lambda
    lambda_set <- lambda_set * sqrt(n)
    run_cpd_ids <- seq_along(lambda_set)

    common_args <- list(
      data = data,
      cpn_max = length(tau),
      dm = 8L,
      cov_rate = 0.8,
      reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
      method = "SN",
      run_cpd_ids = run_cpd_ids,
      lam_set = lambda_set,
      detail = TRUE,
      echo = FALSE
    )

    by_cost_mat <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_cost_mat"))
    )
    by_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_loss_block"))
    )

    expect_equal(by_loss_block$settings$cache_backend, "by_loss_block")
    expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-8)
  })
})

test_that("by_loss_block native runner skips cache state when detail is not requested", {
  with_test_timeout({
    n <- 24L
    d <- 5L
    int_set <- relieverChangepoint::create_relief_itv(n, cov_rate = 0.8, dm = d)
    reg_fun_wrap <- function(l, r, l_end, r_end) {
      matrix(0, nrow = r_end - l_end + 1L, ncol = 1L)
    }

    compact <- relieverChangepoint:::cpd_r_by_loss_block(
      "SN", n, 1L, d, matrix(integer(), 0L, 2L), numeric(), 0,
      1L, reg_fun_wrap,
      int_set$miss_cover_len, int_set$int_len, int_set$layer_point,
      int_set$int_eps,
      cache_state = NULL,
      return_cache_profile = FALSE
    )
    detailed <- relieverChangepoint:::cpd_r_by_loss_block(
      "SN", n, 1L, d, matrix(integer(), 0L, 2L), numeric(), 0,
      1L, reg_fun_wrap,
      int_set$miss_cover_len, int_set$int_len, int_set$layer_point,
      int_set$int_eps,
      cache_state = NULL,
      return_cache_profile = TRUE
    )

    expect_null(compact$cache_profile)
    expect_false(is.null(detailed$cache_profile))
    expect_false(is.null(detailed$cache_profile$cache_state))
    expect_equal(cpp_loss(compact), cpp_loss(detailed))
    expect_equal(cpp_cps_num(compact), cpp_cps_num(detailed))
  })
})

test_that("by_loss_block full-search requests switch to by_cost_mat", {
  with_test_timeout({
    data <- matrix(seq_len(24), ncol = 1)

    res <- NULL
    expect_warning(
      res <- relieverChangepoint::reliever_generic(
        data,
        cpn_max = 1L,
        dm = 4L,
        cov_rate = 1,
        reg_fun = reg_null,
        method = "SN",
        cache_backend = "by_loss_block",
        detail = TRUE
      ),
      "Set cache_backend = \"by_cost_mat\" explicitly"
    )
    expect_equal(res$settings$cache_backend, "by_cost_mat")
    expect_true(is.matrix(res$cache_profile$objects$cost_mat))
    expect_true(all(is.finite(cpd_path_loss(res))))

    loss_block_profile <- relieverChangepoint::reliever_generic(
      data,
      cpn_max = 1L,
      dm = 4L,
      cov_rate = 0.8,
      reg_fun = reg_null,
      method = "SN",
      cache_backend = "by_loss_block",
      detail = TRUE
    )$cache_profile
    cache_warnings <- character()
    reused <- withCallingHandlers(
      relieverChangepoint::reliever_generic(
        data,
        cpn_max = 1L,
        dm = 4L,
        cov_rate = 1,
        reg_fun = reg_null,
        method = "SN",
        cache_backend = "by_loss_block",
        cache_profile = loss_block_profile,
        detail = TRUE
      ),
      warning = function(w) {
        cache_warnings <<- c(cache_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_true(any(grepl("using cache_backend = \"by_cost_mat\"", cache_warnings)))
    expect_true(any(grepl(
      "Omit cache_profile or supply a profile", cache_warnings
    )))
    expect_equal(reused$settings$cache_backend, "by_cost_mat")
    expect_true(is.matrix(reused$cache_profile$objects$cost_mat))
    expect_true(all(is.finite(cpd_path_loss(reused))))
  })
})

test_that("by_loss_block cache backend matches by_cost_mat PELT and OP results", {
  with_test_timeout({
    mean_reg <- function(data, l, r, l_end = l, r_end = r,
                         is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(2L)
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      list(loss = cbind((y - center)^2, abs(y - center)))
    }

    data <- matrix(
      c(rep(0, 10), rep(2, 10), rep(-1, 10)) + seq_len(30) / 100,
      ncol = 1
    )
    method_specs <- list(
      PELT = list(pen_val = c(0.1, 1)),
      OP = list(pen_val = c(0.1, 1))
    )

    for (method in names(method_specs)) {
      common_args <- c(
        list(
          data = data,
          cpn_max = 2L,
          dm = 5L,
          cov_rate = 0.85,
          reg_fun = mean_reg,
          method = method,
          run_cpd_ids = 1:2,
          detail = TRUE,
          echo = FALSE
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

      expect_equal(by_loss_block$settings$cache_backend, "by_loss_block",
                   info = method)
      expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-8)
      expect_equal(by_loss_block$diagnostics$num_pruned,
                   by_cost_mat$diagnostics$num_pruned,
                   info = method)
    }
  })
})

test_that("by_loss_block cache backend matches by_cost_mat WBS-family null results", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1)
    method_specs <- list(
      BS = list(M = 0L),
      WBS = list(M = 8L, wbs_seed = 11L),
      WBS_recursive = list(M = 8L, wbs_seed = 11L),
      SeedBS = list(M = 8L)
    )

    for (method in names(method_specs)) {
      common_args <- c(
        list(
          data = data,
          cpn_max = 2L,
          dm = 5L,
          cov_rate = 0.85,
          reg_fun = reg_null,
          method = method,
          detail = TRUE,
          echo = FALSE
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

    expect_equal(by_loss_block$settings$cache_backend, "by_loss_block")
      expect_same_cpd_path(by_loss_block, by_cost_mat)
    }
  })
})

test_that("WBS-family tie-breaking is backend-stable for nonzero flat losses", {
  with_test_timeout({
    flat_loss_reg <- function(data, l, r, l_end = l, r_end = r,
                              is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      list(loss = matrix(data[l_end:r_end, 1]^2, ncol = 1L))
    }

    set.seed(20260526)
    data <- matrix(stats::rnorm(48), ncol = 1)
    method_specs <- list(
      BS = list(M = 0L),
      WBS = list(M = 8L, wbs_seed = 20260526L),
      SeedBS = list(M = 8L)
    )

    for (method in names(method_specs)) {
      common_args <- c(
        list(
          data = data,
          cpn_max = 2L,
          dm = 6L,
          cov_rate = 0.8,
          reg_fun = flat_loss_reg,
          method = method,
          run_cpd_ids = 1L,
          echo = FALSE
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

      expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-8)
    }
  })
})

test_that("WBS-family by_loss_block uses reliever costs for final non-null loss", {
  with_test_timeout({
    mean_reg <- function(data, l, r, l_end = l, r_end = r,
                         is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(1L)
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      list(loss = matrix((y - center)^2, ncol = 1L))
    }

    data <- matrix(
      c(rep(0, 9), rep(2.5, 10), rep(-1, 9)) + seq_len(28) / 50,
      ncol = 1
    )
    method_specs <- list(
      BS = list(M = 0L),
      WBS = list(M = 8L, wbs_seed = 17L),
      SeedBS = list(M = 8L)
    )

    for (method in names(method_specs)) {
      common_args <- c(
        list(
          data = data,
          cpn_max = 2L,
          dm = 5L,
          cov_rate = 0.85,
          reg_fun = mean_reg,
          method = method,
          detail = TRUE,
          echo = FALSE
        ),
        method_specs[[method]]
      )

      by_cost_mat <- do.call(
        relieverChangepoint::reliever_generic,
        c(common_args, list(cache_backend = "by_cost_mat"))
      )
      by_loss_block <- expect_warning(
        do.call(
          relieverChangepoint::reliever_generic,
          c(common_args, list(cache_backend = "by_loss_block"))
        ),
        NA
      )

      expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-8)
      expect_true(any(cpd_path_loss(by_loss_block) != 0))
    }
  })
})

test_that("by_loss_block warns when an exact interval has no relief owner", {
  with_test_timeout({
    reg_fun_wrap <- function(l, r, l_end, r_end) {
      matrix(0, nrow = r_end - l_end + 1L, ncol = 1L)
    }

    res <- NULL
    expect_warning(
      res <- relieverChangepoint:::cpd_r_by_loss_block(
        "SN", 20L, 1L, 4L, matrix(integer(), 0L, 2L), numeric(), 0,
        1L, reg_fun_wrap,
        as.integer(100L), as.integer(4L), as.integer(1L),
        matrix(as.integer(c(0L, 4L)), ncol = 2L),
        cache_state = NULL,
        return_cache_profile = TRUE,
        use_owner_key = TRUE
      ),
      "please report this warning with a reproducible example"
    )

    expect_gt(res$cache_profile$cache_state$full_update_calls, 0)
  })
})

test_that("by_loss_block results support post-CPD holdout evaluation", {
  with_test_timeout({
    segment_reg <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE, is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(2L)
      }
      center <- mean(data[l:r, 1])
      y <- data[l_end:r_end, 1]
      loss <- cbind(abs(y - center), (y - center)^2)
      model <- if (save_model) list(center = center, interval = c(l, r)) else NULL
      list(loss = loss, model = model)
    }

    data <- matrix(sin(seq_len(24) / 3), ncol = 1)
    eval_data <- matrix(data[, 1] + 0.25, ncol = 1)
    eval_index <- seq_len(nrow(data))
    common_args <- list(
      data = data,
      cpn_max = 2L,
      dm = 5L,
      cov_rate = 0.8,
      reg_fun = segment_reg,
      method = "SN",
      run_cpd_ids = 1:2,
      echo = FALSE
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
    expect_equal(by_loss_block$run_meta, by_cost_mat$run_meta)
    expect_false("models" %in% names(by_loss_block))

    by_cost_mat_eval <- relieverChangepoint::evaluate_reliever_segments(
      by_cost_mat, data = data,
      eval_data = eval_data, eval_index = eval_index, save_model = TRUE
    )
    by_loss_block_eval <- relieverChangepoint::evaluate_reliever_segments(
      by_loss_block, data = data,
      eval_data = eval_data, eval_index = eval_index, save_model = TRUE
    )

    expect_equal(by_loss_block_eval$candidates$eval_loss,
                 by_cost_mat_eval$candidates$eval_loss)
    expect_length(by_loss_block_eval$models,
                  length(by_cost_mat_eval$models))
    expect_true(length(by_loss_block_eval$models) > 0L)
    expect_true(all(vapply(by_loss_block_eval$models, function(x) {
      is.list(x)
    }, logical(1))))

    dc_grid <- seq(2L, 24L, by = 2L)
    by_cost_mat_grid <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_cost_mat", dc_grid = dc_grid))
    )
    by_loss_block_grid <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_loss_block", dc_grid = dc_grid))
    )
    expect_same_cpd_path(by_loss_block_grid, by_cost_mat_grid)

    by_cost_mat_grid_eval <- relieverChangepoint::evaluate_reliever_segments(
      by_cost_mat_grid, data = data,
      eval_data = eval_data, eval_index = eval_index,
      save_model = TRUE
    )
    by_loss_block_grid_eval <- relieverChangepoint::evaluate_reliever_segments(
      by_loss_block_grid, data = data,
      eval_data = eval_data, eval_index = eval_index,
      save_model = TRUE
    )

    expect_equal(by_loss_block_grid_eval$candidates$eval_loss,
                 by_cost_mat_grid_eval$candidates$eval_loss)
    expect_equal(by_cost_mat_grid_eval$models[[1L]]$interval,
                 c(1L, nrow(data)))
    expect_equal(by_loss_block_grid_eval$models[[1L]]$interval,
                 c(1L, nrow(data)))
    expect_true(all(unlist(by_loss_block_grid$cpd_path$candidates$cpd) %in% dc_grid))
  })
})

test_that("SN by_loss_block cache backend matches by_cost_mat reg_fun_lasso_crossfit costs", {
  with_test_timeout({
    set.seed(20240526)

    n <- 40L
    p <- 6L
    x <- matrix(stats::rnorm(n * p), n, p)
    beta <- c(1.2, -0.8, rep(0, p - 2L))
    data <- cbind(as.numeric(x %*% beta + stats::rnorm(n, sd = 0.5)), x)
    lam_set <- c(0.05, 0.2)
    run_cpd_ids <- seq_len(2L + length(lam_set))

    common_args <- list(
      data = data,
      cpn_max = 1L,
      dm = 10L,
      cov_rate = 0.8,
      reg_fun = relieverChangepoint::reg_fun_lasso_crossfit,
      method = "SN",
      run_cpd_ids = run_cpd_ids,
      lam_set = lam_set,
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      nfolds = 2L,
      echo = FALSE
    )

    by_cost_mat <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_cost_mat"))
    )
    by_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(common_args, list(cache_backend = "by_loss_block"))
    )

    expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-8)
    expect_equal(
      by_cost_mat$run_meta$row_type,
      c("recv", "incv", "crossfit_homo_hyper", "crossfit_homo_hyper")
    )
    expect_equal(
      by_loss_block$run_meta,
      by_cost_mat$run_meta
    )
    expect_null(by_loss_block$interval_meta)
  })
})

test_that("SN by_loss_block cache backend matches by_cost_mat generic ReCV losses", {
  with_test_timeout({
    set.seed(20260614)

    mean_data <- matrix(stats::rnorm(32L * 3L), nrow = 32L)
    mean_args <- list(
      data = mean_data,
      cpn_max = 1L,
      dm = 8L,
      cov_rate = 0.85,
      reg_fun = relieverChangepoint::reg_fun_mean_crossfit,
      method = "SN",
      nfolds = 2L,
      fold_type = "blk",
      echo = FALSE
    )

    mean_cost_mat <- do.call(relieverChangepoint::reliever_generic, mean_args)
    mean_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(mean_args, list(cache_backend = "by_loss_block"))
    )

    expect_same_cpd_path(mean_loss_block, mean_cost_mat, tolerance = 1e-10)
    expect_equal(
      mean_loss_block$run_meta,
      mean_cost_mat$run_meta
    )

    x <- matrix(stats::rnorm(26L * 2L), nrow = 26L)
    dist_sq <- as.matrix(stats::dist(x))^2
    kde_args <- list(
      data = dist_sq,
      cpn_max = 1L,
      dm = 7L,
      cov_rate = 0.8,
      reg_fun = relieverChangepoint::reg_fun_kde_nll_crossfit,
      method = "SN",
      nfolds = 2L,
      bandwidth_vec = c(0.8, 1.3),
      var_dim = 2L,
      fold_type = "blk",
      echo = FALSE
    )

    kde_cost_mat <- do.call(relieverChangepoint::reliever_generic, kde_args)
    kde_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(kde_args, list(cache_backend = "by_loss_block"))
    )

    expect_same_cpd_path(kde_loss_block, kde_cost_mat, tolerance = 1e-10)
    expect_equal(
      kde_loss_block$run_meta,
      kde_cost_mat$run_meta
    )
  })
})

test_that("SN cache backends match nonparametric solpath reg_fun", {
  with_test_timeout({
    set.seed(20260711)
    x <- matrix(stats::rnorm(24L * 2L), nrow = 24L)
    dist_sq <- as.matrix(stats::dist(x))^2
    kde_args <- list(
      data = dist_sq,
      cpn_max = 1L,
      dm = 6L,
      cov_rate = 0.85,
      reg_fun = relieverChangepoint::reg_fun_kde_nll_solpath,
      method = "SN",
      run_cpd_ids = 1:2,
      bandwidth_vec = c(0.75, 1.15),
      var_dim = 2L,
      echo = FALSE
    )

    kde_cost_mat <- do.call(relieverChangepoint::reliever_generic, kde_args)
    kde_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(kde_args, list(cache_backend = "by_loss_block"))
    )

    expect_same_cpd_path(kde_loss_block, kde_cost_mat, tolerance = 1e-10)
    expect_equal(
      kde_loss_block$run_meta,
      kde_cost_mat$run_meta
    )

    y <- c(stats::rnorm(12L), stats::rchisq(12L, df = 2L))
    nmcd_args <- list(
      data = y,
      cpn_max = 1L,
      dm = 6L,
      cov_rate = 0.85,
      reg_fun = relieverChangepoint::reg_fun_nmcd,
      method = "SN",
      echo = FALSE
    )

    nmcd_cost_mat <- do.call(relieverChangepoint::reliever_generic, nmcd_args)
    nmcd_loss_block <- do.call(
      relieverChangepoint::reliever_generic,
      c(nmcd_args, list(cache_backend = "by_loss_block"))
    )

    expect_same_cpd_path(nmcd_loss_block, nmcd_cost_mat, tolerance = 1e-10)
  })
})
