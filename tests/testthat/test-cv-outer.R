test_that("compressed evaluation indices equal explicit original-time bounds", {
  n <- 17L
  fold_id <- rep(rep(1:3, each = 2L), length.out = n)

  for (fold in seq_len(3L)) {
    train_id <- which(fold_id != fold)
    eval_id <- which(fold_id == fold)
    eval_pos <- relieverChangepoint:::.cv_reliever_eval_index(
      train_id, eval_id
    )
    cpd_sets <- list(
      integer(),
      2L,
      c(2L, 5L),
      c(1L, length(train_id) - 1L)
    )

    for (cpd_train in cpd_sets) {
      original_right <- c(train_id[cpd_train], n)
      compressed_right <- c(cpd_train, length(train_id))
      original_segment <- vapply(
        eval_id,
        function(i) which(i <= original_right)[1L],
        integer(1L)
      )
      compressed_segment <- vapply(
        eval_pos,
        function(i) which(i <= compressed_right)[1L],
        integer(1L)
      )
      expect_equal(compressed_segment, original_segment)
    }
  }
})

test_that("outer CV constructs proportional grids from one full-data size", {
  expect_equal(
    relieverChangepoint:::.reliever_grid_from_size(
      n = 20L, dc_grid_size = 5L
    ),
    c(5L, 10L, 15L, 20L)
  )
  expect_equal(
    relieverChangepoint:::.reliever_grid_from_size(
      n = 22L, dc_grid_size = 5L
    ),
    c(5L, 10L, 15L, 20L, 22L)
  )
  expect_equal(
    relieverChangepoint:::.reliever_grid_from_size(
      n = 3L, dc_grid_size = 5L
    ),
    3L
  )
  expect_equal(
    relieverChangepoint:::.cv_reliever_fold_grid_size(
      dc_grid_size = 20L, n_train = 23L, n = 31L
    ),
    15L
  )
  expect_equal(
    relieverChangepoint:::.cv_reliever_fold_grid_size(
      dc_grid_size = 20L, n_train = 24L, n = 31L
    ),
    16L
  )
  expect_error(
    relieverChangepoint:::.reliever_grid_from_size(
      n = 20L, dc_grid_size = 0L
    ),
    "dc_grid_size must be a single integer >= 1"
  )
})

test_that("outer CV supports proportional grids on original-index output", {
  with_test_timeout({
    set.seed(20260719)
    n <- 62L
    x <- rbind(
      matrix(rnorm(20L * 4L, mean = 0, sd = 0.4), 20L, 4L),
      matrix(rnorm(20L * 4L, mean = 3, sd = 0.4), 20L, 4L),
      matrix(rnorm(22L * 4L, mean = -3, sd = 0.4), 22L, 4L)
    )
    expected_grid <- c(seq.int(10L, 60L, by = 10L), n)
    common <- list(
      cpn_max = 3L, dm = 8L, cov_rate = 0.6, method = "SN",
      nfolds = 3L, dc_grid_size = 10L
    )

    built_in <- do.call(
      relieverChangepoint::cv.reliever, c(list(X = x), common)
    )
    generic <- do.call(
      relieverChangepoint::cv.reliever_generic,
      c(
        list(data = x, reg_fun = relieverChangepoint::reg_fun_mean),
        common
      )
    )

    expect_equal(built_in$settings$dc_grid_size, 10L)
    expect_equal(built_in$settings$dc_grid, expected_grid)
    expect_equal(built_in$full_data_fit$settings$dc_grid_size, 10L)
    expect_equal(built_in$full_data_fit$settings$dc_grid, expected_grid)
    expect_equal(built_in$cv_loss$K, 0:2)
    expect_true(all(
      unlist(built_in$full_data_fit$cpd_path$candidates$cpd) %in%
        expected_grid
    ))
    expect_true(all(built_in$summary$cpd_hat[[1L]] %in% expected_grid))
    expect_equal(built_in$summary, generic$summary)
    expect_equal(built_in$cv_loss, generic$cv_loss, tolerance = 1e-9)
    expect_equal(
      built_in$full_data_fit$cpd_path$candidates,
      generic$full_data_fit$cpd_path$candidates,
      tolerance = 1e-9
    )
  })
})

test_that("outer CV aggregates non-K selectors even when fold K differs", {
  fold_loss <- data.frame(
    fold_id = 1:3,
    loss_output_id = 1L,
    selector_id = 2L,
    select_value = 5,
    K = 1:3,
    eval_loss = c(3, 4, 5),
    n_eval = 2L
  )
  out <- relieverChangepoint:::.cv_reliever_loss_table(fold_loss, 3L)

  expect_equal(nrow(out), 1L)
  expect_equal(out$select_value, 5)
  expect_equal(out$cv_loss, 12)
  expect_true(out$complete)
})

test_that("outer CV constructs order-preserved folds on the global time axis", {
  with_test_timeout({
    x <- matrix(seq_len(120), 30L, 4L)
    fit <- relieverChangepoint::cv.reliever_generic(
      x,
      reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 2L,
      dm = 5L,
      cov_rate = 0.6,
      nfolds = 3L,
      op_size = 2L,
      detail = TRUE
    )

    expect_equal(
      fit$diagnostics$fold_id,
      rep(rep(1:3, each = 2L), length.out = nrow(x))
    )
    expect_equal(fit$settings$op_size, 2L)
    expect_equal(fit$settings$fold_size, c(10L, 10L, 10L))
  })
})

test_that("outer CV matches adaptive hyperparameters by hyper_id", {
  with_test_timeout({
    dynamic_meta_loss <- function(data, l, r, l_end = l, r_end = r,
                                  save_model = FALSE,
                                  is_virtual_run = FALSE) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 1L,
          loss_output_meta = data.frame(
            loss_output_id = 1L,
            row_type = "custom",
            hyper_id = 1L,
            hyper_value = if (is.null(dim(data))) length(data) else nrow(data),
            auxiliary_n = if (is.null(dim(data))) length(data) else nrow(data)
          )
        ))
      }
      center <- mean(data[l:r])
      list(
        loss = matrix((data[l_end:r_end] - center)^2, ncol = 1L),
        model = NULL
      )
    }

    fit <- relieverChangepoint::cv.reliever_generic(
      data = seq_len(15L),
      reg_fun = dynamic_meta_loss,
      cpn_max = 1L,
      dm = 3L,
      cov_rate = 1,
      method = "SN",
      nfolds = 4L,
      cache_backend = "by_cost_mat"
    )

    expect_s3_class(fit, "cv_reliever_result")
    expect_identical(fit$full_data_fit$run_meta$hyper_id, 1L)
    expect_equal(fit$full_data_fit$run_meta$hyper_value, 15L)
    expect_equal(fit$full_data_fit$run_meta$auxiliary_n, 15L)
    expect_equal(fit$settings$fold_size, c(4L, 4L, 4L, 3L))
  })
})

test_that("reg_fun_mean returns the mean-square individual loss", {
  set.seed(2026)
  x <- matrix(rnorm(40), 20L, 2L)
  out <- relieverChangepoint::reg_fun_mean(
    x, l = 3L, r = 12L, l_end = 1L, r_end = 20L,
    save_model = TRUE
  )
  center <- colMeans(x[3:12, , drop = FALSE])
  expected <- rowMeans(sweep(x, 2L, center, "-")^2)

  expect_equal(as.numeric(out$loss), expected)
  expect_equal(out$model$center, center)
  expect_equal(
    relieverChangepoint::reg_fun_mean(
      x, 1L, 1L, is_virtual_run = TRUE
    )$loss_output_meta$row_type,
    "mean"
  )
})

test_that("outer CV selects two mean changepoints by held-out loss", {
  with_test_timeout({
    set.seed(2026)
    n_seg <- 30L
    p <- 20L
    x <- rbind(
      matrix(rnorm(n_seg * p, mean = 0, sd = 0.5), n_seg, p),
      matrix(rnorm(n_seg * p, mean = 4, sd = 0.5), n_seg, p),
      matrix(rnorm(n_seg * p, mean = -4, sd = 0.5), n_seg, p)
    )

    fit <- relieverChangepoint::cv.reliever_generic(
      x,
      reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 4L,
      dm = 8L,
      cov_rate = 0.6,
      method = "SN",
      nfolds = 3L
    )

    expect_s3_class(fit, "cv_reliever_result")
    printed <- capture.output(print(fit))
    expect_match(printed[1L], "Cross-validated Reliever")
    expect_match(printed[2L], "Method: SN \\| Folds: 3")
    expect_named(
      fit,
      c("summary", "cv_loss", "full_data_fit", "settings")
    )
    expect_equal(fit$summary$K_hat, 2L)
    expect_cpd_error_lte(fit$summary$cpd_hat[[1L]], c(30L, 60L))
    expect_equal(
      fit$full_data_fit$summary,
      relieverChangepoint:::.reliever_empty_summary()
    )
    expect_null(fit$diagnostics)
    expect_equal(fit$settings$selection, "outer_cv")
    expect_equal(fit$settings$method, "SN")
    expect_equal(fit$settings$cpn_max, 4L)
    expect_false("cpn_crit" %in% names(fit$settings))

    selected_path <- subset(fit$cv_loss, K == 2L)
    expect_equal(fit$summary$cv_mean, selected_path$cv_mean)
    expect_named(
      fit$cv_loss,
      c("K", "cv_mean", "cv_se")
    )

    native_cv <- relieverChangepoint::cv.reliever(
      x, cpn_max = 4L, dm = 8L, cov_rate = 0.6, method = "SN",
      nfolds = 3L
    )
    expect_equal(native_cv$summary$K_hat, 2L)
    expect_cpd_error_lte(native_cv$summary$cpd_hat[[1L]], c(30L, 60L))
    expect_equal(
      native_cv$cv_loss[c("K", "cv_mean", "cv_se")],
      fit$cv_loss[c("K", "cv_mean", "cv_se")],
      tolerance = 1e-9
    )

    fold_id <- rep(seq_len(3L), length.out = nrow(x))
    manual_fold_loss <- lapply(seq_len(3L), function(fold) {
      eval_id <- which(fold_id == fold)
      train_id <- which(fold_id != fold)
      train_fit <- relieverChangepoint::reliever_mean(
        x[train_id, , drop = FALSE],
        cpn_max = 4L, dm = 6L, cov_rate = 0.6, method = "SN",
        cpn_crit = "none"
      )
      evaluated <- relieverChangepoint::evaluate_reliever_segments(
        train_fit,
        data = x[train_id, , drop = FALSE],
        eval_data = x[eval_id, , drop = FALSE],
        eval_index = relieverChangepoint:::.cv_reliever_eval_index(
          train_id, eval_id
        )
      )
      data.frame(
        K = evaluated$candidates$K,
        eval_loss = evaluated$candidates$eval_loss
      )
    })
    manual_fold_loss <- do.call(rbind, manual_fold_loss)
    manual_cv_loss <- stats::aggregate(
      eval_loss ~ K, data = manual_fold_loss, FUN = sum
    )
    expect_equal(
      native_cv$cv_loss$cv_mean,
      manual_cv_loss$eval_loss[
        match(native_cv$cv_loss$K, manual_cv_loss$K)
      ] / nrow(x),
      tolerance = 1e-9
    )

    native_fit <- relieverChangepoint::reliever_mean(
      x, cpn_max = 4L, dm = 8L, cov_rate = 0.6, method = "SN",
      cpn_crit = "none"
    )
    expect_equal(
      fit$full_data_fit$cpd_path$candidates[c("K", "cpd")],
      native_fit$cpd_path$candidates[c("K", "cpd")]
    )
    expect_equal(
      fit$full_data_fit$cpd_path$candidates$loss,
      native_fit$cpd_path$candidates$loss,
      tolerance = 1e-9
    )
    expect_equal(
      native_cv$full_data_fit$cpd_path$candidates,
      native_fit$cpd_path$candidates
    )
  })
})

test_that("built-in NMCD and fixed KDE-L2 use outer CV to select K", {
  with_test_timeout({
    set.seed(2026)
    n_nmcd <- 90L
    x_nmcd <- c(
      rnorm(n_nmcd, mean = -4, sd = 0.35),
      rnorm(n_nmcd, mean = 4, sd = 0.35),
      rnorm(n_nmcd, mean = 0, sd = 0.35)
    )
    nmcd_args <- list(
      cpn_max = 4L, dm = 18L, cov_rate = 0.7,
      method = "SN", nfolds = 5L, dc_grid_size = 10L
    )
    nmcd <- do.call(
      relieverChangepoint::cv.reliever,
      c(list(X = x_nmcd, cpd_family = "nmcd"), nmcd_args)
    )
    nmcd_generic <- do.call(
      relieverChangepoint::cv.reliever_generic,
      c(
        list(
          data = x_nmcd,
          reg_fun = relieverChangepoint::reg_fun_nmcd
        ),
        nmcd_args
      )
    )

    expect_identical(nmcd$settings$cpd_family, "nmcd")
    expect_identical(nmcd$summary$K_hat, 2L)
    expect_cpd_error_lte(
      nmcd$summary$cpd_hat[[1L]], c(n_nmcd, 2L * n_nmcd)
    )
    expect_equal(nmcd$summary, nmcd_generic$summary)
    expect_equal(nmcd$cv_loss, nmcd_generic$cv_loss, tolerance = 1e-9)

    set.seed(2026)
    n_kde <- 50L
    x_kde <- c(
      rnorm(n_kde, mean = -3, sd = 0.45),
      rnorm(n_kde, mean = 3, sd = 0.45),
      rnorm(n_kde, mean = 0, sd = 0.45)
    )
    kde_args <- list(
      cpn_max = 4L, dm = 10L, cov_rate = 0.6,
      method = "SN", nfolds = 5L, dc_grid_size = 5L
    )
    kde_l2 <- do.call(
      relieverChangepoint::cv.reliever,
      c(
        list(
          X = x_kde,
          cpd_family = "kde_l2",
          kernel = "gaussian",
          bandwidth = 0.8
        ),
        kde_args
      )
    )
    kde_prepared <- relieverChangepoint:::.kernel_l2_prepare_input(
      x_kde, kernel = "gaussian", bandwidth = 0.8
    )
    kde_generic <- do.call(
      relieverChangepoint::cv.reliever_generic,
      c(
        list(
          data = kde_prepared$data,
          reg_fun = relieverChangepoint::reg_fun_kde_l2
        ),
        kde_args
      )
    )

    expect_identical(kde_l2$settings$cpd_family, "kde_l2")
    expect_identical(kde_l2$summary$K_hat, 2L)
    expect_cpd_error_lte(
      kde_l2$summary$cpd_hat[[1L]], c(n_kde, 2L * n_kde)
    )
    expect_equal(kde_l2$summary, kde_generic$summary)
    expect_equal(kde_l2$cv_loss, kde_generic$cv_loss, tolerance = 1e-9)
    expect_identical(kde_l2$settings$input_spec, kde_prepared$input_spec)
    expect_identical(
      kde_l2$full_data_fit$settings$input_spec,
      kde_prepared$input_spec
    )
  })
})

test_that("outer CV reports its rule and distinct loss-output types", {
  with_test_timeout({
    two_mean_losses <- function(data, l, r, l_end = l, r_end = r,
                                save_model = FALSE,
                                is_virtual_run = FALSE) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 2L,
          loss_output_meta = data.frame(
            loss_output_id = 1:2,
            row_type = c("primary", "secondary"),
            loss_kind = "rss"
          )
        ))
      }
      data <- as.matrix(data)
      center <- colMeans(data[l:r, , drop = FALSE])
      loss <- rowMeans(
        sweep(data[l_end:r_end, , drop = FALSE], 2L, center, "-")^2
      )
      list(loss = cbind(loss, loss + 0.25), model = NULL)
    }
    set.seed(20260724)
    data <- matrix(
      c(stats::rnorm(20L), stats::rnorm(20L, mean = 2)),
      ncol = 1L
    )
    fit <- relieverChangepoint::cv.reliever_generic(
      data,
      reg_fun = two_mean_losses,
      cpn_max = 1L,
      dm = 6L,
      cov_rate = 0.8,
      nfolds = 2L
    )

    expect_identical(fit$summary$rule, "outer_cv")
    expect_identical(fit$summary$row_type, "primary")
    expect_equal(
      unique(fit$cv_loss$row_type),
      c("primary", "secondary")
    )
    expect_false("loss_type" %in% names(fit$cv_loss))
  })
})

test_that("native and generic mean CV agree across K-indexed searches", {
  with_test_timeout({
    set.seed(2026)
    n_seg <- 20L
    p <- 4L
    x <- rbind(
      matrix(rnorm(n_seg * p, mean = 0, sd = 0.4), n_seg, p),
      matrix(rnorm(n_seg * p, mean = 3, sd = 0.4), n_seg, p),
      matrix(rnorm(n_seg * p, mean = -3, sd = 0.4), n_seg, p)
    )
    methods <- c("SN", "WBS", "WBS_recursive", "SeedBS", "BS")

    for (method in methods) {
      common <- list(
        cpn_max = 3L, dm = 6L, cov_rate = 0.6,
        method = method, nfolds = 3L, M = 40L, wbs_seed = 2026L
      )
      generic <- do.call(
        relieverChangepoint::cv.reliever_generic,
        c(
          list(data = x, reg_fun = relieverChangepoint::reg_fun_mean),
          common
        )
      )
      native <- do.call(
        relieverChangepoint::cv.reliever, c(list(X = x), common)
      )
      native_cost_mat <- do.call(
        relieverChangepoint::cv.reliever,
        c(list(X = x), common, list(cache_backend = "by_cost_mat"))
      )

      expect_equal(generic$summary$K_hat, 2L, info = method)
      expect_cpd_error_lte(
        generic$summary$cpd_hat[[1L]], c(20L, 40L), info = method
      )
      expect_equal(
        generic$cv_loss[c("K", "cv_mean", "cv_se")],
        native$cv_loss[c("K", "cv_mean", "cv_se")],
        tolerance = 1e-9,
        info = method
      )
      expect_equal(
        native$cv_loss,
        native_cost_mat$cv_loss,
        tolerance = 1e-9,
        info = method
      )
      expect_equal(
        generic$full_data_fit$cpd_path$candidates,
        native$full_data_fit$cpd_path$candidates,
        tolerance = 1e-9,
        info = method
      )
    }
  })
})

test_that("outer CV selects WBS stopping values through the common selector path", {
  with_test_timeout({
    set.seed(2026)
    n_seg <- 20L
    p <- 4L
    x <- rbind(
      matrix(rnorm(n_seg * p, mean = 0, sd = 0.4), n_seg, p),
      matrix(rnorm(n_seg * p, mean = 3, sd = 0.4), n_seg, p),
      matrix(rnorm(n_seg * p, mean = -3, sd = 0.4), n_seg, p)
    )
    common <- list(
      cpn_max = 4L, dm = 6L, cov_rate = 0.6,
      method = "WBS", nfolds = 3L, M = 40L, wbs_seed = 2026L,
      wbs_stop_crit = c(100, 10, 1)
    )
    generic <- do.call(
      relieverChangepoint::cv.reliever_generic,
      c(
        list(data = x, reg_fun = relieverChangepoint::reg_fun_mean),
        common
      )
    )
    native <- do.call(
      relieverChangepoint::cv.reliever, c(list(X = x), common)
    )
    native_cost_mat <- do.call(
      relieverChangepoint::cv.reliever,
      c(list(X = x), common, list(cache_backend = "by_cost_mat"))
    )

    expect_equal(generic$cv_loss$wbs_stop_crit, c(Inf, 100, 10, 1))
    expect_equal(generic$cv_loss$K, c(0L, 1L, 2L, 2L))
    expect_equal(generic$summary$wbs_stop_crit, 10)
    expect_equal(generic$summary$K_hat, 2L)
    expect_cpd_error_lte(generic$summary$cpd_hat[[1L]], c(20L, 40L))
    expect_named(
      summary(generic),
      c(
        "rule", "wbs_stop_crit", "K_hat", "cpd_hat", "cv_mean", "cv_se"
      )
    )
    expect_equal(generic$cv_loss, native$cv_loss, tolerance = 1e-9)
    expect_equal(native$cv_loss, native_cost_mat$cv_loss, tolerance = 1e-9)
  })
})

test_that("outer CV supports PELT and OP penalty paths in both cache backends", {
  with_test_timeout({
    set.seed(2026)
    n_seg <- 20L
    p <- 4L
    x <- rbind(
      matrix(rnorm(n_seg * p, mean = 0, sd = 0.4), n_seg, p),
      matrix(rnorm(n_seg * p, mean = 3, sd = 0.4), n_seg, p),
      matrix(rnorm(n_seg * p, mean = -3, sd = 0.4), n_seg, p)
    )

    for (method in c("PELT", "OP")) {
      common <- list(
        dm = 6L, cov_rate = 0.6, method = method,
        pen_val = c(0.5, 2, 8), nfolds = 3L
      )
      run_cv <- function(fun, args) do.call(fun, args)
      generic <- run_cv(
        relieverChangepoint::cv.reliever_generic,
        c(
          list(data = x, reg_fun = relieverChangepoint::reg_fun_mean),
          common
        )
      )
      generic_cost_mat <- run_cv(
        relieverChangepoint::cv.reliever_generic,
        c(list(
          data = x,
          reg_fun = relieverChangepoint::reg_fun_mean
        ), common, list(
          cache_backend = "by_cost_mat"
        ))
      )
      native <- run_cv(
        relieverChangepoint::cv.reliever,
        c(list(X = x), common)
      )
      native_cost_mat <- run_cv(
        relieverChangepoint::cv.reliever,
        c(list(X = x), common, list(cache_backend = "by_cost_mat"))
      )

      expect_equal(native$cv_loss$pen_val, c(Inf, 0.5, 2, 8),
                   info = method)
      expect_equal(native$summary$pen_val, 0.5, info = method)
      expect_equal(native$summary$K_hat, 2L, info = method)
      expect_cpd_error_lte(
        native$summary$cpd_hat[[1L]], c(20L, 40L), info = method
      )
      expect_equal(summary(native)$pen_val, 0.5, info = method)
      expect_equal(generic$cv_loss, native$cv_loss,
                   tolerance = 1e-9, info = method)
      expect_equal(generic$cv_loss, generic_cost_mat$cv_loss,
                   tolerance = 1e-9, info = method)
      expect_equal(native$cv_loss, native_cost_mat$cv_loss,
                   tolerance = 1e-9, info = method)
    }
  })
})

test_that("outer CV retains only K values available in every fold", {
  with_test_timeout({
    set.seed(2026)
    x <- matrix(rnorm(31L * 3L), 31L, 3L)
    fit <- relieverChangepoint::cv.reliever(
      x, cpn_max = 6L, dm = 7L, cov_rate = 0.6,
      method = "SN", nfolds = 4L, op_size = 2L,
      detail = TRUE
    )

    expect_equal(fit$settings$fold_size, c(8L, 8L, 8L, 7L))
    expect_equal(fit$settings$fold_dm, rep(6L, 4L))
    expect_equal(fit$cv_loss$K, 0:2)
    expect_equal(
      as.integer(tapply(
        fit$diagnostics$fold_loss$K,
        fit$diagnostics$fold_loss$fold_id,
        max
      )),
      c(2L, 2L, 2L, 3L)
    )
  })
})

test_that("outer CV accepts vector data with a zero-change path", {
  with_test_timeout({
    set.seed(2026)
    x <- rnorm(30L)
    generic <- relieverChangepoint::cv.reliever_generic(
      data = x,
      reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 0L, dm = 5L, cov_rate = 0.6, nfolds = 3L
    )
    native <- relieverChangepoint::cv.reliever(
      x, cpn_max = 0L, dm = 5L, cov_rate = 0.6, nfolds = 3L,
      owner_key = FALSE
    )

    expect_equal(generic$summary$K_hat, 0L)
    expect_length(generic$summary$cpd_hat[[1L]], 0L)
    expect_equal(generic$cv_loss$cv_mean, native$cv_loss$cv_mean)
  })
})

test_that("outer CV jointly selects a lasso loss output and K", {
  with_test_timeout({
    set.seed(2026)
    n <- 120L
    p <- 10L
    b0 <- c(3, -2.5, 2, rep(0, p - 3L))
    delta <- cbind(-2 * b0, 1.8 * b0)
    data <- relieverChangepoint::dgp_linear_regression(
      n, p, tau = c(40L, 80L), b0 = b0, delta = delta, sig = 1
    )$data
    lam_set <- c(20, 8, 3, 1, 0.3)

    fit <- relieverChangepoint::cv.reliever_generic(
      data,
      reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
      cpn_max = 4L,
      dm = 10L,
      cov_rate = 0.6,
      method = "SN",
      nfolds = 3L,
      run_cpd_ids = c(2L, 5L),
      lam_set = lam_set
    )

    expect_equal(fit$summary$K_hat, 2L)
    expect_cpd_error_lte(fit$summary$cpd_hat[[1L]], c(40L, 80L))
    selected_loss <- fit$cv_loss[which.min(fit$cv_loss$cv_mean), , drop = FALSE]
    expect_equal(fit$summary$cv_mean, selected_loss$cv_mean)
    expect_equal(fit$summary$hyper_value, selected_loss$hyper_value)
    expect_true(fit$summary$hyper_value %in% lam_set[c(2L, 5L)])
    expect_equal(attr(fit$summary, "hyper_name"), "lambda")
    expect_true(any(grepl(
      "lambda", capture.output(print(fit)), fixed = TRUE
    )))
    expect_equal(nrow(fit$cv_loss), 2L * 5L)
    expect_named(
      fit$cv_loss,
      c("hyper_value", "K", "cv_mean", "cv_se", "run_id")
    )

    selected_run <- fit$full_data_fit$run_meta$run_id[
      vapply(
        fit$full_data_fit$run_meta$hyper_value,
        identical,
        logical(1L),
        fit$summary$hyper_value
      )
    ]
    final_candidate <- subset(
      fit$full_data_fit$cpd_path$candidates,
      run_id == selected_run & K == fit$summary$K_hat
    )
    expect_equal(final_candidate$cpd[[1L]], fit$summary$cpd_hat[[1L]])
  })
})

test_that("cv.reliever dispatches mean and lasso outer CV", {
  with_test_timeout({
    set.seed(2026)
    n_seg <- 20L
    x_mean <- rbind(
      matrix(rnorm(n_seg * 4L, mean = 0, sd = 0.4), n_seg, 4L),
      matrix(rnorm(n_seg * 4L, mean = 3, sd = 0.4), n_seg, 4L),
      matrix(rnorm(n_seg * 4L, mean = -3, sd = 0.4), n_seg, 4L)
    )
    applied_mean <- relieverChangepoint::cv.reliever(
      x_mean, cpn_max = 3L, dm = 6L, cov_rate = 0.6,
      method = "SN", nfolds = 3L
    )
    expect_equal(applied_mean$settings$cpd_family, "mean")
    expect_equal(applied_mean$settings$ratio, 0.9)
    expect_named(
      summary(applied_mean),
      c("rule", "K_hat", "cpd_hat", "cv_mean", "cv_se")
    )
    expect_identical(summary(applied_mean), applied_mean$summary)
    expect_s3_class(applied_mean$summary, "reliever_summary")
    printed <- capture.output(print(applied_mean))
    expect_false(any(grepl("run_id|loss_output_id", printed)))
    expect_true(any(grepl(
      "Family: mean.*Method: SN.*Folds: 3", printed
    )))

    n <- 3L * n_seg
    p <- 4L
    x <- matrix(rnorm(n * p), nrow = n)
    beta <- rbind(
      matrix(c(1, 0, 0, 0), n_seg, p, byrow = TRUE),
      matrix(c(-1, 0.8, 0, 0), n_seg, p, byrow = TRUE),
      matrix(c(0.6, -0.6, 0, 0), n_seg, p, byrow = TRUE)
    )
    y <- rowSums(x * beta) + rnorm(n, sd = 0.2)
    data <- cbind(y, x)
    lam_set <- c(1, 0.2)
    applied_lasso <- NULL
    expect_warning(
      applied_lasso <- relieverChangepoint::cv.reliever(
        X = x, y = y, cpd_family = "lasso",
        cpn_max = 2L, dm = 7L, cov_rate = 0.6,
        method = "SN", nfolds = 2L, ratio = 0.8, lam_set = lam_set
      ),
      "ratio is ignored"
    )
    generic_lasso <- relieverChangepoint::cv.reliever_generic(
      data, reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
      cpn_max = 2L, dm = 7L, cov_rate = 0.6,
      method = "SN", nfolds = 2L, lam_set = lam_set
    )
    expect_equal(applied_lasso$summary, generic_lasso$summary)
    expect_equal(applied_lasso$cv_loss, generic_lasso$cv_loss)
    expect_equal(applied_lasso$settings$cpd_family, "lasso")
    expect_identical(
      applied_lasso$settings$input_spec,
      list(
        type = "response_predictor",
        n_predictors = p,
        original_form = "separate_xy"
      )
    )
    expect_identical(
      applied_lasso$full_data_fit$settings$input_spec,
      applied_lasso$settings$input_spec
    )
    expect_named(
      summary(applied_lasso),
      c(
        "rule", "hyper_value", "K_hat", "cpd_hat", "cv_mean", "cv_se"
      )
    )
    expect_equal(
      summary(applied_lasso)$hyper_value,
      unlist(applied_lasso$summary$hyper_value, use.names = FALSE)
    )

    expected_lasso_grid <- c(seq.int(7L, 56L, by = 7L), n)
    gridded_lasso <- relieverChangepoint::cv.reliever(
      X = x, y = y, cpd_family = "lasso",
      cpn_max = 0L, dm = 7L, cov_rate = 0.6,
      method = "SN", nfolds = 2L, dc_grid_size = 7L,
      lam_set = lam_set[1L]
    )
    gridded_generic_lasso <- relieverChangepoint::cv.reliever_generic(
      data = data,
      reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
      cpn_max = 0L, dm = 7L, cov_rate = 0.6,
      method = "SN", nfolds = 2L, dc_grid_size = 7L,
      lam_set = lam_set[1L]
    )
    expect_equal(
      gridded_lasso$full_data_fit$settings$dc_grid,
      expected_lasso_grid
    )
    expect_equal(gridded_lasso$summary, gridded_generic_lasso$summary)
    expect_equal(gridded_lasso$cv_loss, gridded_generic_lasso$cv_loss)

    auto_lam <- relieverChangepoint:::.reliever_auto_lam_set(data)
    automatic_lasso <- relieverChangepoint::cv.reliever(
      X = x, y = y, cpd_family = "lasso",
      cpn_max = 0L, dm = 7L, cov_rate = 0.6,
      method = "SN", nfolds = 2L
    )
    expect_equal(
      as.numeric(unlist(automatic_lasso$full_data_fit$run_meta$hyper_value)),
      auto_lam
    )
    short_auto_lam <- relieverChangepoint:::.reliever_auto_lam_set(
      data, nlambda = 8L
    )
    automatic_short_lasso <- relieverChangepoint::cv.reliever(
      X = x, y = y, cpd_family = "lasso",
      cpn_max = 0L, dm = 7L, cov_rate = 0.6,
      method = "SN", nfolds = 2L, nlambda = 8L
    )
    expect_equal(
      as.numeric(unlist(
        automatic_short_lasso$full_data_fit$run_meta$hyper_value
      )),
      short_auto_lam
    )
    expect_error(
      relieverChangepoint::cv.reliever(
        X = x_mean, cpn_max = 2L, dm = 6L, nfolds = 3L,
        dc_grid = seq(5L, nrow(x_mean), by = 5L)
      ),
      "Use dc_grid_size"
    )
    expect_error(
      relieverChangepoint::cv.reliever_generic(
        data = x_mean,
        reg_fun = relieverChangepoint::reg_fun_mean,
        cpn_max = 2L, dm = 6L, nfolds = 3L,
        dc_grid = seq(5L, nrow(x_mean), by = 5L)
      ),
      "Use dc_grid_size"
    )
    expect_error(
      relieverChangepoint::cv.reliever_generic(
        x_mean, reg_fun = relieverChangepoint::reg_fun_mean,
        cpn_max = 2L, dm = 6L, nfolds = 3L,
        cache_profile = list()
      ),
      "does not support cache_profile"
    )
    expect_error(
      relieverChangepoint::cv.reliever(
        x_mean, cpn_max = 2L, dm = 6L, nfolds = 3L,
        cpn_crit = "sic"
      ),
      "does not support cpn_crit"
    )
  })
})

test_that("outer CV validates methods and fold construction", {
  x <- matrix(seq_len(60), 30L, 2L)
  expect_error(
    relieverChangepoint::cv.reliever_generic(
      data = x,
      reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 2L, dm = 5L, method = "invalid", nfolds = 3L
    ),
    "one of"
  )
  expect_error(
    relieverChangepoint::cv.reliever_generic(
      data = x,
      reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 2L, dm = 5L, nfolds = 1L
    ),
    "between 2"
  )
  expect_error(
    relieverChangepoint::cv.reliever_generic(
      data = x,
      reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 2L, dm = 5L, nfolds = 5L, op_size = 8L
    ),
    "every outer fold nonempty"
  )
  expect_error(
    relieverChangepoint::cv.reliever(
      X = x, cpn_max = 2L, dm = 5L, nfolds = 3L,
      dc_grid_size = 1.5
    ),
    "dc_grid_size must be a single integer >= 1"
  )
  expect_warning(
    relieverChangepoint::cv.reliever(
      X = x, cpn_max = 0L, dm = 5L, cov_rate = 0.6,
      nfolds = 3L, ratio = 0.8, dc_grid_size = 5L
    ),
    "ratio is ignored"
  )
})

test_that("outer CV reports a full-search cache switch only once", {
  with_test_timeout({
    x <- matrix(seq_len(60), 30L, 2L)
    messages <- character()
    fit <- withCallingHandlers(
      relieverChangepoint::cv.reliever(
        x,
        cpn_max = 1L,
        dm = 5L,
        cov_rate = 1,
        nfolds = 2L,
        cache_backend = "by_loss_block"
      ),
      reliever_cache_backend_warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    expect_s3_class(fit, "cv_reliever_result")
    expect_length(messages, 1L)
    expect_match(messages, "using cache_backend = \"by_cost_mat\"", fixed = TRUE)
  })
})

test_that("outer CV echo describes folds and the final all-observation fit", {
  messages <- capture.output(
    fit <- relieverChangepoint::cv.reliever(
      X = seq_len(20L),
      cpn_max = 0L,
      dm = 5L,
      nfolds = 2L,
      echo = TRUE
    ),
    type = "message"
  )

  expect_s3_class(fit, "cv_reliever_result")
  expect_identical(
    messages,
    c(
      "Outer CV: fold 1/2",
      "Outer CV: fold 2/2",
      "Outer CV: fitting changepoint candidates on all observations"
    )
  )
})
