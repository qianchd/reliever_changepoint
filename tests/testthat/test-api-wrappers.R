.wrapper_lasso_data <- function(n_seg = 12L, p = 4L) {
  n <- 3L * n_seg
  x <- matrix(stats::rnorm(n * p), nrow = n)
  beta <- rbind(
    matrix(c(1, 0, 0, 0), n_seg, p, byrow = TRUE),
    matrix(c(-1, 0.8, 0, 0), n_seg, p, byrow = TRUE),
    matrix(c(0.6, -0.6, 0, 0), n_seg, p, byrow = TRUE)
  )
  y <- rowSums(x * beta) + stats::rnorm(n, sd = 0.2)
  cbind(y, x)
}

.wrapper_mean_data <- function(n_seg = 12L, p = 3L) {
  rbind(
    matrix(stats::rnorm(n_seg * p, mean = 0), n_seg, p),
    matrix(stats::rnorm(n_seg * p, mean = 2), n_seg, p),
    matrix(stats::rnorm(n_seg * p, mean = -2), n_seg, p)
  )
}

.expect_same_wrapper_result <- function(wrapper_res, direct_res,
                                        tolerance = 1e-8) {
  expect_same_cpd_path(wrapper_res, direct_res, tolerance = tolerance)
  expect_equal(wrapper_res$run_meta, direct_res$run_meta)
  expect_equal(wrapper_res$summary, direct_res$summary)
  expect_equal(wrapper_res$settings$search_n, direct_res$settings$search_n)
  expect_equal(wrapper_res$settings$cache_backend, direct_res$settings$cache_backend)
}

test_that("released reliever wrappers are exported", {
  with_test_timeout({
    exports <- getNamespaceExports("relieverChangepoint")
    expect_true("reliever" %in% exports)
    expect_true("cv.reliever" %in% exports)
    expect_true("cv.reliever_generic" %in% exports)
    expect_false("cv.reliever_mean" %in% exports)
    expect_true("reliever_generic" %in% exports)
    expect_true(all(
      c(
        "reliever_var", "reliever_meanvar", "reliever_lm",
        "reliever_glm", "reliever_em"
      ) %in% exports
    ))
    expect_true("reliever_lasso" %in% exports)
    expect_true("reliever_lasso_crossfit" %in% exports)
    expect_true("reliever_mean_crossfit" %in% exports)
    expect_true("reliever_kde_nll" %in% exports)
    expect_true("reliever_kde_nll_crossfit" %in% exports)
    expect_true("reliever_kde_l2" %in% exports)
    expect_true("reliever_nmcd" %in% exports)
    expect_true("reliever_ranger_crossfit" %in% exports)
    expect_true("reliever_mlp_crossfit" %in% exports)
    expect_true("select_by_run" %in% exports)
    expect_true("select_across_runs" %in% exports)
    expect_false("reselect_by_run" %in% exports)
    expect_false("reselect_across_runs" %in% exports)
    expect_true("select_holdout" %in% exports)
    expect_false("reselect_holdout" %in% exports)
    expect_false("compare_cpn_criteria" %in% exports)
    expect_false(any(grepl("_recv$", exports)))
    expect_false(any(grepl("^model_select_", exports)))
    expect_false(".reliever_select_fixed_crossfit" %in% exports)
  })
})

test_that("built-in entry points identify invalid cpd families", {
  with_test_timeout({
    x <- matrix(seq_len(20L), ncol = 2L)
    expect_identical(
      relieverChangepoint:::.reliever_match_cpd_family(
        "mea", c("mean", "lasso"), "cv.reliever"
      ),
      "mean"
    )
    expect_error(
      relieverChangepoint::reliever(x, cpd_family = "unknown"),
      "cpd_family must be one of .* in reliever\\(\\)"
    )
    expect_error(
      relieverChangepoint::cv.reliever(x, cpd_family = "unknown"),
      paste0(
        'cpd_family must be one of "mean", "var", "meanvar", "lm", ',
        '"glm", "em", "lasso", "kde_l2", "nmcd" in ',
        'cv.reliever\\(\\)'
      )
    )
  })
})

test_that("applied and developer entry points expose distinct API levels", {
  with_test_timeout({
    expect_identical(
      names(formals(relieverChangepoint::reliever)),
      c("X", "y", "cpd_family", "cpn_max", "dm", "cov_rate", "method",
        "cpn_crit", "pen_val", "prune_value", "M", "wbs_seed",
        "wbs_stop_crit", "detail", "cache_backend", "owner_key", "echo",
        "dc_grid_size", "dc_grid", "...")
    )
    generic_formals <- formals(relieverChangepoint::reliever_generic)
    expect_identical(
      names(generic_formals),
      c(
        "data", "reg_fun", "cpn_max", "dm", "cov_rate", "method",
        "cpn_crit", "pen_val", "prune_value", "M", "wbs_seed",
        "wbs_stop_crit", "detail", "cache_backend", "owner_key", "echo",
        "dc_grid_size", "dc_grid", "cache_profile", "run_cpd_ids", "..."
      )
    )
    expect_identical(generic_formals$reg_fun, quote(expr = ))
    expect_identical(generic_formals$cpn_crit, "none")
    expect_identical(generic_formals$dc_grid_size, NULL)
    expect_identical(
      formals(relieverChangepoint::reliever)$cpn_crit,
      "none"
    )
    expect_error(
      relieverChangepoint::reliever(
        matrix(seq_len(12L), ncol = 1L), L = 1L
      ),
      "argument L was renamed to cpn_max"
    )
    expect_error(
      relieverChangepoint::reliever_generic(
        matrix(seq_len(12L), ncol = 1L),
        reg_fun = relieverChangepoint::reg_fun_mean,
        L = 1L
      ),
      "argument L was renamed to cpn_max"
    )
    expect_error(
      relieverChangepoint::reliever_generic(matrix(1:12, ncol = 1L)),
      "reg_fun"
    )
    generic_fit <- relieverChangepoint::reliever_generic(
      matrix(seq_len(18L), ncol = 1L),
      reg_fun = reg_null,
      cpn_max = 1L, dm = 4L, cov_rate = 0.8
    )
    expect_identical(generic_fit$settings$cpn_crit, "none")
    expect_equal(nrow(summary(generic_fit)), 0L)
    uppercase_data <- matrix(seq_len(18L), ncol = 1L)
    uppercase_none <- list(
      relieverChangepoint::reliever(
        uppercase_data, cpn_max = 1L, dm = 4L, cov_rate = 0.8,
        cpn_crit = "NONE"
      ),
      relieverChangepoint::reliever_mean(
        uppercase_data, cpn_max = 1L, dm = 4L, cov_rate = 0.8,
        cpn_crit = "NONE"
      ),
      relieverChangepoint::reliever_generic(
        uppercase_data, reg_fun = reg_null,
        cpn_max = 1L, dm = 4L, cov_rate = 0.8,
        cpn_crit = "NONE"
      )
    )
    expect_true(all(vapply(
      uppercase_none,
      function(fit) {
        identical(fit$settings$cpn_crit, "none") &&
          nrow(summary(fit)) == 0L
      },
      logical(1L)
    )))
    expect_setequal(
      names(formals(relieverChangepoint::select_holdout)),
      c(
        "result", "data", "eval_data", "eval_index", "y", "eval_y", "K",
        "run_ids", "reg_fun", "data_stack_fun", "...", "run_type"
      )
    )
    expect_identical(
      formals(relieverChangepoint::select_by_run)$cpn_crit,
      quote(expr = )
    )
    expect_identical(
      formals(relieverChangepoint::select_across_runs)$cpn_crit,
      quote(expr = )
    )
    expect_setequal(
      names(formals(relieverChangepoint::select_by_run)),
      c("result", "run_ids", "run_type", "cpn_crit", "n")
    )
    expect_setequal(
      names(formals(relieverChangepoint::select_across_runs)),
      c("result", "run_ids", "run_type", "cpn_crit", "n")
    )
    expect_true(all(
      c(
        "result", "data", "eval_data", "eval_index",
        "y", "eval_y", "save_model", "run_ids"
      ) %in% names(formals(
        relieverChangepoint::evaluate_reliever_segments
      ))
    ))
    expect_true(
      "run_type" %in% names(formals(
        relieverChangepoint::evaluate_reliever_segments
      ))
    )
    expect_true(
      "run_type" %in% names(formals(getS3method("plot", "reliever_result")))
    )
    expect_true(
      "run_type" %in% names(formals(getS3method("plot", "cv_reliever_result")))
    )

    focused_wrappers <- list(
      relieverChangepoint::reliever_mean,
      relieverChangepoint::reliever_var,
      relieverChangepoint::reliever_meanvar,
      relieverChangepoint::reliever_lm,
      relieverChangepoint::reliever_glm,
      relieverChangepoint::reliever_em,
      relieverChangepoint::reliever_lasso,
      relieverChangepoint::reliever_lasso_crossfit,
      relieverChangepoint::reliever_mean_crossfit,
      relieverChangepoint::reliever_kde_nll,
      relieverChangepoint::reliever_kde_nll_crossfit,
      relieverChangepoint::reliever_kde_l2,
      relieverChangepoint::reliever_nmcd,
      relieverChangepoint::reliever_ranger_crossfit,
      relieverChangepoint::reliever_mlp_crossfit
    )
    common_arguments <- c(
      "data", "cpn_max", "dm", "cov_rate", "method", "cpn_crit",
      "pen_val", "prune_value", "M", "wbs_seed", "wbs_stop_crit",
      "detail", "cache_backend", "owner_key", "echo", "dc_grid_size"
    )
    expect_true(all(vapply(
      focused_wrappers,
      function(fun) {
        all(common_arguments %in% names(formals(fun)))
      },
      logical(1L)
    )))
    expect_true(all(vapply(
      focused_wrappers,
      function(fun) identical(formals(fun)$method, "SN") &&
        identical(formals(fun)$cpn_crit, "none") &&
        identical(formals(fun)$dc_grid_size, NULL) &&
        identical(formals(fun)$dc_grid, NULL) &&
        identical(formals(fun)$cache_backend, "by_loss_block") &&
        identical(formals(fun)$owner_key, TRUE),
      logical(1L)
    )))
    expect_identical(formals(relieverChangepoint::reliever_mean)$owner_key,
                     TRUE)
    expect_identical(formals(relieverChangepoint::reliever_mean)$cpn_crit,
                     "none")
    expect_identical(formals(relieverChangepoint::cv.reliever)$owner_key,
                     TRUE)
    expect_identical(formals(relieverChangepoint::reliever_lasso_crossfit)$fold_type,
                     "op")
    expect_identical(formals(relieverChangepoint::reg_fun_clf_crossfit_template)$fold_type,
                     "op")
    expect_identical(formals(relieverChangepoint::reg_fun_ranger_crossfit)$fold_type,
                     "op")
    expect_identical(formals(relieverChangepoint::reg_fun_mlp_crossfit)$fold_type,
                     "op")
    expect_true(all(
      c(
        "nfolds", "bandwidth_vec", "var_dim", "fold_type",
        "op_size", "buffer_lag", "fold_stable_const", "kernel"
      ) %in% names(formals(
        relieverChangepoint::reg_fun_kde_nll_crossfit
      ))
    ))
  })
})

test_that("focused wrappers place shared arguments before model arguments", {
  common <- c(
    "data", "cpn_max", "dm", "cov_rate", "method", "cpn_crit",
    "pen_val", "prune_value", "M", "wbs_seed", "wbs_stop_crit",
    "detail", "cache_backend", "owner_key", "echo", "dc_grid_size",
    "cache_profile", "run_cpd_ids", "dc_grid"
  )
  crossfit <- c(
    "nfolds", "fold_type", "op_size", "buffer_lag", "fold_stable_const",
    "loss_output_types"
  )
  expected <- list(
    reliever_var = c(common, "mu"),
    reliever_meanvar = common,
    reliever_lm = c(common, "intercept"),
    reliever_glm = c(common, "family", "intercept", "response_ncol"),
    reliever_em = c(common, "family", "size"),
    reliever_lasso = c(common, "lam_set", "nlambda", "family", "thresh"),
    reliever_lasso_crossfit = c(
      common, crossfit, "lam_set", "nlambda", "family", "thresh"
    ),
    reliever_mean_crossfit = c(common, crossfit),
    reliever_kde_nll = c(
      common, "bandwidth_vec", "var_dim", "kernel", "kernel_args", "input_type"
    ),
    reliever_kde_nll_crossfit = c(
      common, crossfit,
      "bandwidth_vec", "var_dim", "kernel", "kernel_args", "input_type"
    ),
    reliever_kde_l2 = c(common, "kernel", "bandwidth", "kernel_args"),
    reliever_nmcd = c(common, "w_trunc", "sort_X"),
    reliever_ranger_crossfit = c(
      common, crossfit, "hyper_set", "ranger_args"
    ),
    reliever_mlp_crossfit = c(
      common, crossfit, "hyper_set", "nnet_args"
    )
  )
  for (name in names(expected)) {
    expect_identical(
      names(formals(getExportedValue("relieverChangepoint", name))),
      expected[[name]],
      info = name
    )
  }

  native_mean_common <- common[!common %in% c(
    "cache_profile", "run_cpd_ids", "dc_grid"
  )]
  expect_identical(
    names(formals(relieverChangepoint::reliever_mean)),
    c(native_mean_common, "dc_grid", "ratio")
  )
})

test_that("crossfit implementations place template arguments first", {
  inherited <- c(
    "data", "l", "r", "l_end", "r_end",
    "nfolds", "fold_type", "op_size", "buffer_lag", "fold_stable_const",
    "loss_output_types", "save_model", "is_virtual_run"
  )
  expected <- list(
    reg_fun_mean_crossfit = inherited,
    reg_fun_kde_nll_crossfit = c(
      inherited,
      "bandwidth_vec", "var_dim", "kernel", "kernel_args", "distance_power"
    ),
    reg_fun_lasso_crossfit = c(
      inherited, "lam_set", "family", "thresh"
    ),
    reg_fun_ranger_crossfit = c(inherited, "hyper_set", "ranger_args"),
    reg_fun_mlp_crossfit = c(inherited, "hyper_set", "nnet_args")
  )
  for (name in names(expected)) {
    expect_identical(
      names(formals(getExportedValue("relieverChangepoint", name))),
      expected[[name]],
      info = name
    )
  }
})

test_that("reliever dispatches built-in families and handles X/y inputs", {
  with_test_timeout({
    set.seed(20260712)
    data <- .wrapper_lasso_data()
    lam_set <- c(0.5, 0.1)

    direct_lasso <- relieverChangepoint::reliever_lasso(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, echo = FALSE
    )
    dispatched_lasso <- relieverChangepoint::reliever(
      X = data[, -1, drop = FALSE],
      y = data[, 1],
      cpd_family = "lasso",
      cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, echo = FALSE
    )
    .expect_same_wrapper_result(dispatched_lasso, direct_lasso)
    expect_equal(dispatched_lasso$settings$cpd_family, "lasso")
    expect_false("cpd_family" %in% names(dispatched_lasso))
    expect_equal(dispatched_lasso$settings$cpn_crit, "none")

    dispatched_lasso_no_owner_key <- relieverChangepoint::reliever(
      X = data[, -1, drop = FALSE],
      y = data[, 1],
      cpd_family = "lasso",
      cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, owner_key = FALSE, echo = FALSE
    )
    expect_same_cpd_path(dispatched_lasso_no_owner_key, direct_lasso)

    fallback_lasso <- NULL
    expect_warning(
      fallback_lasso <- relieverChangepoint::reliever(
        X = data,
        cpd_family = "lasso",
        cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
        lam_set = lam_set, echo = FALSE
      ),
      "Pass y explicitly"
    )
    .expect_same_wrapper_result(fallback_lasso, direct_lasso)

    direct_lasso_crossfit <- relieverChangepoint::reliever_lasso_crossfit(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    dispatched_lasso_crossfit <- relieverChangepoint::reliever(
      X = data[, -1, drop = FALSE],
      y = data[, 1],
      cpd_family = "lasso_crossfit",
      cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    .expect_same_wrapper_result(dispatched_lasso_crossfit, direct_lasso_crossfit)
    expect_equal(dispatched_lasso_crossfit$settings$cpd_family, "lasso_crossfit")
    expect_equal(dispatched_lasso_crossfit$settings$cpn_crit, "none")
    expect_equal(nrow(summary(dispatched_lasso_crossfit)), 0L)
    expect_false("selection" %in% names(dispatched_lasso_crossfit$settings))

    mean_data <- .wrapper_mean_data()
    dispatched_mean <- relieverChangepoint::reliever(
      X = mean_data,
      cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN", echo = FALSE
    )
    direct_mean <- relieverChangepoint::reliever_mean(
      mean_data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN", echo = FALSE
    )
    .expect_same_wrapper_result(dispatched_mean, direct_mean)
    expect_equal(dispatched_mean$settings$cpd_family, "mean")
    expect_equal(dispatched_mean$settings$cpn_crit, "none")
    expect_equal(nrow(summary(dispatched_mean)), 0L)

    direct_mean_wbs <- relieverChangepoint::reliever_mean(
      mean_data, cpn_max = 2L, dm = 8L, cov_rate = 0.8,
      method = "WBS", M = 8L, wbs_seed = 2026L, echo = FALSE
    )
    dispatched_mean_wbs <- relieverChangepoint::reliever(
      mean_data, cpn_max = 2L, dm = 8L, cov_rate = 0.8,
      method = "WBS", M = 8L, wbs_seed = 2026L, echo = FALSE
    )
    .expect_same_wrapper_result(dispatched_mean_wbs, direct_mean_wbs)
    expect_identical(dispatched_mean_wbs$settings$wbs_seed, 2026L)

    mean_vector <- c(rep(0, 12L), rep(4, 12L), rep(-4, 12L)) +
      stats::rnorm(36L, sd = 0.1)
    direct_mean_crossfit <- relieverChangepoint::reliever_mean_crossfit(
      mean_vector, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    dispatched_mean_crossfit <- relieverChangepoint::reliever(
      mean_vector, cpd_family = "mean_crossfit",
      cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    .expect_same_wrapper_result(dispatched_mean_crossfit, direct_mean_crossfit)
    selected_mean_crossfit <- relieverChangepoint::select_by_run(
      dispatched_mean_crossfit, run_type = "recv", cpn_crit = "loss"
    )
    expect_equal(selected_mean_crossfit$K_hat, 2L)
    expect_cpd_error_lte(
      selected_mean_crossfit$cpd_hat[[1L]],
      c(12L, 24L),
      tolerance = 3L
    )

    expect_warning(
      relieverChangepoint::reliever(
        X = mean_data, y = seq_len(nrow(mean_data)), cpd_family = "mean",
        cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN", echo = FALSE
      ),
      "Omit y"
    )
  })
})

test_that("lasso wrappers construct one reusable data-adaptive lambda path", {
  with_test_timeout({
    set.seed(20260712)
    data <- .wrapper_lasso_data(n_seg = 10L)
    auto_lam <- relieverChangepoint:::.reliever_auto_lam_set(data)
    base_model <- glmnet::glmnet(
      data[, -1, drop = FALSE], data[, 1],
      nlambda = 30L, intercept = FALSE, standardize = FALSE,
      control = list(thresh = 1e-7)
    )
    base_lam <- sort(unique(base_model$lambda * sqrt(nrow(data))),
                     decreasing = TRUE)

    expect_true(all(is.finite(auto_lam) & auto_lam > 0))
    expect_true(all(diff(auto_lam) < 0))
    expect_equal(length(auto_lam), length(base_lam) + 4L)
    expect_gt(max(auto_lam), max(base_lam))
    expect_lt(min(auto_lam), min(base_lam))

    short_nlambda <- 8L
    short_auto_lam <- relieverChangepoint:::.reliever_auto_lam_set(
      data, nlambda = short_nlambda
    )
    short_base_model <- glmnet::glmnet(
      data[, -1, drop = FALSE], data[, 1],
      nlambda = short_nlambda, intercept = FALSE, standardize = FALSE,
      control = list(thresh = 1e-7)
    )
    short_base_lam <- sort(
      unique(short_base_model$lambda * sqrt(nrow(data))),
      decreasing = TRUE
    )
    expect_equal(length(short_auto_lam), length(short_base_lam) + 4L)
    expect_lt(length(short_auto_lam), length(auto_lam))

    automatic <- relieverChangepoint::reliever_lasso(
      data, cpn_max = 1L, dm = 7L, cov_rate = 0.6, method = "SN",
      echo = FALSE
    )
    explicit <- relieverChangepoint::reliever_lasso(
      data, cpn_max = 1L, dm = 7L, cov_rate = 0.6, method = "SN",
      lam_set = auto_lam, echo = FALSE
    )
    .expect_same_wrapper_result(automatic, explicit)
    expect_equal(
      as.numeric(unlist(automatic$run_meta$hyper_value)),
      auto_lam
    )
    expect_equal(automatic$settings$loss_args$lam_set, auto_lam)
    expect_identical(automatic$settings$loss_args$family, "gaussian")
    expect_identical(automatic$settings$loss_args$thresh, 1e-7)

    automatic_short <- relieverChangepoint::reliever_lasso(
      data, cpn_max = 0L, dm = 7L, cov_rate = 0.6, method = "SN",
      nlambda = short_nlambda, echo = FALSE
    )
    explicit_short <- relieverChangepoint::reliever_lasso(
      data, cpn_max = 0L, dm = 7L, cov_rate = 0.6, method = "SN",
      lam_set = short_auto_lam, nlambda = 1L, echo = FALSE
    )
    .expect_same_wrapper_result(automatic_short, explicit_short)
    dispatched_short <- relieverChangepoint::reliever(
      X = data[, -1L, drop = FALSE], y = data[, 1L],
      cpd_family = "lasso",
      cpn_max = 0L, dm = 7L, cov_rate = 0.6, method = "SN",
      nlambda = short_nlambda, echo = FALSE
    )
    .expect_same_wrapper_result(automatic_short, dispatched_short)
    expect_equal(
      automatic_short$settings$loss_args$lam_set,
      short_auto_lam
    )

    automatic_crossfit_short <-
      relieverChangepoint::reliever_lasso_crossfit(
        data, cpn_max = 0L, dm = 7L, cov_rate = 0.6, method = "SN",
        nfolds = 2L, nlambda = short_nlambda, echo = FALSE
      )
    expect_equal(
      automatic_crossfit_short$settings$loss_args$lam_set,
      short_auto_lam
    )
    expect_error(
      relieverChangepoint:::.reliever_auto_lam_set(data, nlambda = 1L),
      "nlambda must be a single integer >= 2"
    )
    expect_named(summary(automatic), c("rule", "K_hat", "cpd_hat"))
    expect_equal(nrow(summary(automatic)), 0L)
    selected_sic <- relieverChangepoint::select_by_run(
      automatic, cpn_crit = "rss_sic"
    )
    expect_equal(selected_sic$hyper_value, auto_lam)

    holdout <- data
    holdout[, 1L] <- holdout[, 1L] + stats::rnorm(nrow(data), sd = 0.05)
    selected_holdout <- relieverChangepoint::select_holdout(
      automatic,
      data = data,
      eval_data = holdout,
      eval_index = seq_len(nrow(data)),
      K = 1L
    )
    expect_equal(selected_holdout$K_hat, 1L)
    expect_true(selected_holdout$hyper_value %in% auto_lam)
    expect_true(is.finite(selected_holdout$score))

    expect_error(
      relieverChangepoint::reg_fun_lasso_solpath(
        data, 1L, 10L, is_virtual_run = TRUE
      ),
      "lam_set is required"
    )
    expect_error(
      relieverChangepoint::reg_fun_lasso_crossfit(
        data, 1L, 10L, is_virtual_run = TRUE
      ),
      "lam_set is required"
    )
  })
})

test_that("summary methods expose the compact selected result", {
  with_test_timeout({
    set.seed(20260712)
    data <- .wrapper_mean_data()
    fit <- relieverChangepoint::reliever(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.6, method = "SN",
      cpn_crit = "rss_sic"
    )
    expect_named(summary(fit), c("rule", "K_hat", "cpd_hat"))
    expect_identical(summary(fit), fit$summary)
    printed <- capture.output(print(fit))
    expect_false(any(grepl("run_id|candidate_id", printed)))
    expect_true(any(grepl("Family: mean.*Method: SN", printed)))

    path_fit <- relieverChangepoint::reliever(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.6, method = "SN"
    )
    empty_print <- paste(
      capture.output(print(summary(path_fit))), collapse = "\n"
    )
    expect_match(empty_print, "fits candidate changepoint models")
    expect_match(empty_print, "cv.reliever\\(\\)")
    expect_match(empty_print, "select_by_run", fixed = TRUE)
    expect_match(empty_print, 'run_type = "recv"', fixed = TRUE)

    fit$summary$cpd_hat[[1L]] <- c(5L, 10L, 15L, 20L)
    full_cpd_print <- capture.output(print(fit))
    expect_true(any(grepl("5, 10, 15, 20", full_cpd_print, fixed = TRUE)))
    summary_print <- capture.output(print(summary(fit)))
    expect_true(any(grepl("5, 10, 15, 20", summary_print, fixed = TRUE)))
  })
})

test_that("lasso reliever wrappers match direct reg_fun calls on two changes", {
  with_test_timeout({
    set.seed(20260712)
    data <- .wrapper_lasso_data()
    lam_set <- c(0.5, 0.1)

    wrapper_path <- relieverChangepoint::reliever_lasso(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, cpn_crit = "rss_sic", echo = FALSE
    )
    direct_path <- relieverChangepoint::reliever_generic(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_lasso_solpath, lam_set = lam_set,
      cpn_crit = "rss_sic", echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_path, direct_path)
    expect_equal(wrapper_path$settings$cpn_crit, "rss_sic")
    selection_note <- "K was selected separately within each lambda value"
    expect_equal(sum(grepl(
      selection_note, capture.output(print(wrapper_path)), fixed = TRUE
    )), 1L)
    expect_equal(sum(grepl(
      selection_note, capture.output(print(summary(wrapper_path))),
      fixed = TRUE
    )), 1L)
    expect_true(any(grepl(
      "lambda", capture.output(print(summary(wrapper_path))), fixed = TRUE
    )))

    wrapper_crossfit <- relieverChangepoint::reliever_lasso_crossfit(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, nfolds = 2L, fold_type = "blk",
      loss_output_types = c("recv", "crossfit_homo_hyper"),
      cpn_crit = "loss", echo = FALSE
    )
    direct_crossfit <- relieverChangepoint::reliever_generic(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_lasso_crossfit, lam_set = lam_set,
      nfolds = 2L, fold_type = "blk",
      loss_output_types = c("recv", "crossfit_homo_hyper"),
      cpn_crit = "loss", echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_crossfit, direct_crossfit)
    expect_equal(summary(wrapper_crossfit)$rule, "loss")
    expect_false("hyper_value" %in% names(summary(wrapper_crossfit)))
    expect_equal(wrapper_crossfit$settings$cpn_crit, "loss")
    expect_false("selection" %in% names(wrapper_crossfit$settings))
    default_meta <- wrapper_crossfit$run_meta[
      wrapper_crossfit$run_meta$default_selection, , drop = FALSE
    ]
    expect_equal(nrow(default_meta), 1L)
    expect_equal(default_meta$row_type, "recv")

    fixed_cf_runs <- with(
      wrapper_crossfit$run_meta,
      run_id[row_type == "crossfit_homo_hyper"]
    )
    homogeneous_tuning <- relieverChangepoint::select_across_runs(
      wrapper_crossfit,
      run_type = "crossfit_homo_hyper",
      cpn_crit = "loss"
    )
    expect_true(homogeneous_tuning$run_id %in% fixed_cf_runs)
    expect_true("hyper_value" %in% names(homogeneous_tuning))
    expect_equal(nrow(homogeneous_tuning), 1L)
    expect_false(any(grepl(
      selection_note, capture.output(print(summary(wrapper_crossfit))),
      fixed = TRUE
    )))
    expect_equal(
      nrow(wrapper_crossfit$run_meta),
      1L + length(lam_set)
    )

    no_selection <- relieverChangepoint::reliever_lasso_crossfit(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      lam_set = lam_set, nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    expect_equal(nrow(summary(no_selection)), 0L)
    expect_false("selection" %in% names(no_selection$settings))
  })
})

test_that("mean and KDE-L2 reliever wrappers match direct calls on two changes", {
  with_test_timeout({
    set.seed(20260712)
    data <- .wrapper_mean_data()

    wrapper_mean_crossfit <- relieverChangepoint::reliever_mean_crossfit(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    direct_mean_crossfit <- relieverChangepoint::reliever_generic(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_mean_crossfit,
      nfolds = 2L, fold_type = "blk", cpn_crit = "none", echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_mean_crossfit, direct_mean_crossfit)
    expect_equal(wrapper_mean_crossfit$settings$cpn_crit, "none")
    expect_identical(wrapper_mean_crossfit$settings$loss_args$nfolds, 2L)
    expect_identical(wrapper_mean_crossfit$settings$loss_args$fold_type, "blk")

    mean_crossfit_unpenalized <- relieverChangepoint::reliever_mean_crossfit(
      data, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "blk", cpn_crit = "loss", echo = FALSE
    )
    expect_equal(mean_crossfit_unpenalized$summary$rule, "loss")

    dist_sq <- as.matrix(stats::dist(data))^2
    kernel_mat <- exp(-0.4 * dist_sq)
    wrapper_kernel <- relieverChangepoint::reliever_kde_l2(
      kernel_mat, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      echo = FALSE
    )
    direct_kernel <- relieverChangepoint::reliever_generic(
      kernel_mat, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_kde_l2,
      cpn_crit = "none", echo = FALSE
    )
    dispatched_kernel <- relieverChangepoint::reliever(
      kernel_mat, cpd_family = "kde_l2",
      cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN", echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_kernel, direct_kernel)
    .expect_same_wrapper_result(dispatched_kernel, wrapper_kernel)
    expect_equal(dispatched_kernel$settings$cpd_family, "kde_l2")
    expect_equal(wrapper_kernel$settings$cpn_crit, "none")
  })
})

test_that("KDE reliever wrappers match direct calls on two changes", {
  with_test_timeout({
    set.seed(20260712)
    x <- .wrapper_mean_data(n_seg = 10L, p = 2L)
    dist_sq <- as.matrix(stats::dist(x))^2
    bandwidth_vec <- c(0.8, 1.3)

    wrapper_path <- relieverChangepoint::reliever_kde_nll(
      dist_sq, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      bandwidth_vec = bandwidth_vec, var_dim = 2L, echo = FALSE
    )
    direct_path <- relieverChangepoint::reliever_generic(
      dist_sq, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_kde_nll_solpath,
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      cpn_crit = "none", echo = FALSE
    )
    dispatched_path <- relieverChangepoint::reliever(
      dist_sq, cpd_family = "kde_nll",
      cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      bandwidth_vec = bandwidth_vec, var_dim = 2L, echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_path, direct_path)
    .expect_same_wrapper_result(dispatched_path, wrapper_path)
    expect_equal(dispatched_path$settings$cpd_family, "kde_nll")
    expect_equal(wrapper_path$settings$cpn_crit, "none")
    expect_equal(nrow(summary(wrapper_path)), 0L)

    wrapper_crossfit <- relieverChangepoint::reliever_kde_nll_crossfit(
      dist_sq, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    direct_crossfit <- relieverChangepoint::reliever_generic(
      dist_sq, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_kde_nll_crossfit,
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      nfolds = 2L, fold_type = "blk", cpn_crit = "none", echo = FALSE
    )
    dispatched_crossfit <- relieverChangepoint::reliever(
      dist_sq, cpd_family = "kde_nll_crossfit",
      cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_crossfit, direct_crossfit)
    .expect_same_wrapper_result(dispatched_crossfit, wrapper_crossfit)
    expect_equal(dispatched_crossfit$settings$cpd_family, "kde_nll_crossfit")
  })
})

test_that("KDE NLL raw-data input matches precomputed squared distances", {
  with_test_timeout({
    set.seed(20260719)
    x <- .wrapper_mean_data(n_seg = 8L, p = 2L)
    dist_sq <- as.matrix(stats::dist(x))^2
    bandwidth_vec <- c(0.7, 1.4)

    raw_gaussian <- relieverChangepoint::reliever(
      x, cpd_family = "kde_nll",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN",
      bandwidth_vec = bandwidth_vec
    )
    distance_gaussian <- relieverChangepoint::reliever(
      dist_sq, cpd_family = "kde_nll",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN",
      bandwidth_vec = bandwidth_vec, var_dim = 2L
    )
    expect_same_cpd_path(raw_gaussian, distance_gaussian, tolerance = 1e-12)
    expect_identical(
      raw_gaussian$settings$input_spec$original_form, "raw_data"
    )
    expect_identical(
      distance_gaussian$settings$input_spec$original_form,
      "squared_distance"
    )

    raw_laplace <- relieverChangepoint::reliever(
      x, cpd_family = "kde_nll_crossfit",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN",
      kernel = "laplace", bandwidth_vec = bandwidth_vec,
      nfolds = 2L, fold_type = "blk"
    )
    distance_laplace <- relieverChangepoint::reliever_generic(
      dist_sq,
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_kde_nll_crossfit,
      kernel = "laplace", bandwidth_vec = bandwidth_vec, var_dim = 2L,
      nfolds = 2L, fold_type = "blk"
    )
    expect_same_cpd_path(raw_laplace, distance_laplace, tolerance = 1e-12)
    expect_identical(raw_laplace$settings$input_spec$distance_power, 1L)
  })
})

test_that("KDE-L2 builds exactly one fixed Gram matrix", {
  with_test_timeout({
    set.seed(20260719)
    x <- .wrapper_mean_data(n_seg = 8L, p = 2L)
    calls <- 0L
    custom_kernel <- function(x, y, scale) {
      calls <<- calls + 1L
      exp(-scale * relieverChangepoint:::.kernel_dist_sq(x, y))
    }

    from_raw <- relieverChangepoint::reliever(
      x, cpd_family = "kde_l2",
      kernel = custom_kernel, kernel_args = list(scale = 0.4),
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN"
    )
    expect_identical(calls, 1L)

    kernel_mat <- exp(
      -0.4 * relieverChangepoint:::.kernel_dist_sq(x)
    )
    from_gram <- relieverChangepoint::reliever(
      kernel_mat, cpd_family = "kde_l2",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN"
    )
    expect_same_cpd_path(from_raw, from_gram, tolerance = 1e-12)
    expect_identical(from_raw$settings$input_spec$kernel_name, "custom")

    matern <- relieverChangepoint::reliever(
      x, cpd_family = "kde_l2",
      kernel = "matern32", bandwidth = 1.2,
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN"
    )
    matern_mat <- relieverChangepoint:::.kernel_l2_matrix(
      x, kernel = "matern32", bandwidth = 1.2
    )
    matern_precomputed <- relieverChangepoint::reliever(
      matern_mat, cpd_family = "kde_l2",
      cpn_max = 1L, dm = 6L, cov_rate = 0.85, method = "SN"
    )
    expect_same_cpd_path(matern, matern_precomputed, tolerance = 1e-12)
  })
})

test_that("kernel wrappers reject ambiguous or incompatible inputs", {
  with_test_timeout({
    x <- matrix(seq_len(18L), nrow = 9L, ncol = 2L)
    dist_sq <- as.matrix(stats::dist(x))^2

    expect_error(
      relieverChangepoint::reliever_kde_nll(
        dist_sq, cpn_max = 1L, dm = 3L, bandwidth_vec = 1
      ),
      "looks like a squared-distance matrix"
    )
    expect_error(
      relieverChangepoint::reliever_kde_nll(
        dist_sq, cpn_max = 1L, dm = 3L, bandwidth_vec = 1,
        input_type = "distance", var_dim = 1.5
      ),
      "positive integer"
    )
    expect_error(
      relieverChangepoint::reliever_kde_nll(
        x, cpn_max = 1L, dm = 3L, bandwidth_vec = 1,
        kernel = "polynomial"
      ),
      "should be one of"
    )
    expect_error(
      relieverChangepoint::reliever_kde_l2(
        x, cpn_max = 1L, dm = 3L, kernel = "gaussian"
      ),
      "bandwidth must be"
    )
    expect_error(
      relieverChangepoint::reliever_kde_l2(
        x, cpn_max = 1L, dm = 3L, bandwidth = 1
      ),
      "bandwidth requires kernel"
    )
    expect_silent(
      raw_as_features <- relieverChangepoint::reliever(
        x, cpd_family = "kde_l2",
        cpn_max = 1L, dm = 3L, cov_rate = 0.85, method = "SN"
      )
    )
    expect_identical(
      raw_as_features$settings$input_spec$original_form,
      "precomputed_features"
    )
    expect_silent(
      relieverChangepoint::reliever(
        exp(-dist_sq / 2), cpd_family = "kde_l2",
        cpn_max = 1L, dm = 3L, cov_rate = 0.85, method = "SN"
      )
    )
  })
})

test_that("reg_fun_nmcd reliever wrapper matches direct call on two changes", {
  with_test_timeout({
    set.seed(20260712)
    x <- c(
      stats::rnorm(12L, mean = 0),
      stats::rnorm(12L, mean = 2),
      stats::rnorm(12L, mean = -2)
    )

    wrapper_res <- relieverChangepoint::reliever_nmcd(
      x, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN", echo = FALSE
    )
    direct_res <- relieverChangepoint::reliever_generic(
      x, cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_nmcd,
      cpn_crit = "none", echo = FALSE
    )
    dispatched_res <- relieverChangepoint::reliever(
      x, cpd_family = "nmcd",
      cpn_max = 2L, dm = 8L, cov_rate = 0.8, method = "SN", echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_res, direct_res)
    .expect_same_wrapper_result(dispatched_res, wrapper_res)
    expect_equal(dispatched_res$settings$cpd_family, "nmcd")
    expect_equal(wrapper_res$settings$cpn_crit, "none")

    eval_data <- x + stats::rnorm(length(x), sd = 0.1)
    evaluated <- relieverChangepoint::evaluate_reliever_segments(
      wrapper_res,
      data = x,
      eval_data = eval_data,
      eval_index = seq_along(x)
    )
    k0 <- evaluated$candidates$K == 0L
    stacked <- c(x, eval_data)
    expected_k0 <- sum(relieverChangepoint::reg_fun_nmcd(
      stacked,
      l = 1L, r = length(x),
      l_end = length(x) + 1L, r_end = length(stacked),
      sort_X = sort(x)
    )$loss)
    expect_equal(evaluated$candidates$eval_loss[k0], expected_k0)
  })
})

test_that("classifier ReCV reliever wrappers match direct calls on two changes", {
  skip_if_not_full_tests("optional classifier wrapper equivalence checks")
  with_test_timeout({
    set.seed(20260712)
    data <- .wrapper_mean_data(n_seg = 10L, p = 3L)

    testthat::skip_if_not_installed("ranger")
    wrapper_ranger <- suppressWarnings(relieverChangepoint::reliever_ranger_crossfit(
      data, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "op",
      hyper_set = data.frame(num.trees = 5L, min.node.size = 2L),
      ranger_args = list(seed = 20260712L),
      echo = FALSE
    ))
    direct_ranger <- suppressWarnings(relieverChangepoint::reliever_generic(
      data, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_ranger_crossfit,
      nfolds = 2L, fold_type = "op",
      hyper_set = data.frame(num.trees = 5L, min.node.size = 2L),
      ranger_args = list(seed = 20260712L),
      cpn_crit = "none",
      echo = FALSE
    ))
    dispatched_ranger <- suppressWarnings(relieverChangepoint::reliever(
      data, cpd_family = "ranger_crossfit",
      cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "op",
      hyper_set = data.frame(num.trees = 5L, min.node.size = 2L),
      ranger_args = list(seed = 20260712L), echo = FALSE
    ))
    .expect_same_wrapper_result(wrapper_ranger, direct_ranger)
    .expect_same_wrapper_result(dispatched_ranger, wrapper_ranger)
    expect_equal(dispatched_ranger$settings$cpd_family, "ranger_crossfit")
  })
})

test_that("MLP ReCV reliever wrapper matches direct call on two changes", {
  skip_if_not_full_tests("optional MLP wrapper equivalence check")
  with_test_timeout({
    testthat::skip_if_not_installed("nnet")
    set.seed(20260712)
    data <- .wrapper_mean_data(n_seg = 10L, p = 3L)

    set.seed(20260712)
    wrapper_mlp <- relieverChangepoint::reliever_mlp_crossfit(
      data, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "op",
      hyper_set = data.frame(size = 2L, maxit = 20L),
      echo = FALSE
    )
    set.seed(20260712)
    direct_mlp <- relieverChangepoint::reliever_generic(
      data, cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      reg_fun = relieverChangepoint::reg_fun_mlp_crossfit,
      nfolds = 2L, fold_type = "op",
      hyper_set = data.frame(size = 2L, maxit = 20L),
      cpn_crit = "none",
      echo = FALSE
    )
    set.seed(20260712)
    dispatched_mlp <- relieverChangepoint::reliever(
      data, cpd_family = "mlp_crossfit",
      cpn_max = 2L, dm = 7L, cov_rate = 0.8, method = "SN",
      nfolds = 2L, fold_type = "op",
      hyper_set = data.frame(size = 2L, maxit = 20L), echo = FALSE
    )
    .expect_same_wrapper_result(wrapper_mlp, direct_mlp, tolerance = 1e-7)
    .expect_same_wrapper_result(dispatched_mlp, wrapper_mlp, tolerance = 1e-7)
    expect_equal(dispatched_mlp$settings$cpd_family, "mlp_crossfit")
  })
})
