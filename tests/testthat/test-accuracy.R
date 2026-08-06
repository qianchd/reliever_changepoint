.accuracy_cp_error <- function(est, truth) {
  relieverChangepoint::cp_error(
    as.integer(est[is.finite(est) & est > 0]),
    as.integer(truth)
  )
}

.accuracy_best_k_error <- function(result, truth) {
  tab <- cpd_path_candidates(result)
  tab <- tab[tab$K == length(truth), , drop = FALSE]
  if (nrow(tab) == 0L) {
    return(Inf)
  }
  min(vapply(tab$cpd, .accuracy_cp_error, numeric(1L), truth = truth))
}

.expect_selected_cpd_close <- function(result, truth, tolerance = 10L) {
  errs <- vapply(
    seq_along(result$summary$K_hat),
    function(run_id) .accuracy_cp_error(result$summary$cpd_hat[[run_id]], truth),
    numeric(1L)
  )
  expect_true(any(result$summary$K_hat == length(truth) & errs <= tolerance))
}

.expect_path_has_cpd_close <- function(result, truth, tolerance = 10L) {
  expect_lte(.accuracy_best_k_error(result, truth), tolerance)
}

.accuracy_mean_data <- function(n_seg = 40L, p = 3L, sd = 0.25) {
  rbind(
    matrix(stats::rnorm(n_seg * p, mean = 0, sd = sd), n_seg, p),
    matrix(stats::rnorm(n_seg * p, mean = 4, sd = sd), n_seg, p),
    matrix(stats::rnorm(n_seg * p, mean = -4, sd = sd), n_seg, p)
  )
}

.accuracy_lasso_data <- function(n_seg = 40L, p = 20L, sd = 0.5) {
  n <- 3L * n_seg
  x <- matrix(stats::rnorm(n * p), n)
  beta <- rbind(
    matrix(c(2.5, -2, 1.5, -1, 1, rep(0, p - 5L)), n_seg, p,
           byrow = TRUE),
    matrix(c(-2.5, 2, -1.5, 1, -1, rep(0, p - 5L)), n_seg, p,
           byrow = TRUE),
    matrix(c(2, -2.2, 1.8, -1.2, 1.2, rep(0, p - 5L)), n_seg, p,
           byrow = TRUE)
  )
  y <- rowSums(x * beta) + stats::rnorm(n, sd = sd)
  cbind(y, x)
}

.accuracy_kde_dist_sq <- function(n_seg = 300L, p = 5L) {
  x <- rbind(
    matrix(stats::rnorm(n_seg * p, mean = 0, sd = 0.5), n_seg, p),
    matrix(stats::rnorm(n_seg * p, mean = 4, sd = 0.5), n_seg, p),
    matrix(stats::rnorm(n_seg * p, mean = -4, sd = 0.5), n_seg, p)
  )
  as.matrix(stats::dist(x))^2
}

test_that("applied mean workflows select two changes at n = 900", {
  with_test_timeout({
    set.seed(2026)
    n_seg <- 300L
    truth <- c(300L, 600L)
    data <- rbind(
      matrix(stats::rnorm(n_seg * 5L, mean = 0, sd = 0.5), n_seg, 5L),
      matrix(stats::rnorm(n_seg * 5L, mean = 4, sd = 0.5), n_seg, 5L),
      matrix(stats::rnorm(n_seg * 5L, mean = -4, sd = 0.5), n_seg, 5L)
    )

    path_fit <- relieverChangepoint::reliever(data)
    expect_equal(nrow(summary(path_fit)), 0L)

    for (fit in list(
      relieverChangepoint::cv.reliever(data),
      relieverChangepoint::cv.reliever(
        data, cpn_max = 20L, dm = 30L, cov_rate = 0.8,
        method = "WBS", M = 200L,
        wbs_stop_crit = c(5, 10, 15, 20), nfolds = 3L
      )
    )) {
      selected <- summary(fit)
      expect_equal(selected$K_hat, 2L)
      expect_cpd_error_lte(selected$cpd_hat[[1L]], truth)
    }
  }, seconds = 10)
})

test_that("explicit lasso selection gives accurate applied results at n = 900", {
  skip_if_not_full_tests("n = 900 lasso selection accuracy check")
  with_test_timeout({
    set.seed(2026)
    n <- 900L
    p <- 100L
    truth <- c(300L, 600L)
    b0 <- c(3, -2.5, 2, -1.5, 1.5, rep(0, p - 5L))
    delta <- cbind(-2 * b0, 1.8 * b0)
    data <- relieverChangepoint::dgp_linear_regression(
      n, p, truth, b0, delta, sig = 1
    )$data

    crossfit <- relieverChangepoint::reliever(
      X = data[, -1L, drop = FALSE], y = data[, 1L],
      cpd_family = "lasso_crossfit", cpn_max = 7L, dm = 30L, cov_rate = 0.8,
      nfolds = 2L, method = "SN",
      loss_output_types = c("recv", "crossfit_homo_hyper")
    )
    selected_recv <- relieverChangepoint::select_by_run(
      crossfit, run_type = "recv", cpn_crit = "loss"
    )
    outer_cv <- relieverChangepoint::cv.reliever(
      X = data[, -1L, drop = FALSE], y = data[, 1L],
      cpd_family = "lasso", cpn_max = 7L, dm = 30L, cov_rate = 0.8,
      nfolds = 2L, method = "SN"
    )

    for (selected in list(selected_recv, summary(outer_cv))) {
      expect_equal(selected$K_hat, 2L)
      expect_cpd_error_lte(selected$cpd_hat[[1L]], truth)
    }

    fixed_cf_runs <- with(
      crossfit$run_meta,
      run_id[row_type == "crossfit_homo_hyper"]
    )
    homogeneous_tuning <- relieverChangepoint::select_across_runs(
      crossfit, run_type = "crossfit_homo_hyper", cpn_crit = "loss"
    )
    expect_equal(homogeneous_tuning$K_hat, 2L)
    expect_cpd_error_lte(homogeneous_tuning$cpd_hat[[1L]], truth)
    expect_true("hyper_value" %in% names(homogeneous_tuning))
  }, seconds = 120)
})

test_that("solution-path reg_fun wrappers have accurate fixed-K paths", {
  with_test_timeout({
    set.seed(20260712)
    lasso_data <- .accuracy_lasso_data()
    lasso_truth <- c(40L, 80L)

    lasso_res <- relieverChangepoint::reliever_lasso(
      lasso_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "WBS", M = 60L, wbs_seed = 20260712L,
      lam_set = c(10, 5, 2, 1, 0.5, 0.2), echo = FALSE
    )
    .expect_path_has_cpd_close(lasso_res, lasso_truth, tolerance = 3L)

    mean_data <- .accuracy_mean_data()
    mean_truth <- c(40L, 80L)
    dist_sq <- as.matrix(stats::dist(mean_data))^2

    kde_res <- relieverChangepoint::reliever_kde_nll(
      dist_sq, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "WBS", M = 60L, wbs_seed = 20260712L,
      bandwidth_vec = exp(
        seq(log(1 / sqrt(20)), log(sqrt(3 / 2)), length.out = 6L)
      ),
      var_dim = 3L, echo = FALSE
    )
    .expect_path_has_cpd_close(kde_res, mean_truth, tolerance = 1L)
  }, seconds = 25)
})

test_that("ReCV reg_fun wrappers have accurate fixed-K paths", {
  with_test_timeout({
    set.seed(20260712)
    lasso_data <- .accuracy_lasso_data()
    lasso_truth <- c(40L, 80L)

    lasso_crossfit_res <- relieverChangepoint::reliever_lasso_crossfit(
      lasso_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "BS", M = 0L,
      lam_set = c(5, 1), nfolds = 2L, fold_type = "blk",
      echo = FALSE
    )
    .expect_path_has_cpd_close(lasso_crossfit_res, lasso_truth, tolerance = 4L)

    mean_data <- .accuracy_mean_data()
    mean_truth <- c(40L, 80L)

    mean_crossfit_res <- relieverChangepoint::reliever_mean_crossfit(
      mean_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "SeedBS", M = 60L, wbs_seed = 20260712L,
      nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    .expect_path_has_cpd_close(mean_crossfit_res, mean_truth, tolerance = 1L)

    dist_sq <- as.matrix(stats::dist(mean_data))^2
    kde_cv_res <- relieverChangepoint::reliever_kde_nll_crossfit(
      dist_sq, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "SN",
      bandwidth_vec = exp(
        seq(log(1 / sqrt(20)), log(sqrt(3 / 2)), length.out = 6L)
      ),
      var_dim = 3L,
      nfolds = 2L, fold_type = "blk", echo = FALSE
    )
    .expect_path_has_cpd_close(kde_cv_res, mean_truth, tolerance = 1L)
  }, seconds = 30)
})

test_that("mean-like and nonparametric wrappers have accurate fixed-K paths", {
  with_test_timeout({
    set.seed(20260712)
    mean_data <- .accuracy_mean_data()
    mean_truth <- c(40L, 80L)
    dist_sq <- as.matrix(stats::dist(mean_data))^2
    kernel_mat <- exp(-0.05 * dist_sq)

    mean_res <- relieverChangepoint::reliever_mean(
      mean_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "SN", cache_backend = "by_loss_block", echo = FALSE
    )
    .expect_path_has_cpd_close(mean_res, mean_truth, tolerance = 1L)

    kernel_res <- relieverChangepoint::reliever_kde_l2(
      kernel_mat, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "BS", M = 0L, echo = FALSE
    )
    .expect_path_has_cpd_close(kernel_res, mean_truth, tolerance = 1L)

    nmcd_data <- c(
      stats::rnorm(40L, mean = 0, sd = 0.1),
      stats::rnorm(40L, mean = 4, sd = 0.1),
      stats::rnorm(40L, mean = -4, sd = 0.1)
    )
    nmcd_res <- relieverChangepoint::reliever_nmcd(
      nmcd_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "SN", echo = FALSE
    )
    .expect_path_has_cpd_close(nmcd_res, mean_truth, tolerance = 1L)
  }, seconds = 20)
})

test_that("PELT and OP penalty paths include accurate two-change fits", {
  with_test_timeout({
    set.seed(20260712)
    mean_data <- .accuracy_mean_data()
    mean_truth <- c(40L, 80L)
    pen_val <- c(1, 3, 10, 30, 100, 300, 1000)

    for (method in c("PELT", "OP")) {
      fit_path <- function() {
        relieverChangepoint::reliever_generic(
          mean_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
          method = method, pen_val = pen_val,
          reg_fun = relieverChangepoint::reg_fun_kde_l2,
          cache_backend = "by_loss_block", echo = FALSE
        )
      }
      res <- fit_path()
      .expect_path_has_cpd_close(res, mean_truth, tolerance = 1L)
    }
  }, seconds = 20)
})

test_that("optional classifier ReCV wrappers have accurate fixed-K paths", {
  with_test_timeout({
    skip_if_not_installed("ranger")
    set.seed(20260712)
    mean_data <- .accuracy_mean_data(n_seg = 30L)
    mean_truth <- c(30L, 60L)

    ranger_res <- suppressWarnings(relieverChangepoint::reliever_ranger_crossfit(
      mean_data, cpn_max = 7L, dm = 8L, cov_rate = 0.55,
      method = "SN", nfolds = 2L, fold_type = "op",
      hyper_set = data.frame(num.trees = 10L, min.node.size = 2L),
      ranger_args = list(seed = 20260712L), echo = FALSE
    ))
    .expect_path_has_cpd_close(ranger_res, mean_truth, tolerance = 2L)
  }, seconds = 25)
})

test_that("optional MLP ReCV wrapper has an accurate fixed-K path", {
  with_test_timeout({
    skip_if_not_installed("nnet")
    set.seed(20260712)
    mean_data <- .accuracy_mean_data(n_seg = 30L)
    mean_truth <- c(30L, 60L)

    mlp_res <- suppressWarnings(relieverChangepoint::reliever_mlp_crossfit(
      mean_data, cpn_max = 7L, dm = 8L, cov_rate = 0.55,
      method = "SN", nfolds = 2L, fold_type = "op",
      loss_output_types = c("recv", "incv"),
      hyper_set = data.frame(size = 2L, maxit = 20L), echo = FALSE
    ))
    .expect_path_has_cpd_close(mlp_res, mean_truth, tolerance = 2L)
  }, seconds = 25)
})

test_that("RSS-SIC selects reasonable K on easy squared-loss examples", {
  with_test_timeout({
    set.seed(20260712)

    mean_data <- .accuracy_mean_data()
    mean_truth <- c(40L, 80L)
    mean_res <- relieverChangepoint::reliever_mean(
      mean_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "SN", cpn_crit = "rss_sic", echo = FALSE
    )
    .expect_selected_cpd_close(mean_res, mean_truth, tolerance = 1L)

    lasso_data <- .accuracy_lasso_data()
    lasso_truth <- c(40L, 80L)
    lasso_res <- relieverChangepoint::reliever_lasso(
      lasso_data, cpn_max = 7L, dm = 10L, cov_rate = 0.55,
      method = "SN", lam_set = c(10, 5, 2, 1, 0.5, 0.2),
      cpn_crit = "rss_sic", echo = FALSE
    )
    .expect_selected_cpd_close(lasso_res, lasso_truth, tolerance = 3L)
  }, seconds = 25)
})

test_that("explicit KDE SIC paths contain two strong mean changes", {
  skip_if_not_full_tests("KDE NLL bandwidth-path accuracy check")
  with_test_timeout({
    set.seed(2026)
    p <- 5L
    dist_sq <- .accuracy_kde_dist_sq(p = p)
    path_fit <- relieverChangepoint::reliever_kde_nll(
      dist_sq, cpn_max = 7L, dm = 30L, cov_rate = 0.8,
      var_dim = p, method = "SN", echo = FALSE
    )
    path_selected <- relieverChangepoint::select_by_run(
      path_fit, cpn_crit = "sic"
    )
    path_accurate <- path_selected$K_hat == 2L & vapply(
      path_selected$cpd_hat,
      function(cpd) {
        relieverChangepoint::cp_error(cpd, c(300L, 600L)) <= 10L
      },
      logical(1L)
    )
    expect_true(any(path_accurate))
  }, seconds = 120)
})

test_that("KDE ReCV with SIC selects two strong mean changes", {
  skip_if_not_full_tests("KDE NLL ReCV accuracy check")
  with_test_timeout({
    set.seed(2026)
    p <- 5L
    dist_sq <- .accuracy_kde_dist_sq(p = p)
    fit <- relieverChangepoint::reliever_kde_nll_crossfit(
      dist_sq, cpn_max = 7L, dm = 30L, cov_rate = 0.8,
      var_dim = p, nfolds = 2L, method = "SN",
      echo = FALSE
    )
    adaptive_run <- with(
      fit$run_meta,
      run_id[row_type == "recv"]
    )
    selected_sic <- relieverChangepoint::select_by_run(
      fit, run_ids = adaptive_run, cpn_crit = "sic"
    )
    expect_equal(selected_sic$K_hat, 2L)
    expect_cpd_error_lte(selected_sic$cpd_hat[[1L]], c(300L, 600L))
  }, seconds = 120)
})

test_that("KDE NLL ReCV selects changes in a mixed distribution DGP", {
  skip_if_not_full_tests("KDE NLL mixed-distribution accuracy check")
  with_test_timeout({
    set.seed(2026)
    n <- 1000L
    p <- 10L
    signal <- 0.8
    truth <- floor(c(0.35, 0.5, 0.88) * n)
    x <- matrix(stats::rnorm(n * p), n, p)
    sigma_x <- stats::toeplitz(signal^(0:(p - 1L)))
    x[seq_len(truth[1L]), ] <-
      x[seq_len(truth[1L]), ] %*% chol(sigma_x)
    x[(truth[2L] + 1L):truth[3L], 1:5] <-
      x[(truth[2L] + 1L):truth[3L], 1:5] + signal
    fit <- relieverChangepoint::reliever_kde_nll_crossfit(
      as.matrix(stats::dist(x))^2,
      cpn_max = 10L, dm = 50L, cov_rate = 0.55,
      var_dim = p, nfolds = 5L, method = "SN",
      cpn_crit = "sic", echo = FALSE
    )
    selected <- summary(fit)
    expect_equal(selected$K_hat, 3L)
    expect_cpd_error_lte(selected$cpd_hat[[1L]], truth, tolerance = 2L)
  }, seconds = 60)
})
