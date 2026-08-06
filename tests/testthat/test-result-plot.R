.capture_reliever_plot <- function(expr) {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  force(expr)
}

.plot_test_mean_data <- function() {
  rbind(
    matrix(stats::rnorm(20L * 3L), 20L, 3L),
    matrix(stats::rnorm(20L * 3L, mean = 3), 20L, 3L),
    matrix(stats::rnorm(20L * 3L, mean = -3), 20L, 3L)
  )
}

.plot_test_lasso_data <- function() {
  n_seg <- 30L
  p <- 5L
  x <- matrix(stats::rnorm(3L * n_seg * p), 3L * n_seg, p)
  beta <- rbind(
    matrix(c(2, -1, 0, 0, 0), n_seg, p, byrow = TRUE),
    matrix(c(-2, 1, 0, 0, 0), n_seg, p, byrow = TRUE),
    matrix(c(2, -1, 0, 0, 0), n_seg, p, byrow = TRUE)
  )
  cbind(rowSums(x * beta) + stats::rnorm(3L * n_seg, sd = 0.3), x)
}

test_that("Reliever plots reuse candidate losses and selection scores", {
  with_test_timeout({
    set.seed(2026)
    mean_fit <- relieverChangepoint::reliever_mean(
      .plot_test_mean_data(), cpn_max = 4L, dm = 8L, cov_rate = 0.6,
      method = "SN", cpn_crit = "rss_sic", echo = FALSE
    )
    mean_path <- .capture_reliever_plot(plot(mean_fit))
    expect_named(
      mean_path,
      c("run_id", "K", "cpd", "loss", "candidate_id", "score", "selected")
    )
    expect_equal(mean_path$K[mean_path$selected], summary(mean_fit)$K_hat)

    no_selection_fit <- relieverChangepoint::reliever_generic(
      .plot_test_mean_data(), reg_fun = relieverChangepoint::reg_fun_mean,
      cpn_max = 4L, dm = 8L, cov_rate = 0.6, method = "SN"
    )
    no_selection_path <- .capture_reliever_plot(plot(no_selection_fit))
    expect_false(any(no_selection_path$selected))

    set.seed(2026)
    lasso_data <- .plot_test_lasso_data()
    lam_set <- c(8, 3, 1, 0.3)
    lasso_fit <- relieverChangepoint::reliever_lasso(
      lasso_data, cpn_max = 4L, dm = 10L, cov_rate = 0.55,
      method = "SN", lam_set = lam_set, echo = FALSE
    )
    hyper_path <- .capture_reliever_plot(plot(
      lasso_fit, x_axis = "hyperparameter", K = 2L,
      cpn_crit = "loss", log = "x"
    ))
    expect_equal(hyper_path$hyper_value, sort(lam_set))
    expect_equal(hyper_path$score, hyper_path$loss)
    expect_false(any(hyper_path$selected))

    boundary_fit <- lasso_fit
    boundary_rows <- boundary_fit$cpd_path$candidates$K == 2L
    boundary_meta <- match(
      boundary_fit$cpd_path$candidates$run_id[boundary_rows],
      boundary_fit$run_meta$run_id
    )
    boundary_fit$cpd_path$candidates$loss[boundary_rows] <- as.numeric(unlist(
      boundary_fit$run_meta$hyper_value[boundary_meta], use.names = FALSE
    ))
    expect_silent(.capture_reliever_plot(plot(
      boundary_fit, x_axis = "hyperparameter", K = 2L,
      cpn_crit = "loss", log = "x"
    )))

    no_selection_lasso <- lasso_fit
    no_selection_lasso$settings$cpn_crit <- "none"
    no_selection_hyper <- .capture_reliever_plot(plot(
      no_selection_lasso, x_axis = "hyperparameter", K = 2L,
      cpn_crit = "none", log = "x"
    ))
    expect_false(any(no_selection_hyper$selected))
    expect_error(
      plot(no_selection_lasso, x_axis = "hyperparameter"),
      "K must be supplied"
    )

    expected <- merge(
      lasso_fit$cpd_path$candidates[
        lasso_fit$cpd_path$candidates$K == 2L,
        c("run_id", "loss")
      ],
      lasso_fit$run_meta[c("run_id", "hyper_value")],
      by = "run_id"
    )
    expected$hyper_value <- unlist(expected$hyper_value, use.names = FALSE)
    expected <- expected[order(expected$hyper_value), ]
    expect_equal(hyper_path$loss, expected$loss)

    expect_error(
      .capture_reliever_plot(plot(mean_fit, x_axis = "hyperparameter", K = 2L)),
      "hyper_value"
    )
    expect_error(
      plot(mean_fit, x_axis = "search_value"),
      "indexed by K"
    )
    expect_error(plot(mean_fit, K = 2L), "only")
  }, seconds = 12)
})

test_that("fit and outer-CV plots support categorical hyperparameters", {
  with_test_timeout({
    categorical_loss <- function(data, l, r, l_end = l, r_end = r,
                                  save_model = FALSE,
                                  is_virtual_run = FALSE) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 2L,
          loss_output_meta = data.frame(
            row_type = "custom",
            hyper_id = 1:2,
            hyper_value = c("small", "large"),
            loss_kind = "squared"
          )
        ))
      }
      center <- mean(data[l:r])
      loss <- (data[l_end:r_end] - center)^2
      list(loss = cbind(loss, loss + 1), model = NULL)
    }

    data <- seq_len(24L)
    fit <- relieverChangepoint::reliever_generic(
      data, reg_fun = categorical_loss,
      cpn_max = 1L, dm = 4L, cov_rate = 1,
      cpn_crit = "none", cache_backend = "by_cost_mat"
    )
    fit_path <- .capture_reliever_plot(plot(
      fit, x_axis = "hyperparameter", K = 0L,
      run_ids = 1:2, cpn_crit = "loss"
    ))
    expect_identical(fit_path$hyper_value, c("small", "large"))

    cv_fit <- relieverChangepoint::cv.reliever_generic(
      data, reg_fun = categorical_loss,
      cpn_max = 0L, dm = 4L, cov_rate = 1,
      nfolds = 3L, cache_backend = "by_cost_mat"
    )
    cv_path <- .capture_reliever_plot(plot(
      cv_fit, x_axis = "hyperparameter", K = 0L,
      run_ids = 1:2, show_se = FALSE
    ))
    expect_identical(cv_path$hyper_value, c("small", "large"))
  }, seconds = 8)
})

test_that("plot_reliever_data draws an explicit selected segmentation", {
  with_test_timeout({
    set.seed(2026)
    data <- .plot_test_mean_data()
    colnames(data) <- c("first", "second", "third")
    fit <- relieverChangepoint::reliever_mean(
      data, cpn_max = 4L, dm = 8L, cov_rate = 0.6,
      method = "SN", cpn_crit = "rss_sic", echo = FALSE
    )

    plotted <- .capture_reliever_plot(
      relieverChangepoint::plot_reliever_data(
        fit, data, columns = c("first", "second")
      )
    )
    expect_equal(names(plotted), c("index", "first", "second"))
    expect_equal(plotted$index, seq_len(nrow(data)))
    expect_equal(
      attr(plotted, "cpd_hat"),
      as.integer(summary(fit)$cpd_hat[[1L]])
    )

    cv_fit <- relieverChangepoint::cv.reliever(
      data, cpn_max = 4L, dm = 8L, cov_rate = 0.6,
      method = "SN", nfolds = 2L, echo = FALSE
    )
    cv_plot <- .capture_reliever_plot(
      relieverChangepoint::plot_reliever_data(cv_fit, data[, 1L])
    )
    expect_equal(
      attr(cv_plot, "cpd_hat"),
      as.integer(summary(cv_fit)$cpd_hat[[1L]])
    )

    conservative <- relieverChangepoint::select_by_run(
      fit, cpn_crit = "rss_sic"
    )
    selected_plot <- .capture_reliever_plot(
      relieverChangepoint::plot_reliever_data(
        fit, data[, 1L], selection = conservative
      )
    )
    expect_equal(
      attr(selected_plot, "cpd_hat"),
      as.integer(conservative$cpd_hat[[1L]])
    )

    multi_summary <- fit
    multi_summary$summary <- rbind(fit$summary, fit$summary)
    expect_error(
      relieverChangepoint::plot_reliever_data(multi_summary, data),
      "2 summary rows"
    )
    empty_summary <- fit
    empty_summary$summary <- fit$summary[0L, , drop = FALSE]
    expect_error(
      relieverChangepoint::plot_reliever_data(empty_summary, data),
      "select_by_run.*select_across_runs.*select_holdout"
    )
    expect_silent(.capture_reliever_plot(
      relieverChangepoint::plot_reliever_data(
        multi_summary, data, selection = 2L
      )
    ))
    expect_error(
      relieverChangepoint::plot_reliever_data(fit, data[-1L, ]),
      "must match"
    )
    expect_error(
      relieverChangepoint::plot_reliever_data(fit, data, columns = "missing"),
      "not present"
    )
  }, seconds = 12)
})

test_that("outer-CV plots identify the same held-out-loss minimum", {
  with_test_timeout({
    set.seed(2026)
    mean_fit <- relieverChangepoint::cv.reliever(
      .plot_test_mean_data(), cpn_max = 4L, dm = 8L, cov_rate = 0.6,
      method = "SN", nfolds = 2L, echo = FALSE
    )
    mean_path <- .capture_reliever_plot(plot(mean_fit))
    expect_true(all(c("K", "cv_mean", "cv_se", "selected") %in%
                      names(mean_path)))
    expect_equal(mean_path$K[mean_path$selected], summary(mean_fit)$K_hat)
    expect_match(
      capture.output(print(mean_fit))[2L],
      "Method: SN.*Folds: 2"
    )

    set.seed(2026)
    lasso_data <- .plot_test_lasso_data()
    lam_set <- c(8, 3, 1, 0.3)
    lasso_fit <- relieverChangepoint::cv.reliever(
      X = lasso_data[, -1], y = lasso_data[, 1],
      cpd_family = "lasso", cpn_max = 4L, dm = 10L, cov_rate = 0.55,
      method = "SN", nfolds = 2L, lam_set = lam_set, echo = FALSE
    )
    hyper_path <- .capture_reliever_plot(plot(
      lasso_fit, x_axis = "hyperparameter", K = 2L,
      show_se = FALSE, log = "x"
    ))
    expect_equal(hyper_path$hyper_value, sort(lam_set))
    expected_best <- which.min(hyper_path$cv_mean)
    expect_identical(which(hyper_path$selected), expected_best)

    boundary_fit <- lasso_fit
    k2 <- boundary_fit$cv_loss$K == 2L
    k2_hyper <- as.numeric(unlist(
      boundary_fit$cv_loss$hyper_value[k2], use.names = FALSE
    ))
    boundary_fit$cv_loss$cv_mean[k2] <- k2_hyper
    expect_warning(
      .capture_reliever_plot(plot(
        boundary_fit, x_axis = "hyperparameter", K = 2L,
        show_se = FALSE, log = "x"
      )),
      "lower endpoint"
    )
    expect_error(plot(lasso_fit, K = 2L), "only")
    expect_error(
      plot(lasso_fit, x_axis = "hyperparameter", K = -1L),
      "non-negative"
    )
  }, seconds = 15)
})

test_that("plots support ReCV defaults and selector-indexed search paths", {
  with_test_timeout({
    set.seed(2026)
    lasso_data <- .plot_test_lasso_data()
    crossfit_fit <- relieverChangepoint::reliever_lasso_crossfit(
      lasso_data, cpn_max = 3L, dm = 10L, cov_rate = 0.55,
      method = "SN", lam_set = c(3, 1), nfolds = 2L,
      fold_type = "blk",
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      echo = FALSE
    )
    crossfit_path <- .capture_reliever_plot(plot(crossfit_fit))
    explicit_recv_path <- .capture_reliever_plot(plot(
      crossfit_fit, run_type = "recv"
    ))
    default_runs <- crossfit_fit$run_meta$run_id[
      crossfit_fit$run_meta$default_selection
    ]
    expect_identical(unique(crossfit_path$run_id), default_runs)
    expect_identical(crossfit_path, explicit_recv_path)
    expect_equal(
      crossfit_fit$run_meta$row_type[
        match(default_runs, crossfit_fit$run_meta$run_id)
      ],
      "recv"
    )

    unmarked_fit <- crossfit_fit
    unmarked_fit$run_meta$default_selection <- FALSE
    unmarked_path <- .capture_reliever_plot(plot(unmarked_fit))
    expect_setequal(
      unique(unmarked_path$run_id), unmarked_fit$run_meta$run_id
    )
    incv_runs <- crossfit_fit$run_meta$run_id[
      crossfit_fit$run_meta$row_type == "incv"
    ]
    incv_path <- .capture_reliever_plot(plot(
      crossfit_fit, run_type = "incv", cpn_crit = "loss"
    ))
    expect_identical(unique(incv_path$run_id), incv_runs)
    expect_error(
      .capture_reliever_plot(plot(
        crossfit_fit,
        run_ids = incv_runs,
        run_type = "incv",
        cpn_crit = "loss"
      )),
      "only one of run_ids and run_type"
    )

    fixed_cf_runs <- crossfit_fit$run_meta$run_id[
      crossfit_fit$run_meta$row_type == "crossfit_homo_hyper"
    ]
    expect_error(
      .capture_reliever_plot(plot(
        crossfit_fit, x_axis = "hyperparameter", K = 2L,
        cpn_crit = "loss", log = "x"
      )),
      "non-missing hyper_value"
    )
    crossfit_hyper_path <- suppressWarnings(.capture_reliever_plot(plot(
      crossfit_fit, x_axis = "hyperparameter", K = 2L,
      run_type = "crossfit_homo_hyper", cpn_crit = "loss", log = "x"
    )))
    expect_setequal(unique(crossfit_hyper_path$run_id), fixed_cf_runs)
    expect_equal(crossfit_hyper_path$hyper_value, c(1, 3))
    expect_equal(sum(crossfit_hyper_path$selected), 1L)

    boundary_fit <- crossfit_fit
    fixed_meta <- boundary_fit$run_meta[
      boundary_fit$run_meta$row_type == "crossfit_homo_hyper", , drop = FALSE
    ]
    fixed_k2 <- boundary_fit$cpd_path$candidates$run_id %in%
      fixed_meta$run_id & boundary_fit$cpd_path$candidates$K == 2L
    meta_id <- match(
      boundary_fit$cpd_path$candidates$run_id[fixed_k2], fixed_meta$run_id
    )
    boundary_fit$cpd_path$candidates$loss[fixed_k2] <- as.numeric(unlist(
      fixed_meta$hyper_value[meta_id], use.names = FALSE
    ))
    expect_warning(
      .capture_reliever_plot(plot(
        boundary_fit, x_axis = "hyperparameter", K = 2L,
        run_type = "crossfit_homo_hyper", cpn_crit = "loss", log = "x"
      )),
      "lower endpoint"
    )

    set.seed(2026)
    mean_data <- .plot_test_mean_data()
    wbs_fit <- relieverChangepoint::reliever_mean(
      mean_data, cpn_max = 5L, dm = 8L, cov_rate = 0.6,
      method = "WBS", M = 30L, wbs_seed = 2026L,
      wbs_stop_crit = c(0, 1, 10, 100), cpn_crit = "none",
      echo = FALSE
    )
    wbs_path <- .capture_reliever_plot(plot(wbs_fit))
    active_ids <- unique(wbs_fit$cpd_path$selector_map$candidate_id)
    expect_setequal(wbs_path$candidate_id, active_ids)
    expect_false(any(wbs_path$selected))
    wbs_search_path <- .capture_reliever_plot(plot(
      wbs_fit, x_axis = "search_value"
    ))
    expect_equal(
      sort(wbs_search_path$wbs_stop_crit),
      c(0, 1, 10, 100)
    )
    expect_false(any(wbs_search_path$selected))

    pelt_fit <- relieverChangepoint::reliever_mean(
      mean_data, cpn_max = 5L, dm = 8L, cov_rate = 0.6,
      method = "PELT", pen_val = c(1, 10, 100),
      cpn_crit = "none", echo = FALSE
    )
    pelt_path <- .capture_reliever_plot(plot(pelt_fit))
    expect_true(all(pelt_path$candidate_id %in%
                      pelt_fit$cpd_path$selector_map$candidate_id))
    expect_false(any(pelt_path$selected))
    pelt_search_path <- .capture_reliever_plot(plot(
      pelt_fit, x_axis = "search_value"
    ))
    expect_equal(sort(pelt_search_path$pen_val), c(1, 10, 100))
    expect_false(any(pelt_search_path$selected))

    cv_wbs <- relieverChangepoint::cv.reliever(
      mean_data, cpn_max = 5L, dm = 8L, cov_rate = 0.6,
      method = "WBS", M = 30L, wbs_seed = 2026L,
      wbs_stop_crit = c(0, 1, 10, 100), nfolds = 2L,
      echo = FALSE
    )
    cv_wbs_path <- .capture_reliever_plot(plot(cv_wbs))
    expect_equal(cv_wbs_path$K[cv_wbs_path$selected],
                 summary(cv_wbs)$K_hat)
    expect_equal(nrow(cv_wbs_path), nrow(cv_wbs$cv_loss))
    cv_wbs_search <- .capture_reliever_plot(plot(
      cv_wbs, x_axis = "search_value"
    ))
    expect_equal(sort(cv_wbs_search$wbs_stop_crit), c(0, 1, 10, 100))
    expect_equal(
      cv_wbs_search$K[cv_wbs_search$selected],
      summary(cv_wbs)$K_hat
    )

    cv_pelt <- relieverChangepoint::cv.reliever(
      mean_data, cpn_max = 5L, dm = 8L, cov_rate = 0.6,
      method = "PELT", pen_val = c(1, 10, 100), nfolds = 2L,
      echo = FALSE
    )
    cv_pelt_search <- .capture_reliever_plot(plot(
      cv_pelt, x_axis = "search_value"
    ))
    expect_equal(sort(cv_pelt_search$pen_val), c(1, 10, 100))
    expect_equal(
      cv_pelt_search$K[cv_pelt_search$selected],
      summary(cv_pelt)$K_hat
    )
  }, seconds = 15)
})

test_that("outer-CV hyperparameter plots default to comparable crossfit rows", {
  with_test_timeout({
    run_meta <- data.frame(
      run_id = 1:6,
      row_type = c(
        "recv", "incv", "crossfit_homo_hyper", "crossfit_homo_hyper",
        "lasso", "lasso"
      ),
      hyper_value = c(NA, NA, 0.1, 0.2, 0.1, 0.2)
    )
    fit <- list(
      cv_loss = data.frame(
        row_type = run_meta$row_type,
        hyper_value = run_meta$hyper_value,
        K = 1L,
        cv_mean = c(9, 8, 3, 2, 4, 5),
        cv_se = 0.1,
        run_id = run_meta$run_id
      ),
      full_data_fit = list(
        run_meta = run_meta,
        cpd_path = list(
          candidates = data.frame(run_id = run_meta$run_id)
        )
      )
    )
    class(fit) <- c("cv_reliever_result", "list")

    path <- suppressWarnings(.capture_reliever_plot(plot(
      fit, x_axis = "hyperparameter", K = 1L, show_se = FALSE
    )))
    expect_identical(unique(path$run_id), 3:4)
    expect_identical(which(path$selected), 2L)
    in_path <- suppressWarnings(.capture_reliever_plot(plot(
      fit,
      x_axis = "hyperparameter",
      K = 1L,
      run_type = "lasso",
      show_se = FALSE
    )))
    expect_identical(unique(in_path$run_id), 5:6)
    expect_error(
      plot(
        fit,
        run_ids = 5:6,
        run_type = "lasso",
        show_se = FALSE
      ),
      "only one of run_ids and run_type"
    )
  }, seconds = 4)
})

test_that("two-step results reject a misleading loss-path plot", {
  with_test_timeout({
    set.seed(2026)
    fit <- relieverChangepoint::twostep(
      .plot_test_mean_data(), reg_fun = relieverChangepoint::reg_fun_kde_l2,
      cpn_max = 3L, dm = 8L, method = "BS", echo = FALSE
    )
    expect_error(plot(fit), "whole-segmentation losses")
  }, seconds = 4)
})

test_that("result printing stays compact for long model paths", {
  with_test_timeout({
    set.seed(2026)
    fit <- relieverChangepoint::reliever_mean(
      .plot_test_mean_data(), cpn_max = 3L, dm = 8L, cov_rate = 0.6,
      method = "SN", echo = FALSE
    )
    fit$summary <- fit$summary[rep(1L, 12L), , drop = FALSE]
    rownames(fit$summary) <- NULL

    compact <- capture.output(print(fit))
    expect_true(any(grepl("Showing first 10 of 12 rows", compact)))
    full <- capture.output(print(fit, max_rows = Inf))
    expect_false(any(grepl("Showing first", full)))
    expect_error(print(fit, max_rows = 0L), "positive integer or Inf")
  }, seconds = 4)
})
