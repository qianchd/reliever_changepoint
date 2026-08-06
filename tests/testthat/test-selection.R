test_that("reliever_mean defaults to no selection and supports RSS-SIC", {
  with_test_timeout({
    set.seed(1)
    for (sigma in c(0.5, 1, 2)) {
      data <- matrix(stats::rnorm(120, sd = sigma), ncol = 1L)
      res <- relieverChangepoint::reliever_mean(
        data,
        cpn_max = 3L,
        dm = 20L,
        cov_rate = 1,
        method = "SN",
        cache_backend = "by_cost_mat",
        echo = FALSE
      )
      selected <- relieverChangepoint::select_by_run(
        res, cpn_crit = "rss_sic"
      )

      expect_equal(res$settings$cpn_crit, "none")
      expect_equal(res$settings$cpn_penalty, 0)
      expect_equal(res$settings$cpn_criterion, "loss")
      expect_equal(nrow(res$summary), 0L)
      expect_equal(selected$K_hat, 0L, info = paste("sigma =", sigma))
      expect_equal(selected$rule, "rss_sic")
      expect_equal(selected$cpd_hat[[1L]], integer())
    }
    expect_error(
      relieverChangepoint::reliever_mean(
        data, cpn_max = 1L, dm = 20L, cpn_crit = "default"
      ),
      "Unknown penalty preset"
    )
    expect_silent(
      additive_fit <- relieverChangepoint::reliever_mean(
        data, cpn_max = 1L, dm = 20L, cpn_crit = "sic"
      )
    )
    expect_identical(additive_fit$settings$cpn_criterion, "loss")
    expect_silent(
      additive_selected <- relieverChangepoint::select_by_run(
        res, cpn_crit = "sic"
      )
    )
    expect_identical(additive_selected$rule, "sic")
  })
})

test_that("post-fit RSS-SIC selects strong mean changes", {
  with_test_timeout({
    set.seed(2)
    data <- matrix(
      c(
        stats::rnorm(50, mean = 0, sd = 1),
        stats::rnorm(50, mean = 3, sd = 1),
        stats::rnorm(50, mean = -3, sd = 1)
      ),
      ncol = 1L
    )
    res <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 4L,
      dm = 15L,
      cov_rate = 1,
      method = "SN",
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    selected <- relieverChangepoint::select_by_run(
      res, cpn_crit = "rss_sic"
    )

    expect_equal(selected$K_hat, 2L)
    expect_cpd_error_lte(selected$cpd_hat[[1L]], c(50L, 100L))
    expect_equal(nrow(res$summary), 0L)
  })
})

test_that("named criteria are not gated by declared loss_kind", {
  result <- list(
    cpd_path = list(
      select_by = "K",
      selector_map = NULL,
      candidates = data.frame(
        run_id = 1L,
        candidate_id = 1:2,
        K = 0:1,
        loss = c(10, 4),
        cpd = I(list(integer(), 20L))
      )
    ),
    run_meta = data.frame(
      run_id = 1L,
      loss_output_id = 1L,
      row_type = "test",
      loss_kind = "rss"
    ),
    settings = list(n = 60L, search_n = 60L)
  )

  for (loss_kind in c("rss", "negative_log_likelihood", "crossfit_loss")) {
    result$run_meta$loss_kind <- loss_kind
    expect_silent(
      relieverChangepoint::select_by_run(result, cpn_crit = "sic")
    )
    expect_silent(
      relieverChangepoint::select_by_run(result, cpn_crit = "rss_sic")
    )
  }
})

test_that("reliever returns default row-wise changepoint estimates", {
  with_test_timeout({
    data <- matrix(seq_len(40), ncol = 1L)
    res <- relieverChangepoint::reliever_generic(
      data,
      cpn_max = 2L,
      dm = 8L,
      cov_rate = 1,
      method = "SN",
      reg_fun = reg_null,
      cpn_crit = "sic",
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    selected <- relieverChangepoint::select_by_run(res, cpn_crit = "sic")

    expect_equal(res$summary$K_hat, selected$K_hat)
    expect_equal(res$summary$cpd_hat, selected$cpd_hat)
    expect_equal(length(res$summary$K_hat), nrow(res$run_meta))
    expect_equal(length(res$summary$cpd_hat), nrow(res$run_meta))
  })
})

test_that("dc_grid criteria use the original observation count", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1L)
    dc_grid <- seq(2L, 30L, by = 2L)
    res <- relieverChangepoint::reliever_generic(
      data = data,
      reg_fun = reg_null,
      cpn_max = 1L,
      dm = 2L,
      cov_rate = 0.8,
      method = "SN",
      cpn_crit = "sic",
      echo = FALSE,
      dc_grid = dc_grid
    )
    selected <- relieverChangepoint::select_by_run(res, cpn_crit = "sic")

    expect_equal(res$settings$n, nrow(data))
    expect_equal(res$settings$search_n, length(dc_grid))
    expect_equal(res$settings$cpn_penalty, log(nrow(data)))
    expect_equal(res$summary$K_hat, selected$K_hat)
    expect_equal(res$summary$cpd_hat, selected$cpd_hat)
  })
})

test_that("additive and RSS changepoint penalties use explicit scales", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1L)

    rss_sic <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 1L,
      dm = 8L,
      cov_rate = 1,
      method = "SN",
      cpn_crit = "rss_sic",
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    numeric_penalty <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 1L,
      dm = 8L,
      cov_rate = 1,
      method = "SN",
      cpn_crit = 4,
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    expect_identical(numeric_penalty$settings$cpn_crit, 4)
    loss_penalty <- relieverChangepoint::reliever_mean(
      data,
      cpn_max = 1L,
      dm = 8L,
      cov_rate = 1,
      method = "SN",
      cpn_crit = 1,
      cache_backend = "by_cost_mat",
      echo = FALSE
    )

    n <- rss_sic$settings$n
    expect_equal(rss_sic$settings$cpn_penalty, log(n))
    expect_equal(rss_sic$settings$cpn_crit, "rss_sic")
    expect_equal(rss_sic$settings$cpn_criterion, "rss_log_loss")
    expect_equal(numeric_penalty$settings$cpn_penalty, 4)
    expect_equal(numeric_penalty$settings$cpn_crit, 4)
    expect_equal(numeric_penalty$settings$cpn_criterion, "loss")
    expect_false("ic" %in% names(numeric_penalty$cpd_path$candidates))
    expect_equal(
      loss_penalty$settings$cpn_penalty,
      1
    )
    expect_equal(loss_penalty$settings$cpn_crit, 1)
    expect_equal(loss_penalty$settings$cpn_criterion, "loss")
    numeric_selected <- relieverChangepoint::select_by_run(
      numeric_penalty, cpn_crit = 4
    )
    numeric_candidate <- subset(
      numeric_penalty$cpd_path$candidates,
      candidate_id == numeric_selected$candidate_id
    )
    expect_equal(
      numeric_selected$score,
      numeric_candidate$loss + numeric_penalty$settings$cpn_penalty *
        numeric_candidate$K
    )
    expect_error(
      relieverChangepoint::select_by_run(rss_sic, cpn_crit = "bad"),
      "Unknown penalty preset"
    )
    expect_error(
      relieverChangepoint::select_by_run(rss_sic, cpn_crit = -1),
      "non-negative"
    )
    expect_error(
      relieverChangepoint::select_by_run(rss_sic, cpn_crit = "none"),
      "does not select"
    )

    additive <- lapply(
      c("aic", "hqc", "sic"),
      relieverChangepoint:::.cpn_penalty,
      n = n
    )
    rss <- lapply(
      c("rss_aic", "rss_hqc", "rss_sic"),
      relieverChangepoint:::.cpn_penalty,
      n = n
    )
    expected_penalty <- c(2, 2 * log(log(n)), log(n))
    expect_equal(vapply(additive, `[[`, numeric(1), "value"), expected_penalty)
    expect_equal(vapply(rss, `[[`, numeric(1), "value"), expected_penalty)
    expect_equal(vapply(additive, `[[`, character(1), "criterion"),
                 rep("loss", 3L))
    expect_equal(vapply(rss, `[[`, character(1), "criterion"),
                 rep("rss_log_loss", 3L))
    expect_error(relieverChangepoint:::.cpn_penalty("4", n), "Unknown")
    expect_error(relieverChangepoint:::.cpn_penalty("hq", n), "Unknown")
    expect_error(relieverChangepoint:::.cpn_penalty("rss_hq", n), "Unknown")
    expect_equal(
      relieverChangepoint:::.selection_criterion(
        loss = c(100, 80), K = 0:1, n = n,
        penalty = log(n), criterion = "loss"
      ),
      c(100, 80 + log(n))
    )
    expect_equal(
      relieverChangepoint:::.selection_criterion(
        loss = c(100, 80), K = 0:1, n = n,
        penalty = log(n), criterion = "rss_log_loss"
      ),
      n / 2 * log(c(100, 80) / n) + log(n) * c(0, 1)
    )
    expect_error(
      relieverChangepoint:::.selection_criterion(
        loss = -1, K = 0, n = n,
        penalty = log(n), criterion = "rss_log_loss"
      ),
      "non-negative loss"
    )
  })
})

test_that("changepoint selection is done separately for each loss output", {
  with_test_timeout({
    result <- list(
      cpd_path = list(
        select_by = "K",
        selector_map = NULL,
        candidates = data.frame(
          run_id = rep(1:2, each = 3L),
          candidate_id = rep(1:3, times = 2L),
          K = rep(0:2, times = 2L),
          loss = c(100, 10, 9.9, 100, 80, 60),
          cpd = I(list(
            integer(), 11L, c(11L, 22L),
            integer(), 33L, c(33L, 44L)
          ))
        )
      ),
      run_meta = data.frame(run_id = 1:2, loss_output_id = 1:2),
      settings = list(n = 100L, search_n = 100L, cpn_crit = "sic")
    )

    selected <- relieverChangepoint::select_by_run(
      result, cpn_crit = "sic"
    )

    expect_s3_class(selected, "reliever_model_selection")
    expect_equal(length(selected$K_hat), 2L)
    expect_equal(selected$run_id, 1:2)
    expect_equal(selected$K_hat, c(1L, 2L))
    expect_equal(selected$cpd_hat[[1L]], 11L)
    expect_equal(selected$cpd_hat[[2L]], c(33L, 44L))
    expect_equal(selected$rule, rep("sic", 2L))
    compact_print <- capture.output(print(selected))
    detailed_print <- capture.output(print(selected, details = TRUE))
    expect_false(any(grepl("run_id|candidate_id", compact_print)))
    expect_true(any(grepl("run_id", detailed_print)))
    expect_true(any(grepl("candidate_id", detailed_print)))

    selected_second <- relieverChangepoint::select_by_run(
      result, run_ids = 2L, cpn_crit = "sic"
    )
    expect_equal(selected_second$run_id, 2L)
    expect_equal(selected_second$K_hat, 2L)
    expect_equal(selected_second$cpd_hat[[1L]], c(33L, 44L))

    partial_nonfinite <- result
    partial_nonfinite$cpd_path$candidates$loss[
      partial_nonfinite$cpd_path$candidates$run_id == 2L
    ] <- Inf
    selected_partial <- relieverChangepoint::select_by_run(
      result = partial_nonfinite, cpn_crit = "sic"
    )
    expect_equal(nrow(selected_partial), 2L)
    expect_equal(selected_partial$run_id, 1:2)
    expect_equal(selected_partial$K_hat[[1L]], 1L)
    expect_true(is.na(selected_partial$K_hat[[2L]]))
    expect_identical(selected_partial$cpd_hat[[2L]], NA_integer_)
    expect_true(is.na(selected_partial$score[[2L]]))
    expect_true(is.na(selected_partial$candidate_id[[2L]]))

    all_nonfinite <- result
    all_nonfinite$cpd_path$candidates$loss[] <- Inf
    selected_all <- relieverChangepoint::select_by_run(
      result = all_nonfinite, cpn_crit = "sic"
    )
    expect_equal(nrow(selected_all), 2L)
    expect_equal(selected_all$run_id, 1:2)
    expect_true(all(is.na(selected_all$K_hat)))
    expect_true(all(is.na(selected_all$score)))
    expect_true(all(is.na(selected_all$candidate_id)))

    expect_error(
      relieverChangepoint::select_by_run(result),
      "cpn_crit must be supplied"
    )
    expect_error(
      relieverChangepoint::select_by_run(result, cpn_crit = NULL),
      "must be supplied"
    )
    expect_error(
      relieverChangepoint::select_across_runs(
        result, run_ids = 1:2
      ),
      "cpn_crit must be supplied"
    )
  })
})

test_that("reliever results expose compact candidates and estimates", {
  with_test_timeout({
    data <- matrix(c(rep(0, 30L), rep(3, 30L), rep(-3, 30L)), ncol = 1L)
    result <- relieverChangepoint::reliever_mean(
      data, cpn_max = 3L, dm = 10L, cov_rate = 1, method = "SN",
      cache_backend = "by_cost_mat", echo = FALSE
    )

    expect_s3_class(result, "reliever_result")
    expect_named(
      result,
      c("summary", "cpd_path", "run_meta", "settings", "timing")
    )
    expect_named(
      result$cpd_path,
      c("select_by", "candidates", "selector_map")
    )
    expect_named(
      result$cpd_path$candidates,
      c("run_id", "K", "cpd", "loss", "candidate_id")
    )
    expect_named(
      result$summary,
      c("rule", "K_hat", "cpd_hat")
    )
    expect_s3_class(result$summary, "reliever_summary")
    expect_match(capture.output(print(result))[1L], "Reliever changepoint result")
    expect_named(summary(result), c("rule", "K_hat", "cpd_hat"))
    expect_identical(summary(result), result$summary)
    expect_false(any(grepl(
      "run_id|candidate_id", capture.output(print(result))
    )))
    expect_equal(result$settings$cpn_max, 3L)
    expect_equal(result$settings$ratio, 0.9)
    expect_equal(result$settings$cache_backend, "by_cost_mat")
    expect_false("owner_key" %in% names(result$settings))
  })
})

test_that("selector maps restrict model selection without duplicating candidates", {
  result <- list(
    cpd_path = list(
      select_by = "wbs_stop_crit",
      candidates = data.frame(
        run_id = 1L,
        candidate_id = 1:3,
        K = 0:2,
        loss = c(10, 0, 5),
        cpd = I(list(integer(), 10L, c(10L, 20L)))
      ),
      selector_map = data.frame(
        run_id = 1L,
        select_value = c(Inf, 3),
        candidate_id = c(1L, 3L)
      )
    ),
    run_meta = data.frame(run_id = 1L, loss_output_id = 1L),
    settings = list(n = 30L, search_n = 30L)
  )

  selected <- relieverChangepoint::select_by_run(
    result, cpn_crit = "loss"
  )
  expect_equal(selected$candidate_id, 3L)
  expect_equal(selected$K_hat, 2L)
  expect_error(
    relieverChangepoint::select_by_run(result, cpn_crit = "recv"),
    "Unknown penalty preset"
  )
  expect_equal(nrow(result$cpd_path$candidates), 3L)
})

test_that("default run metadata identifies preferred selection paths", {
  run_meta <- data.frame(
    run_id = 1:2,
    loss_output_id = 1:2,
    row_type = c("recv", "incv"),
    default_selection = c(TRUE, FALSE)
  )
  expect_identical(
    relieverChangepoint:::.reliever_default_run_ids(run_meta),
    1L
  )

  unmarked <- run_meta
  unmarked$default_selection <- FALSE
  expect_null(
    relieverChangepoint:::.reliever_default_run_ids(unmarked)
  )
  expect_identical(
    relieverChangepoint:::.reliever_default_run_ids(
      unmarked, unmarked = "none"
    ),
    integer()
  )
})

test_that("search-value summaries identify statistically distinct run types", {
  result <- list(
    cpd_path = list(
      select_by = "wbs_stop_crit",
      candidates = data.frame(
        run_id = 1:2,
        candidate_id = 1L,
        K = c(1L, 2L),
        loss = c(3, 2),
        cpd = I(list(10L, c(10L, 20L)))
      ),
      selector_map = data.frame(
        run_id = 1:2,
        select_value = c(4, 5),
        candidate_id = 1L
      )
    ),
    run_meta = data.frame(
      run_id = 1:2,
      loss_output_id = 1:2,
      row_type = c("first", "second"),
      default_selection = TRUE
    ),
    settings = list(cpn_crit = "none")
  )

  selected <- relieverChangepoint:::.reliever_build_summary(result)
  expect_named(
    selected,
    c("row_type", "wbs_stop_crit", "K_hat", "cpd_hat")
  )
  expect_identical(selected$row_type, c("first", "second"))
})

test_that("select_across_runs supports homogeneous crossfit tuning", {
  with_test_timeout({
    cpd <- function(a, b) list(integer(), a, c(a, b))
    cpd_rows <- c(
      cpd(10L, 11L), cpd(20L, 21L), cpd(30L, 31L),
      cpd(40L, 41L), cpd(50L, 51L), cpd(60L, 61L)
    )
    losses <- c(
      100, 8, 7,
      100, 9, 8,
      100, 5, 4,
      100, 3, 10,
      100, 50, 40,
      100, 30, 100
    )
    result <- list(
      cpd_path = list(
        select_by = "K",
        selector_map = NULL,
        candidates = data.frame(
          run_id = rep(1:6, each = 3L),
          candidate_id = rep(1:3, times = 6L),
          K = rep(0:2, times = 6L),
          loss = losses,
          cpd = I(cpd_rows)
        )
      ),
      run_meta = data.frame(
        run_id = 1:6,
        loss_output_id = 1:6,
        row_type = c(
          "recv", "incv", "crossfit_homo_hyper", "crossfit_homo_hyper",
          "lasso", "lasso"
        ),
        default_selection = c(TRUE, rep(FALSE, 5L)),
        hyper_id = c(NA, NA, 1L, 2L, 1L, 2L),
        hyper_value = c(NA, NA, 0.1, 0.2, 0.1, 0.2)
      )
    )
    result$settings <- list(n = 100L, search_n = 100L)

    by_run <- relieverChangepoint::select_by_run(
      result, run_ids = c(1L, 2L), cpn_crit = "loss"
    )
    expect_equal(by_run$run_id, 1:2)
    expect_equal(by_run$K_hat, c(2L, 2L))
    expect_identical(by_run$row_type, c("recv", "incv"))

    by_type <- relieverChangepoint::select_by_run(
      result, run_type = c("recv", "incv"), cpn_crit = "loss"
    )
    expect_equal(by_type, by_run)

    across_cf <- relieverChangepoint::select_across_runs(
      result, run_type = "crossfit_homo_hyper", cpn_crit = "loss"
    )
    expect_equal(across_cf$run_id, 4L)
    expect_equal(across_cf$K_hat, 1L)
    expect_equal(across_cf$cpd_hat[[1L]], 40L)
    expect_equal(across_cf$hyper_value, 0.2)
    expect_identical(across_cf$row_type, "crossfit_homo_hyper")
    expect_match(
      paste(capture.output(print(across_cf)), collapse = "\n"),
      "crossfit_homo_hyper"
    )

    lasso_result <- result
    lasso_result$run_meta$loss_kind <- c(
      NA, NA, NA, NA, "rss", "rss"
    )
    expect_warning(
      across_lasso <- relieverChangepoint::select_across_runs(
        lasso_result, run_type = "lasso", cpn_crit = "sic"
      ),
      "run_type = \"lasso\" may select incomparable losses.*cv.reliever",
      class = "reliever_across_run_comparability_warning"
    )
    expect_equal(across_lasso$run_id, 6L)
    expect_identical(across_lasso$rule, "sic")

    expect_warning(
      across_recv <- relieverChangepoint::select_across_runs(
        result, run_type = "recv", cpn_crit = "loss"
      ),
      "selects one run.*Use select_by_run",
      class = "reliever_single_run_selection_warning"
    )
    expect_equal(across_recv$run_id, 1L)

    expect_silent(
      relieverChangepoint::select_across_runs(
        result, run_ids = 5:6, cpn_crit = "loss"
      )
    )

    expect_error(
      relieverChangepoint::select_by_run(
        result, run_ids = 1L, run_type = "recv", cpn_crit = "loss"
      ),
      "only one of run_ids and run_type"
    )
    expect_error(
      relieverChangepoint::select_by_run(
        result, run_type = c("recv", "missing"), cpn_crit = "loss"
      ),
      "does not match.*missing.*Available values.*crossfit_homo_hyper"
    )
    expect_error(
      relieverChangepoint::select_by_run(
        result, run_type = "", cpn_crit = "loss"
      ),
      "non-empty character"
    )
    no_type <- result
    no_type$run_meta$row_type <- NULL
    expect_error(
      relieverChangepoint::select_by_run(
        no_type, run_type = "recv", cpn_crit = "loss"
      ),
      "must contain row_type"
    )
    expect_error(
      relieverChangepoint::select_across_runs(
        result, cpn_crit = "loss"
      ),
      "run_ids or run_type must be supplied"
    )
    mixed_scale <- result
    mixed_scale$run_meta$loss_kind <- c(
      "crossfit_loss", "rss",
      "crossfit_loss", "crossfit_loss", "rss", "rss"
    )
    expect_warning(
      mixed_selection <- relieverChangepoint::select_across_runs(
        mixed_scale, run_ids = 1:2, cpn_crit = "loss"
      ),
      "different loss_kind values.*losses may not be comparable.*continue",
      class = "reliever_loss_kind_comparability_warning"
    )
    expect_equal(mixed_selection$run_id, 1L)

    crossfit_result <- result
    crossfit_result$settings$cpn_crit <- "loss"
    estimated <- relieverChangepoint:::.reliever_finalize_result(crossfit_result)
    expect_false("sic" %in% estimated$summary$rule)
    expect_equal(estimated$summary$rule, "loss")
    expect_identical(estimated$summary$row_type, "recv")

    result_sic <- crossfit_result
    result_sic$settings$cpn_crit <- "sic"
    estimated_sic <- relieverChangepoint:::.reliever_finalize_result(result_sic)
    expect_equal(estimated_sic$summary$rule, "sic")
  })
})
