test_that("twostep and local_refine run on the null reg_fun path", {
  with_test_timeout({
    data <- matrix(seq_len(24), ncol = 1)

    step <- relieverChangepoint::twostep(
      data = data,
      cpn_max = 1L,
      dm = 5L,
      reg_fun = reg_null,
      method = "BS",
      M = 0L,
      echo = FALSE,
      num_init = 1L
    )

    expect_equal(nrow(step$run_meta), 1L)
    expect_equal(step$cpd_path$select_by, "K")
    expect_equal(cpd_path_candidate(step, K = 1L), 5L)
    expect_equal(length(step$timing$n_model_fit), 1L)
    expect_true(is.finite(step$timing$total_time))
    expect_equal(step$settings$cpn_max, 1L)
    expect_equal(step$settings$init_cand, 0.5)
    expect_false(any(c("M", "wbs_seed") %in% names(step$settings)))

    refined <- relieverChangepoint::local_refine(
      data = data,
      tau_init = 12L,
      dm = 4L,
      r = 0.2,
      reg_fun = reg_null
    )

    expect_equal(
      refined$run_meta,
      data.frame(run_id = 1L, loss_output_id = 1L)
    )
    expect_equal(refined$tau_est, matrix(5, nrow = 1L, ncol = 1L))
    expect_true(is.finite(refined$total_time))
  })
})

test_that("local_refine preserves custom loss-output metadata", {
  with_test_timeout({
    reg_two_outputs <- function(data, l, r, l_end = l, r_end = r,
                                save_model = FALSE,
                                is_virtual_run = FALSE, ...) {
      if (is_virtual_run) {
        return(list(
          n_loss_outputs = 2L,
          loss_output_meta = data.frame(
            row_type = c("first", "second"),
            hyper_value = c(1, 2)
          )
        ))
      }
      y <- data[l_end:r_end, 1L]
      center <- mean(data[l:r, 1L])
      list(loss = cbind(abs(y - center), (y - center)^2))
    }

    refined <- relieverChangepoint::local_refine(
      matrix(seq_len(24), ncol = 1L),
      tau_init = 12L, dm = 4L, r = 0.2,
      reg_fun = reg_two_outputs
    )

    expect_equal(nrow(refined$tau_est), 2L)
    expect_equal(refined$run_meta$row_type, c("first", "second"))
    expect_equal(refined$run_meta$hyper_value, c(1, 2))
  })
})

test_that("twostep supports BS, WBS-family, and SeedBS detail paths", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1)
    method_specs <- list(
      BS = list(M = 0L),
      WBS = list(M = 8L, wbs_seed = 123L),
      WBS_recursive = list(M = 8L, wbs_seed = 123L),
      SeedBS = list(M = 8L)
    )

    for (method in names(method_specs)) {
      res <- do.call(
        relieverChangepoint::twostep,
        c(
          list(
            data = data,
            cpn_max = 1L,
            dm = 5L,
            reg_fun = reg_null,
            method = method,
            num_init = 1L,
            detail = TRUE,
            echo = FALSE
          ),
          method_specs[[method]]
        )
      )

      expect_equal(nrow(res$run_meta), 1L)
      expect_equal(res$cpd_path$select_by, "K")
      expect_equal(nrow(res$cpd_path$candidates), 2L)
      expect_true(is.matrix(res$diagnostics$gain_mat))
      expect_true(is.matrix(res$diagnostics$split_mat))
      expect_equal(res$settings$cpn_max, 1L)
      expect_equal(res$settings$init_cand, 0.5)
      if (method %in% c("WBS", "WBS_recursive")) {
        expect_equal(res$settings$M, 8L)
        expect_equal(res$settings$wbs_seed, 123L)
      } else {
        expect_false(any(c("M", "wbs_seed") %in% names(res$settings)))
      }
      search_intervals <- if (method == "SeedBS") {
        relieverChangepoint::create_seed_itv(30L, 10L)
      } else {
        relieverChangepoint::create_wbs_itv(
          30L, 10L, method_specs[[method]]$M,
          method_specs[[method]]$wbs_seed
        )
      }
      expect_equal(ncol(res$diagnostics$gain_mat),
                   nrow(search_intervals) + 1L)
      expect_true(all(is.finite(res$timing$total_time)))
    }
  })
})

test_that("twostep WBS-family tie handling matches stable snapshots", {
  with_test_timeout({
    data <- matrix(seq_len(30), ncol = 1)
    snapshots <- list(
      BS = matrix(c(5, 5, NaN, 10), nrow = 2L, ncol = 2L),
      WBS = matrix(c(19, 9, NaN, 19), nrow = 2L, ncol = 2L),
      SeedBS = matrix(c(5, 5, NaN, 15), nrow = 2L, ncol = 2L)
    )
    method_specs <- list(
      BS = list(M = 0L),
      WBS = list(M = 8L, wbs_seed = 123L),
      SeedBS = list(M = 8L)
    )

    for (method in names(method_specs)) {
      res <- do.call(
        relieverChangepoint::twostep,
        c(
          list(
            data = data,
            cpn_max = 2L,
            dm = 5L,
            reg_fun = reg_null,
            method = method,
            num_init = 1L,
            detail = TRUE,
            echo = FALSE
          ),
          method_specs[[method]]
        )
      )

      expect_equal(cpd_path_matrix(res), snapshots[[method]],
                   info = method)
      expect_equal(sum(res$timing$n_model_fit),
                   c(BS = 6, WBS = 24, SeedBS = 15)[[method]],
                   info = method)
    }
  })
})

test_that("twostep and local_refine validate user-facing controls", {
  with_test_timeout({
    data <- matrix(seq_len(24), ncol = 1)

    expect_error(
      relieverChangepoint::twostep(data, cpn_max = 1L, dm = 5L),
      "reg_fun must be supplied"
    )
    expect_error(
      relieverChangepoint::local_refine(data, tau_init = 12L, dm = 4L),
      "reg_fun must be supplied"
    )
    expect_error(
      relieverChangepoint::twostep(data, cpn_max = 1L, dm = 5L, reg_fun = reg_null,
                    init_cand = c(0.25, 1.2), echo = FALSE),
      "init_cand"
    )
    expect_error(
      relieverChangepoint::twostep(data, cpn_max = 1L, dm = 5L, reg_fun = reg_null,
                    method = "BS", wbs_stop_crit = NA_real_, echo = FALSE),
      "wbs_stop_crit"
    )
    expect_error(
      relieverChangepoint::local_refine(data, tau_init = 0L, dm = 4L,
                         reg_fun = reg_null),
      "tau_init"
    )
    expect_error(
      relieverChangepoint::local_refine(data, tau_init = 12.5, dm = 4L,
                         reg_fun = reg_null),
      "tau_init"
    )
    expect_error(
      relieverChangepoint::local_refine(data, tau_init = 12L, dm = 4L, r = 1.5,
                         reg_fun = reg_null),
      "r must"
    )

    refined <- relieverChangepoint::local_refine(
      data,
      tau_init = c(15L, 8L),
      dm = 3L,
      r = 0.2,
      reg_fun = reg_null,
      echo = FALSE
    )
    expect_equal(dim(refined$tau_est), c(1L, 2L))
    expect_true(all(is.finite(refined$tau_est)))

    full_window <- relieverChangepoint::local_refine(
      data, tau_init = 12L, dm = 4L, r = 0, reg_fun = reg_null
    )
    unchanged <- relieverChangepoint::local_refine(
      data, tau_init = c(15L, 8L), dm = 4L, r = 1, reg_fun = reg_null
    )
    expect_equal(full_window$tau_est, matrix(4L, 1L, 1L))
    expect_equal(unchanged$tau_est, matrix(c(8L, 15L), 1L, 2L))
  })
})

test_that("twostep initial fits are disjoint and respect dm", {
  with_test_timeout({
    calls <- list()
    recording_loss <- function(data, l, r, l_end = l, r_end = r,
                               save_model = FALSE,
                               is_virtual_run = FALSE) {
      if (is_virtual_run) {
        return(1L)
      }
      calls[[length(calls) + 1L]] <<- c(l = l, r = r)
      list(loss = matrix(0, r_end - l_end + 1L, 1L))
    }
    relieverChangepoint::twostep(
      matrix(seq_len(20L), ncol = 1L),
      reg_fun = recording_loss,
      cpn_max = 1L,
      dm = 5L,
      init_cand = 0.01,
      method = "BS"
    )

    expect_gte(calls[[1L]][["r"]] - calls[[1L]][["l"]] + 1L, 5L)
    expect_gte(calls[[2L]][["r"]] - calls[[2L]][["l"]] + 1L, 5L)
    expect_equal(calls[[1L]][["r"]] + 1L, calls[[2L]][["l"]])
  })
})

test_that("twostep supports WBS stopping criteria through cpd_path", {
  with_test_timeout({
    data <- matrix(c(rep(0, 15L), rep(3, 15L), rep(-2, 15L)), ncol = 1L)
    base <- relieverChangepoint::twostep(
      data = data,
      cpn_max = 2L,
      dm = 5L,
      reg_fun = reg_null,
      method = "WBS",
      M = 8L,
      wbs_seed = 123L,
      num_init = 1L,
      detail = TRUE,
      echo = FALSE
    )
    gains <- base$diagnostics$path_score[[1L]]
    stop_crit <- c(max(gains) + 1, min(gains) - 1)
    stopped <- relieverChangepoint::twostep(
      data = data,
      cpn_max = 2L,
      dm = 5L,
      reg_fun = reg_null,
      method = "WBS",
      M = 8L,
      wbs_seed = 123L,
      wbs_stop_crit = stop_crit,
      num_init = 1L,
      echo = FALSE
    )

    expect_equal(stopped$cpd_path$select_by, "wbs_stop_crit")
    expect_equal(stopped$cpd_path$selector_map$select_value, c(Inf, stop_crit))
    expect_equal(stopped$cpd_path$selector_map$candidate_id, c(1L, 1L, 3L))
    expect_equal(stopped$cpd_path$candidates$K, 0:2)
    expect_true(all(is.na(stopped$cpd_path$candidates$loss)))
    expect_equal(stopped$summary$wbs_stop_crit, stop_crit)
    expect_equal(stopped$summary$K_hat, c(0L, 2L))
    expect_null(stopped$diagnostics)
  })
})
