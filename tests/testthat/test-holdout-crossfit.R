test_that("lasso-crossfit holdout tuning uses training rows only", {
  with_test_timeout({
    set.seed(20260721)
    n <- 24L
    p <- 3L
    x <- matrix(stats::rnorm(n * p), nrow = n)
    beta <- c(1.2, -0.7, 0.4)
    y <- drop(x %*% beta + stats::rnorm(n, sd = 0.25))
    eval_x <- x + matrix(
      stats::rnorm(n * p, sd = 0.2), nrow = n
    )
    eval_y <- drop(
      eval_x %*% beta + stats::rnorm(n, sd = 0.25)
    )
    lam_set <- c(1.5, 0.35)

    fit <- relieverChangepoint::reliever(
      X = x, y = y,
      cpd_family = "lasso_crossfit",
      cpn_max = 1L,
      dm = 6L,
      cov_rate = 0.85,
      method = "SN",
      lam_set = lam_set,
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      nfolds = 2L,
      fold_type = "blk"
    )
    evaluated <- relieverChangepoint::evaluate_reliever_segments(
      fit,
      data = x,
      y = y,
      eval_data = eval_x,
      eval_y = eval_y
    )

    train_data <- cbind(y, x)
    holdout_data <- cbind(eval_y, eval_x)
    training_fit <- relieverChangepoint::reg_fun_lasso_crossfit(
      train_data,
      l = 1L,
      r = n,
      l_end = 1L,
      r_end = n,
      lam_set = lam_set,
      loss_output_types = c("recv", "crossfit_homo_hyper"),
      nfolds = 2L,
      fold_type = "blk"
    )
    fixed_cf <- training_fit$loss[
      , 1L + seq_along(lam_set), drop = FALSE
    ]
    best_lambda_id <- which.min(colSums(fixed_cf))[1L]

    external_fixed <- relieverChangepoint::reg_fun_lasso_solpath(
      rbind(train_data, holdout_data),
      l = 1L,
      r = n,
      l_end = n + 1L,
      r_end = 2L * n,
      lam_set = lam_set
    )$loss
    expected_by_output <- unname(c(
      sum(external_fixed[, best_lambda_id]),
      sum(external_fixed[, best_lambda_id]),
      colSums(external_fixed)
    ))

    no_change <- evaluated$candidates[
      evaluated$candidates$K == 0L, , drop = FALSE
    ]
    output_id <- fit$run_meta$loss_output_id[
      match(no_change$run_id, fit$run_meta$run_id)
    ]
    expect_equal(
      unname(no_change$eval_loss),
      expected_by_output[output_id],
      tolerance = 1e-10
    )

    response_first <- relieverChangepoint::evaluate_reliever_segments(
      fit,
      data = train_data,
      eval_data = holdout_data
    )
    expect_equal(
      response_first$candidates$eval_loss,
      evaluated$candidates$eval_loss,
      tolerance = 1e-10
    )

    recv_run <- fit$run_meta$run_id[fit$run_meta$row_type == "recv"]
    selected <- relieverChangepoint::select_holdout(
      fit,
      data = x,
      y = y,
      eval_data = eval_x,
      eval_y = eval_y,
      K = 0L,
      run_ids = recv_run
    )
    expect_s3_class(selected, "reliever_model_selection")
    expect_identical(selected$K_hat, 0L)
    expect_identical(selected$row_type, "recv")
    expect_equal(
      selected$score,
      no_change$eval_loss[no_change$run_id == recv_run],
      tolerance = 1e-10
    )

    selected_by_type <- relieverChangepoint::select_holdout(
      fit,
      data = x,
      y = y,
      eval_data = eval_x,
      eval_y = eval_y,
      K = 0L,
      run_type = "recv"
    )
    expect_equal(selected_by_type, selected)

    selected_fixed <- relieverChangepoint::select_holdout(
      fit,
      data = x,
      y = y,
      eval_data = eval_x,
      eval_y = eval_y,
      K = 0L,
      run_type = "crossfit_homo_hyper"
    )
    best_fixed_id <- which.min(colSums(external_fixed))[1L]
    expect_identical(selected_fixed$row_type, "crossfit_homo_hyper")
    expect_equal(selected_fixed$hyper_value, lam_set[best_fixed_id])
    expect_equal(
      selected_fixed$score,
      sum(external_fixed[, best_fixed_id]),
      tolerance = 1e-10
    )

    expect_error(
      relieverChangepoint::select_holdout(
        fit,
        data = x,
        y = y,
        eval_data = eval_x,
        eval_y = eval_y,
        run_ids = recv_run,
        run_type = "recv"
      ),
      "only one of run_ids and run_type"
    )
  }, seconds = 15)
})

test_that("lasso holdout validates response domains from the fitted family", {
  with_test_timeout({
    set.seed(20260724)
    n <- 30L
    x <- matrix(stats::rnorm(n * 3L), nrow = n)
    lam_set <- c(0.8, 0.2)
    y_binomial <- rep(c(0, 1), length.out = n)

    binomial_fit <- relieverChangepoint::reliever(
      X = x, y = y_binomial,
      cpd_family = "lasso",
      cpn_max = 0L, dm = 8L, cov_rate = 1,
      cache_backend = "by_cost_mat",
      family = "binomial",
      lam_set = lam_set
    )
    expect_error(
      relieverChangepoint::select_holdout(
        binomial_fit,
        data = x,
        y = y_binomial,
        eval_data = x,
        eval_y = rep(2, n)
      ),
      "eval_y must contain only 0 and 1"
    )
    expect_error(
      relieverChangepoint::select_holdout(
        binomial_fit,
        data = cbind(y_binomial, x),
        eval_data = cbind(rep(-1, n), x)
      ),
      "response column of eval_data must contain only 0 and 1"
    )

    y_poisson <- stats::rpois(n, lambda = 2)
    poisson_fit <- relieverChangepoint::reliever(
      X = x, y = y_poisson,
      cpd_family = "lasso_crossfit",
      cpn_max = 0L, dm = 8L, cov_rate = 1,
      cache_backend = "by_cost_mat",
      family = "poisson",
      lam_set = lam_set,
      nfolds = 2L,
      fold_type = "blk"
    )
    fractional_poisson <- relieverChangepoint::select_holdout(
      poisson_fit,
      data = x,
      y = y_poisson,
      eval_data = x,
      eval_y = y_poisson + 0.25,
      run_type = "recv"
    )
    expect_s3_class(fractional_poisson, "reliever_model_selection")
    expect_error(
      relieverChangepoint::select_holdout(
        poisson_fit,
        data = x,
        y = y_poisson,
        eval_data = x,
        eval_y = rep(-1, n),
        run_type = "recv"
      ),
      "eval_y must contain non-negative values"
    )
  }, seconds = 15)
})

test_that("classifier crossfit stops before external rows can enter fitting", {
  with_test_timeout({
    set.seed(20260722)
    n <- 16L
    x <- matrix(stats::rnorm(n * 2L), nrow = n)
    holdout <- matrix(1e6, nrow = n, ncol = 2L)
    tracker <- new.env(parent = emptyenv())
    tracker$called <- FALSE

    prob_fun <- function(data, y, train_id, eval_id, hyper, ...) {
      tracker$called <- TRUE
      rep(0.75, length(eval_id))
    }
    classifier_reg_fun <- function(data, l, r, l_end = l, r_end = r,
                                   save_model = FALSE,
                                   is_virtual_run = FALSE) {
      relieverChangepoint::reg_fun_clf_crossfit_template(
        data = data,
        l = l,
        r = r,
        l_end = l_end,
        r_end = r_end,
        prob_fun = prob_fun,
        hyper_set = list(list()),
        loss_output_types = c(
          "recv", "incv", "crossfit_homo_hyper"
        ),
        nfolds = 2L,
        fold_type = "op",
        save_model = save_model,
        is_virtual_run = is_virtual_run
      )
    }
    fit_reg_fun <- function(data, l, r, l_end = l, r_end = r,
                            save_model = FALSE,
                            is_virtual_run = FALSE) {
      if (is_virtual_run) {
        return(classifier_reg_fun(
          data, l, r, l_end, r_end,
          save_model = save_model, is_virtual_run = TRUE
        ))
      }
      list(
        loss = matrix(0, nrow = r_end - l_end + 1L, ncol = 3L),
        model = NULL
      )
    }

    fit <- relieverChangepoint::reliever_generic(
      x,
      reg_fun = fit_reg_fun,
      cpn_max = 0L,
      dm = 4L,
      cov_rate = 0.85,
      method = "SN"
    )
    tracker$called <- FALSE

    expect_error(
      relieverChangepoint::evaluate_reliever_segments(
        fit,
        data = x,
        eval_data = holdout,
        reg_fun = classifier_reg_fun
      ),
      "No evaluation rows were used for fitting",
      fixed = TRUE
    )
    expect_false(tracker$called)
  }, seconds = 5)
})
