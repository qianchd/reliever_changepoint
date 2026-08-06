test_that("lasso solution-path caches are tied to lam_set", {
  with_test_timeout({
    set.seed(2026)
    x <- matrix(stats::rnorm(90), 30L, 3L)
    data <- cbind(x[, 1L] - x[, 2L] + stats::rnorm(30L), x)

    for (backend in c("by_cost_mat", "by_loss_block")) {
      fit <- relieverChangepoint::reliever_generic(
        data = data,
        reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
        cpn_max = 1L, dm = 5L, cov_rate = 0.8,
        detail = TRUE, cache_backend = backend,
        run_cpd_ids = 1:2,
        lam_set = c(2, 1)
      )
      expect_error(
        relieverChangepoint::reliever_generic(
          data = data,
          reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
          cpn_max = 1L, dm = 5L, cov_rate = 0.8,
          detail = TRUE, cache_backend = backend,
          cache_profile = fit$cache_profile,
          run_cpd_ids = 1:2,
          lam_set = c(1, 0.5)
        ),
        "different data, reg_fun, dc_grid, or loss-function arguments",
        info = backend
      )
    }
  })
})

test_that("lasso interval loss explains glmnet's predictor requirement", {
  with_test_timeout({
    set.seed(2026)
    x <- matrix(stats::rnorm(40L), ncol = 1L)
    data <- cbind(2 * x[, 1L] + stats::rnorm(40L), x)
    expect_error(
      relieverChangepoint::reg_fun_lasso_solpath(
        data, l = 1L, r = 30L, l_end = 1L, r_end = 40L,
        lam_set = c(2, 1)
      ),
      "at least two observations and two predictors"
    )
  })
})

test_that("direct lasso validation is limited to involved rows", {
  set.seed(20260721)
  n <- 12L
  data <- cbind(
    response = stats::rnorm(n),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
  data[n, 1L] <- NA_real_

  direct <- relieverChangepoint::reg_fun_lasso_solpath(
    data = data, l = 1L, r = 8L, l_end = 1L, r_end = 8L,
    lam_set = 1, family = "gaussian"
  )
  expect_true(all(is.finite(direct$loss)))

  make_loss_fun <- getFromNamespace(
    ".reliever_make_individual_loss_fun", "relieverChangepoint"
  )
  expect_error(
    make_loss_fun(
      data = data,
      reg_fun = relieverChangepoint::reg_fun_lasso_solpath,
      para_list = list(lam_set = 1, family = "gaussian"),
      n_loss_outputs = 1L
    ),
    "finite"
  )
})
