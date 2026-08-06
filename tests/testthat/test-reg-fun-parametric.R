test_that("Gaussian variance losses match explicit MLE likelihoods", {
  set.seed(20260720)
  x <- matrix(rnorm(60L), 30L, 2L)
  fit_id <- 3:20
  eval_id <- 21:28
  mu <- c(0.2, -0.3)

  var_fit <- relieverChangepoint::reg_fun_var(
    x, min(fit_id), max(fit_id), min(eval_id), max(eval_id),
    save_model = TRUE, mu = mu
  )
  centered_fit <- sweep(x[fit_id, , drop = FALSE], 2L, mu, "-")
  covariance <- crossprod(centered_fit) / length(fit_id)
  centered_eval <- sweep(x[eval_id, , drop = FALSE], 2L, mu, "-")
  expected <- ncol(x) * log(2 * pi) +
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
    rowSums((centered_eval %*% solve(covariance)) * centered_eval)

  expect_equal(as.numeric(var_fit$loss), expected, tolerance = 1e-10)
  expect_equal(var_fit$model$center, mu)
  expect_equal(var_fit$model$covariance, covariance)
  expect_equal(
    relieverChangepoint::reg_fun_var(
      x, min(fit_id), max(fit_id), save_model = TRUE, mu = 0
    )$model$center,
    rep(0, ncol(x))
  )

  meanvar_fit <- relieverChangepoint::reg_fun_meanvar(
    x, min(fit_id), max(fit_id), min(eval_id), max(eval_id),
    save_model = TRUE
  )
  fitted_mu <- colMeans(x[fit_id, , drop = FALSE])
  centered_fit <- sweep(x[fit_id, , drop = FALSE], 2L, fitted_mu, "-")
  expected_covariance <- crossprod(centered_fit) / length(fit_id)
  expect_equal(meanvar_fit$model$center, fitted_mu, tolerance = 1e-12)
  expect_equal(
    meanvar_fit$model$covariance, expected_covariance, tolerance = 1e-12
  )
})

test_that("LM and GLM losses match stats fits", {
  set.seed(20260721)
  n <- 50L
  x <- matrix(rnorm(n * 2L), n, 2L)
  y <- 0.5 + x[, 1L] - 0.7 * x[, 2L] + rnorm(n, sd = 0.4)
  fit_id <- 1:30
  eval_id <- 31:50
  design <- cbind(1, x)

  lm_loss <- relieverChangepoint::reg_fun_lm(
    cbind(y, x), min(fit_id), max(fit_id), min(eval_id), max(eval_id),
    save_model = TRUE
  )
  lm_fit <- stats::lm.fit(design[fit_id, , drop = FALSE], y[fit_id])
  expected <- (
    y[eval_id] -
      drop(design[eval_id, , drop = FALSE] %*% lm_fit$coefficients)
  )^2
  expect_equal(as.numeric(lm_loss$loss), expected, tolerance = 1e-12)

  trials <- rep(5L, n)
  probability <- stats::plogis(-0.2 + 0.8 * x[, 1L])
  success <- stats::rbinom(n, trials, probability)
  response <- cbind(success, trials - success)
  glm_loss <- relieverChangepoint::reg_fun_glm(
    cbind(response, x),
    min(fit_id), max(fit_id), min(eval_id), max(eval_id),
    family = stats::binomial(), response_ncol = 2L,
    save_model = TRUE
  )
  glm_fit <- suppressWarnings(stats::glm.fit(
    design[fit_id, , drop = FALSE],
    success[fit_id] / trials[fit_id],
    weights = trials[fit_id],
    family = stats::binomial()
  ))
  fitted <- stats::binomial()$linkinv(
    drop(design[eval_id, , drop = FALSE] %*% glm_fit$coefficients)
  )
  expected <- stats::binomial()$dev.resids(
    success[eval_id] / trials[eval_id],
    fitted,
    trials[eval_id]
  )
  expect_equal(as.numeric(glm_loss$loss), expected, tolerance = 1e-10)
  expect_identical(glm_loss$model$response_ncol, 2L)
})

test_that("built-in searches validate once while direct losses stay local", {
  make_loss_fun <- getFromNamespace(
    ".reliever_make_individual_loss_fun", "relieverChangepoint"
  )

  x <- seq_len(8L)
  lm_data <- cbind(y = 1 + 2 * x, x = x)
  lm_data[8L, 2L] <- NA_real_
  direct_lm <- relieverChangepoint::reg_fun_lm(
    data = lm_data, l = 1L, r = 4L, l_end = 5L, r_end = 6L
  )
  expect_true(all(is.finite(direct_lm$loss)))
  expect_error(
    make_loss_fun(
      data = lm_data,
      reg_fun = relieverChangepoint::reg_fun_lm,
      para_list = list(),
      n_loss_outputs = 1L
    ),
    "finite numeric observations"
  )

  expect_false(exists(
    ".reliever_parametric_wrapper",
    envir = asNamespace("relieverChangepoint"),
    inherits = FALSE
  ))
})

test_that("exponential-family losses use interval MLE parameters", {
  set.seed(20260722)
  x <- stats::rpois(40L, lambda = 3)
  fit_id <- 1:25
  eval_id <- 26:40
  out <- relieverChangepoint::reg_fun_em(
    x, min(fit_id), max(fit_id), min(eval_id), max(eval_id),
    family = "pois", save_model = TRUE
  )
  lambda <- mean(x[fit_id])

  expect_equal(out$model$param, lambda)
  expect_equal(
    as.numeric(out$loss),
    -2 * stats::dpois(x[eval_id], lambda = lambda, log = TRUE)
  )
  expect_equal(
    relieverChangepoint::reg_fun_em(
      x, 1L, 1L, family = "poisson", is_virtual_run = TRUE
    )$loss_output_meta$row_type,
    "em_pois"
  )
  expect_error(
    relieverChangepoint::reg_fun_em(
      c(0.2, 1.2), 1L, 2L, family = "beta"
    ),
    "strictly between"
  )
})

test_that("cv.reliever dispatches all single-path parametric families", {
  with_test_timeout({
    set.seed(20260723)
    n <- 36L
    observations <- matrix(rnorm(n * 2L), n, 2L)
    predictors <- matrix(rnorm(n * 2L), n, 2L)
    response <- 1 + predictors[, 1L] + rnorm(n)
    binary <- stats::rbinom(
      n, 1L, stats::plogis(0.3 + predictors[, 1L])
    )
    counts <- stats::rpois(n, 2.5)
    common <- list(
      cpn_max = 0L, dm = 8L, cov_rate = 0.6,
      method = "SN", nfolds = 2L
    )
    calls <- list(
      var = c(list(X = observations, cpd_family = "var"), common),
      meanvar = c(list(X = observations, cpd_family = "meanvar"), common),
      lm = c(
        list(X = predictors, y = response, cpd_family = "lm"),
        common
      ),
      glm = c(
        list(
          X = predictors, y = binary, cpd_family = "glm",
          family = stats::binomial()
        ),
        common
      ),
      em = c(
        list(X = counts, cpd_family = "em", family = "pois"),
        common
      )
    )

    for (family in names(calls)) {
      fit <- do.call(relieverChangepoint::cv.reliever, calls[[family]])
      expect_s3_class(fit, "cv_reliever_result")
      expect_identical(fit$settings$cpd_family, family, info = family)
      expect_identical(fit$summary$K_hat, 0L, info = family)
      expect_length(fit$summary$cpd_hat[[1L]], 0L)
      expect_equal(nrow(fit$full_data_fit$run_meta), 1L, info = family)
      expect_true(is.finite(fit$summary$cv_mean), info = family)
    }

    trials <- rep(4L, n)
    success <- stats::rbinom(n, trials, stats::plogis(predictors[, 1L]))
    binomial_counts <- relieverChangepoint::cv.reliever(
      X = predictors,
      y = cbind(success, trials - success),
      cpd_family = "glm",
      family = stats::binomial(),
      cpn_max = 0L, dm = 8L, cov_rate = 0.6, nfolds = 2L
    )
    expect_identical(
      binomial_counts$settings$input_spec$n_response_columns, 2L
    )
    expect_true(is.finite(binomial_counts$summary$cv_mean))

    expect_warning(
      combined_binomial <- relieverChangepoint::cv.reliever(
        X = cbind(success, trials - success, predictors),
        cpd_family = "glm",
        family = stats::binomial(),
        response_ncol = 2L,
        cpn_max = 0L, dm = 8L, cov_rate = 0.6, nfolds = 2L
      ),
      "expects a response y"
    )
    expect_identical(
      combined_binomial$settings$input_spec$n_response_columns, 2L
    )
    expect_identical(
      combined_binomial$settings$input_spec$n_predictors, 2L
    )
  })
})

test_that("parametric wrappers match the primary reliever dispatcher", {
  with_test_timeout({
    set.seed(20260724)
    n <- 30L
    observations <- matrix(rnorm(n * 2L), n, 2L)
    predictors <- matrix(rnorm(n * 2L), n, 2L)
    response <- 1 + predictors[, 1L] + rnorm(n)
    binary <- rbinom(n, 1L, plogis(predictors[, 1L]))
    counts <- rpois(n, 2)
    common <- list(
      cpn_max = 0L, dm = 7L, cov_rate = 0.6, method = "SN"
    )
    calls <- list(
      var = list(
        wrapper = c(list(data = observations), common),
        dispatched = c(
          list(X = observations, cpd_family = "var"), common
        )
      ),
      meanvar = list(
        wrapper = c(list(data = observations), common),
        dispatched = c(
          list(X = observations, cpd_family = "meanvar"), common
        )
      ),
      lm = list(
        wrapper = c(list(data = cbind(response, predictors)), common),
        dispatched = c(
          list(X = predictors, y = response, cpd_family = "lm"), common
        )
      ),
      glm = list(
        wrapper = c(
          list(data = cbind(binary, predictors), family = binomial()),
          common
        ),
        dispatched = c(
          list(
            X = predictors, y = binary, cpd_family = "glm",
            family = binomial()
          ),
          common
        )
      ),
      em = list(
        wrapper = c(list(data = counts, family = "pois"), common),
        dispatched = c(
          list(X = counts, cpd_family = "em", family = "pois"), common
        )
      )
    )
    wrappers <- list(
      var = relieverChangepoint::reliever_var,
      meanvar = relieverChangepoint::reliever_meanvar,
      lm = relieverChangepoint::reliever_lm,
      glm = relieverChangepoint::reliever_glm,
      em = relieverChangepoint::reliever_em
    )

    for (family in names(calls)) {
      wrapper_fit <- do.call(wrappers[[family]], calls[[family]]$wrapper)
      dispatched_fit <- do.call(
        relieverChangepoint::reliever, calls[[family]]$dispatched
      )
      expect_equal(
        wrapper_fit$cpd_path$candidates,
        dispatched_fit$cpd_path$candidates,
        info = family
      )
      expect_equal(wrapper_fit$run_meta, dispatched_fit$run_meta, info = family)
      expect_identical(
        dispatched_fit$settings$cpd_family, family, info = family
      )
    }
  })
})

test_that("parametric cv families validate their required arguments", {
  x <- matrix(seq_len(40L), 20L, 2L)
  expect_error(
    relieverChangepoint::cv.reliever(
      x[, 2L, drop = FALSE], y = x[, 1L], cpd_family = "glm",
      cpn_max = 0L, dm = 5L, nfolds = 2L
    ),
    "family is required"
  )
  expect_error(
    relieverChangepoint::cv.reliever(
      x[, 1L], cpd_family = "em",
      cpn_max = 0L, dm = 5L, nfolds = 2L
    ),
    "family is required"
  )
  expect_error(
    relieverChangepoint::reliever_glm(
      cbind(x[, 1L], x[, 2L]), cpn_max = 0L, dm = 5L
    ),
    "family is required"
  )
  expect_error(
    relieverChangepoint::reliever_em(
      x[, 1L], cpn_max = 0L, dm = 5L
    ),
    "family is required"
  )
})
