test_that("reg_fun_lasso_solpath returns one loss column per lambda", {
  with_test_timeout({
    set.seed(20240524)
    x <- matrix(stats::rnorm(30L * 4L), nrow = 30L)
    y <- drop(x %*% c(1, -0.5, 0, 0))
    data <- cbind(y, x)
    lam_set <- c(0.5, 0.05)

    virtual <- relieverChangepoint::reg_fun_lasso_solpath(
      data, 1L, 1L, is_virtual_run = TRUE, lam_set = lam_set
    )
    expect_equal(virtual$n_loss_outputs, length(lam_set))
    expect_equal(virtual$loss_output_meta$row_type, rep("lasso", length(lam_set)))
    expect_equal(virtual$loss_output_meta$hyper_id, seq_along(lam_set))
    expect_equal(unlist(virtual$loss_output_meta$hyper_value), lam_set)
    expect_equal(virtual$loss_output_meta$hyper_name, rep("lambda", 2L))
    expect_equal(virtual$loss_output_meta$loss_kind, rep("rss", 2L))

    out <- relieverChangepoint::reg_fun_lasso_solpath(
      data,
      l = 1L,
      r = 20L,
      l_end = 1L,
      r_end = 25L,
      lam_set = lam_set
    )

    expect_equal(dim(out$loss), c(25L, length(lam_set)))
    expect_true(all(is.finite(out$loss)))
    expect_null(out$model)

    reference <- glmnet::glmnet(
      x = x[1:20, , drop = FALSE],
      y = y[1:20],
      family = "gaussian",
      lambda = lam_set / sqrt(20),
      intercept = FALSE,
      standardize = FALSE,
      control = list(thresh = 1e-7)
    )
    expected <- sweep(
      stats::predict(reference, newx = x[1:25, , drop = FALSE]),
      1L, y[1:25], "-"
    )^2
    expect_equal(unname(out$loss), unname(expected), tolerance = 1e-10)

    expect_error(
      relieverChangepoint::reg_fun_lasso_solpath(
        data, 1L, 20L, lam_set = lam_set, family = "cox"
      ),
      'family must be one of "gaussian", "binomial", or "poisson"'
    )
  })
})

test_that("lasso losses validate non-Gaussian response domains", {
  with_test_timeout({
    set.seed(20260723)
    n <- 60L
    x <- matrix(stats::rnorm(n * 3L), nrow = n)
    lam_set <- c(0.8, 0.2)

    binomial_data <- cbind(rep(c(0, 1), length.out = n), x)
    binomial_virtual <- relieverChangepoint::reg_fun_lasso_solpath(
      binomial_data, 1L, 1L,
      is_virtual_run = TRUE,
      lam_set = lam_set,
      family = "binomial"
    )
    expect_equal(
      binomial_virtual$loss_output_meta$loss_kind,
      rep("negative_log_likelihood", length(lam_set))
    )
    binomial_loss <- relieverChangepoint::reg_fun_lasso_solpath(
      binomial_data, 1L, 40L, 1L, n,
      lam_set = lam_set,
      family = "binomial"
    )$loss
    expect_equal(dim(binomial_loss), c(n, length(lam_set)))
    expect_true(all(is.finite(binomial_loss)))

    invalid_binomial <- binomial_data
    invalid_binomial[1L, 1L] <- 2
    expect_error(
      relieverChangepoint::reg_fun_lasso_solpath(
        invalid_binomial, 1L, 40L,
        lam_set = lam_set,
        family = "binomial"
      ),
      "only 0 and 1"
    )

    poisson_data <- cbind(stats::rpois(n, lambda = 2), x)
    poisson_loss <- relieverChangepoint::reg_fun_lasso_solpath(
      poisson_data, 1L, 40L, 1L, n,
      lam_set = lam_set,
      family = "poisson"
    )$loss
    expect_equal(dim(poisson_loss), c(n, length(lam_set)))
    expect_true(all(is.finite(poisson_loss)))

    fractional_poisson <- poisson_data
    fractional_poisson[1L, 1L] <- 0.5
    fractional_loss <- relieverChangepoint::reg_fun_lasso_solpath(
      fractional_poisson, 1L, 40L, 1L, n,
      lam_set = lam_set,
      family = "poisson"
    )
    expect_true(all(is.finite(fractional_loss$loss)))

    negative_poisson <- poisson_data
    negative_poisson[1L, 1L] <- -0.5
    expect_error(
      relieverChangepoint::reg_fun_lasso_crossfit(
        negative_poisson, 1L, 40L,
        is_virtual_run = TRUE,
        lam_set = lam_set,
        family = "poisson"
      ),
      "non-negative"
    )
  })
})

test_that("reg_fun_lasso_solpath evaluates both sides of its fitted interval", {
  with_test_timeout({
    set.seed(20240525)
    x <- matrix(stats::rnorm(32L * 4L), nrow = 32L)
    y <- drop(x %*% c(1, -0.5, 0, 0))
    data <- cbind(y, x)
    lam_set <- c(0.4, 0.04)

    out <- relieverChangepoint::reg_fun_lasso_solpath(
      data,
      l = 4L,
      r = 24L,
      l_end = 2L,
      r_end = 28L,
      lam_set = lam_set
    )

    expect_equal(dim(out$loss), c(27L, length(lam_set)))
  })
})

test_that("reg_fun_lasso_crossfit exposes requested ReCV outputs", {
  with_test_timeout({
    set.seed(20240612)
    x <- matrix(stats::rnorm(36L * 5L), nrow = 36L)
    y <- drop(x %*% c(0.8, -0.5, 0.3, 0, 0) + stats::rnorm(36L, sd = 0.2))
    data <- cbind(y, x)
    lam_set <- c(0.3, 0.1, 0.03)
    nfolds <- 3L

    virtual <- relieverChangepoint::reg_fun_lasso_crossfit(
      data, 1L, 1L, is_virtual_run = TRUE, lam_set = lam_set
    )
    expect_equal(virtual$n_loss_outputs, 1L)
    expect_identical(virtual$loss_output_meta$row_type, "recv")

    output_types <- c("crossfit_homo_hyper", "incv", "recv")
    virtual <- relieverChangepoint::reg_fun_lasso_crossfit(
      data, 1L, 1L, is_virtual_run = TRUE, lam_set = lam_set,
      loss_output_types = output_types
    )
    expect_equal(virtual$n_loss_outputs, 2L + length(lam_set))
    expect_equal(
      virtual$loss_output_meta$row_type,
      c(
        "recv", "incv",
        rep("crossfit_homo_hyper", length(lam_set))
      )
    )
    expect_equal(virtual$loss_output_meta$hyper_id, c(NA, NA, 1:3))
    expect_equal(unlist(virtual$loss_output_meta$hyper_value[3:5]), lam_set)
    expect_equal(
      virtual$loss_output_meta$hyper_name,
      c(NA, NA, rep("lambda", length(lam_set)))
    )
    expect_identical(
      virtual$loss_output_meta$default_selection,
      c(TRUE, rep(FALSE, 1L + length(lam_set)))
    )
    expect_equal(
      virtual$loss_output_meta$loss_kind,
      c(
        "crossfit_loss", "rss",
        rep("crossfit_loss", length(lam_set))
      )
    )

    out <- relieverChangepoint::reg_fun_lasso_crossfit(
      data,
      l = 8L,
      r = 28L,
      l_end = 5L,
      r_end = 31L,
      lam_set = lam_set,
      loss_output_types = output_types,
      nfolds = nfolds,
      fold_type = "stable_blk",
      fold_stable_const = 2L
    )

    expect_equal(dim(out$loss), c(27L, 2L + length(lam_set)))
    expect_identical(
      colnames(out$loss),
      c("recv", "incv", rep("crossfit_homo_hyper", length(lam_set)))
    )
    expect_true(all(is.finite(out$loss)))

    expect_error(
      relieverChangepoint::reg_fun_lasso_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE, lam_set = lam_set,
        loss_output_types = "crossfit_homo_hyper"
      ),
      'must include "recv"'
    )
    expect_error(
      relieverChangepoint::reg_fun_lasso_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE, lam_set = lam_set,
        loss_output_types = c("recv", "fixed_hyper_in")
      ),
      "Unsupported loss_output_types"
    )
  })
})

test_that("reg_fun_mean_crossfit returns mean-square ReCV rows", {
  with_test_timeout({
    set.seed(20260613)
    data <- matrix(stats::rnorm(30L * 3L), nrow = 30L)

    virtual <- relieverChangepoint::reg_fun_mean_crossfit(
      data, 1L, 1L, is_virtual_run = TRUE
    )
    expect_equal(virtual$n_loss_outputs, 1L)
    expect_identical(virtual$loss_output_meta$row_type, "recv")

    virtual <- relieverChangepoint::reg_fun_mean_crossfit(
      data, 1L, 1L, is_virtual_run = TRUE,
      loss_output_types = c("recv", "incv")
    )
    expect_equal(virtual$n_loss_outputs, 2L)
    expect_equal(
      virtual$loss_output_meta$row_type,
      c("recv", "incv")
    )
    expect_equal(
      virtual$loss_output_meta$loss_kind,
      c("crossfit_loss", "rss")
    )

    out <- relieverChangepoint::reg_fun_mean_crossfit(
      data, l = 8L, r = 24L, l_end = 5L, r_end = 27L,
      nfolds = 3L, fold_type = "blk",
      loss_output_types = c("recv", "incv")
    )
    mu <- colMeans(data[8:24, , drop = FALSE])
    expected_in <- rowMeans(sweep(data[5:27, , drop = FALSE], 2L, mu, "-")^2)

    expect_equal(dim(out$loss), c(23L, 2L))
    expect_equal(out$loss[, 2L], expected_in, tolerance = 1e-12)

    vector_out <- relieverChangepoint::reg_fun_mean_crossfit(
      data[, 1L], l = 8L, r = 24L, l_end = 5L, r_end = 27L,
      nfolds = 3L, fold_type = "blk"
    )
    matrix_out <- relieverChangepoint::reg_fun_mean_crossfit(
      data[, 1L, drop = FALSE], l = 8L, r = 24L,
      l_end = 5L, r_end = 27L, nfolds = 3L, fold_type = "blk"
    )
    expect_equal(vector_out$loss, matrix_out$loss, tolerance = 1e-12)
  })
})

test_that("crossfit folds allow interval lengths not divisible by nfolds", {
  with_test_timeout({
    set.seed(20260720)
    data <- matrix(stats::rnorm(23L * 2L), nrow = 23L)

    for (fold_type in c("op", "blk", "stable_blk")) {
      out <- relieverChangepoint::reg_fun_mean_crossfit(
        data,
        l = 3L, r = 19L,
        l_end = 1L, r_end = 23L,
        nfolds = 4L,
        fold_type = fold_type,
        fold_stable_const = 1
      )
      expect_equal(dim(out$loss), c(23L, 1L), info = fold_type)
      expect_true(all(is.finite(out$loss)), info = fold_type)
    }
  })
})

test_that("KDE ReCV losses reproduce direct NLL formulas", {
  with_test_timeout({
    set.seed(20260613)
    x <- matrix(stats::rnorm(18L * 2L), nrow = 18L)
    dist_sq <- as.matrix(stats::dist(x))^2
    bandwidth_vec <- c(0.8, 1.3)

    nll <- relieverChangepoint::reg_fun_kde_nll_crossfit(
      dist_sq, l = 5L, r = 14L, l_end = 4L, r_end = 16L,
      bandwidth_vec = bandwidth_vec, var_dim = 2L,
      nfolds = 2L, fold_type = "blk",
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      )
    )
    fixed_cf <- nll$loss[, 2L + seq_along(bandwidth_vec), drop = FALSE]
    best_bandwidth_id <- which.min(colSums(fixed_cf))[1L]
    dist_sub <- dist_sq[4:16, 5:14, drop = FALSE]
    row_min_dist <- apply(dist_sub, 1L, min)
    centered_dist <- sweep(dist_sub, 1L, row_min_dist, "-")
    direct_nll <- sapply(bandwidth_vec, function(bandwidth) {
      two_bandwidth_sq <- 2 * bandwidth^2
      log_fac <- -(2 / 2) * log(pi * two_bandwidth_sq)
      log_kernel_mean <- -row_min_dist / two_bandwidth_sq +
        log(rowMeans(exp(-centered_dist / two_bandwidth_sq)))
      -(log_fac + log_kernel_mean)
    })
    expect_equal(
      unname(nll$loss[, 2L]),
      unname(direct_nll[, best_bandwidth_id]),
      tolerance = 1e-12
    )
  })
})

test_that("KDE NLL default bandwidth grid is scale-adaptive", {
  with_test_timeout({
    set.seed(20260716)
    x <- matrix(stats::rnorm(12L * 2L), nrow = 12L)
    dist_sq <- as.matrix(stats::dist(x))^2
    max_lag <- ceiling(sqrt(nrow(dist_sq)))
    local_dist_sq <- unlist(lapply(seq_len(max_lag), function(lag) {
      n <- nrow(dist_sq)
      dist_sq[cbind(seq_len(n - lag), seq.int(lag + 1L, n))]
    }), use.names = FALSE)
    distance_scale <- stats::median(local_dist_sq)
    default_bandwidth <- exp(seq(
      log(sqrt(distance_scale / 120)),
      log(sqrt(distance_scale)),
      length.out = 20L
    ))
    virtual <- relieverChangepoint::reg_fun_kde_nll_crossfit(
      dist_sq, 1L, 12L,
      var_dim = 2L, is_virtual_run = TRUE,
      loss_output_types = c("recv", "crossfit_homo_hyper")
    )
    fixed_meta <- virtual$loss_output_meta$row_type == "crossfit_homo_hyper"
    observed_bandwidth <- unlist(
      virtual$loss_output_meta$hyper_value[fixed_meta]
    )
    expect_equal(
      observed_bandwidth, default_bandwidth, tolerance = 1e-12
    )
    expect_equal(
      unique(virtual$loss_output_meta$hyper_name[fixed_meta]), "bandwidth"
    )
    expect_equal(
      virtual$loss_output_meta$loss_kind,
      rep("crossfit_loss", 21L)
    )

    scaled <- relieverChangepoint::reg_fun_kde_nll_crossfit(
      16 * dist_sq, 1L, 12L,
      var_dim = 2L, is_virtual_run = TRUE,
      loss_output_types = c("recv", "crossfit_homo_hyper")
    )
    scaled_meta <- scaled$loss_output_meta$row_type == "crossfit_homo_hyper"
    scaled_bandwidth <- unlist(
      scaled$loss_output_meta$hyper_value[scaled_meta]
    )
    expect_equal(
      scaled_bandwidth, 4 * observed_bandwidth, tolerance = 1e-12
    )

    explicit_bandwidth <- c(0.2, 0.7)
    explicit <- relieverChangepoint::reg_fun_kde_nll_crossfit(
      matrix(0, 4L, 4L), 1L, 4L,
      bandwidth_vec = explicit_bandwidth,
      var_dim = 2L, is_virtual_run = TRUE,
      loss_output_types = c("recv", "crossfit_homo_hyper")
    )
    explicit_meta <- explicit$loss_output_meta$row_type == "crossfit_homo_hyper"
    expect_equal(
      unlist(explicit$loss_output_meta$hyper_value[explicit_meta]),
      explicit_bandwidth
    )
    expect_error(
      relieverChangepoint::reg_fun_kde_nll_crossfit(
        matrix(0, 4L, 4L), 1L, 4L,
        var_dim = 2L, is_virtual_run = TRUE
      ),
      "positive finite distances"
    )
  })
})

test_that("automatic lasso and KDE NLL grids bracket cross-fitted minima", {
  with_test_timeout({
    fixed_segment_cf_loss <- function(reg_fun, data, changepoints,
                                      n_hyper, ...) {
      boundaries <- c(0L, changepoints, nrow(data))
      total_loss <- numeric(n_hyper)
      for (segment in seq_len(length(boundaries) - 1L)) {
        l <- boundaries[segment] + 1L
        r <- boundaries[segment + 1L]
        fit <- reg_fun(
          data, l, r, l, r,
          loss_output_types = c("recv", "crossfit_homo_hyper"), ...
        )
        fixed_cf <- fit$loss[, 1L + seq_len(n_hyper), drop = FALSE]
        total_loss <- total_loss + colSums(fixed_cf)
      }
      total_loss
    }

    set.seed(2026)
    n_seg <- 80L
    p <- 20L
    n <- 3L * n_seg
    x <- matrix(stats::rnorm(n * p), n, p)
    beta <- rbind(
      matrix(c(2, -1.5, 1, rep(0, p - 3L)), n_seg, p, byrow = TRUE),
      matrix(c(-2, 1.5, -1, rep(0, p - 3L)), n_seg, p, byrow = TRUE),
      matrix(c(1.5, -1, 2, rep(0, p - 3L)), n_seg, p, byrow = TRUE)
    )
    lasso_data <- cbind(rowSums(x * beta) + stats::rnorm(n), x)
    lam_set <- relieverChangepoint:::.reliever_auto_lam_set(lasso_data)
    expect_true(all(diff(lam_set) < 0))
    lasso_loss <- fixed_segment_cf_loss(
      relieverChangepoint::reg_fun_lasso_crossfit,
      lasso_data, c(n_seg, 2L * n_seg), length(lam_set),
      nfolds = 3L, fold_type = "op", lam_set = lam_set
    )
    lasso_training_loss <- numeric(length(lam_set))
    boundaries <- c(0L, n_seg, 2L * n_seg, n)
    for (segment in seq_len(3L)) {
      l <- boundaries[segment] + 1L
      r <- boundaries[segment + 1L]
      lasso_training_loss <- lasso_training_loss + colSums(
        relieverChangepoint::reg_fun_lasso_solpath(
          lasso_data, l, r, l, r, lam_set = lam_set
        )$loss
      )
    }
    lasso_order <- order(lam_set)
    lasso_loss <- lasso_loss[lasso_order]
    lasso_training_loss <- lasso_training_loss[lasso_order]
    lasso_best <- which.min(lasso_loss)
    expect_gte(length(lam_set), 20L)
    expect_true(all(diff(lasso_training_loss) >= -1e-7))
    expect_equal(unname(which.min(lasso_training_loss)), 1L)
    expect_gt(lasso_best, 1L)
    expect_lt(lasso_best, length(lasso_loss))
    expect_gt(lasso_loss[1L], lasso_loss[lasso_best])
    expect_gt(lasso_loss[length(lasso_loss)], lasso_loss[lasso_best])

    set.seed(2026)
    p <- 5L
    kde_data <- rbind(
      matrix(stats::rnorm(n_seg * p), n_seg, p),
      matrix(stats::rnorm(n_seg * p, mean = 1), n_seg, p),
      matrix(stats::rnorm(n_seg * p, mean = -1), n_seg, p)
    )
    dist_sq <- as.matrix(stats::dist(kde_data))^2
    bandwidth_vec <- relieverChangepoint:::.kde_resolve_bandwidth_vec(
      NULL, p, dist_sq
    )
    expect_length(bandwidth_vec, 20L)
    expect_true(all(diff(bandwidth_vec) > 0))
    kde_loss <- fixed_segment_cf_loss(
      relieverChangepoint::reg_fun_kde_nll_crossfit,
      dist_sq, c(n_seg, 2L * n_seg), length(bandwidth_vec),
      nfolds = 3L, fold_type = "op",
      bandwidth_vec = bandwidth_vec, var_dim = p
    )
    kde_best <- which.min(kde_loss)
    expect_gt(kde_best, 1L)
    expect_lt(kde_best, length(kde_loss))
    expect_gt(kde_loss[1L], kde_loss[kde_best])
    expect_gt(kde_loss[length(kde_loss)], kde_loss[kde_best])
  }, seconds = 8)
})

test_that("KDE L2 total-loss profiles cannot select an RBF bandwidth", {
  with_test_timeout({
    set.seed(2026)
    n_seg <- 40L
    p <- 4L
    x <- rbind(
      matrix(stats::rnorm(n_seg * p), n_seg, p),
      matrix(stats::rnorm(n_seg * p, mean = 1), n_seg, p),
      matrix(stats::rnorm(n_seg * p, mean = -1), n_seg, p)
    )
    dist_sq <- as.matrix(stats::dist(x))^2
    total_l2_loss <- function(bandwidth) {
      kernel_mat <- exp(-dist_sq / (2 * bandwidth^2))
      boundaries <- c(0L, n_seg, 2L * n_seg, 3L * n_seg)
      sum(vapply(seq_len(3L), function(segment) {
        l <- boundaries[segment] + 1L
        r <- boundaries[segment + 1L]
        sum(relieverChangepoint::reg_fun_kde_l2(
          kernel_mat, l, r, l, r
        )$loss)
      }, numeric(1L)))
    }
    holdout_l2_loss <- function(bandwidth) {
      kernel_mat <- exp(-dist_sq / (2 * bandwidth^2))
      fold_id <- rep(1:3, length.out = nrow(kernel_mat))
      boundaries <- c(0L, n_seg, 2L * n_seg, 3L * n_seg)
      total <- 0
      for (segment in seq_len(3L)) {
        rows <- seq.int(boundaries[segment] + 1L, boundaries[segment + 1L])
        for (fold in 1:3) {
          train_id <- rows[fold_id[rows] != fold]
          eval_id <- rows[fold_id[rows] == fold]
          center <- colMeans(kernel_mat[train_id, , drop = FALSE])
          residual <- sweep(
            kernel_mat[eval_id, , drop = FALSE], 2L, center, "-"
          )
          total <- total + sum(rowMeans(residual^2))
        }
      }
      total
    }

    distance_scale <- stats::median(dist_sq[upper.tri(dist_sq)])
    bandwidth_vec <- exp(seq(
      log(sqrt(distance_scale / 60)),
      log(sqrt(distance_scale * 5e5)),
      length.out = 12L
    ))
    training_loss <- vapply(bandwidth_vec, total_l2_loss, numeric(1L))
    holdout_loss <- vapply(bandwidth_vec, holdout_l2_loss, numeric(1L))

    expect_true(all(is.finite(training_loss)))
    expect_true(all(is.finite(holdout_loss)))
    expect_equal(which.min(training_loss), length(training_loss))
    expect_equal(which.min(holdout_loss), length(holdout_loss))
    expect_lt(training_loss[length(training_loss)], training_loss[1L] * 1e-8)
    expect_lt(holdout_loss[length(holdout_loss)], holdout_loss[1L] * 1e-8)
  }, seconds = 3)
})

test_that("KDE NLL remains stable when kernel values underflow", {
  with_test_timeout({
    dist_sub <- matrix(c(1e6, 1.21e6, 0.81e6, 1e6), nrow = 2L)
    loss <- relieverChangepoint:::.kde_nll_loss(
      dist_sub, bandwidth_vec = c(1 / sqrt(2), 0.5), var_dim = 1L
    )
    expect_true(all(is.finite(loss)))
    expect_gt(min(loss), 8e5)
    expect_true(all(loss[, 2L] > loss[, 1L]))
  })
})

test_that("KDE NLL stable evaluation matches the direct formula", {
  with_test_timeout({
    dist_sub <- matrix(c(
      0.1, 0.4, 0.8,
      0.3, 0.2, 0.6,
      0.7, 0.5, 0.15
    ), nrow = 3L, byrow = TRUE)
    bandwidth_vec <- c(0.25, 0.7, 1.5)
    observed <- relieverChangepoint:::.kde_nll_loss(
      dist_sub, bandwidth_vec, var_dim = 2L
    )
    expected <- vapply(bandwidth_vec, function(bandwidth) {
      two_bandwidth_sq <- 2 * bandwidth^2
      -(-log(pi * two_bandwidth_sq) +
          log(rowMeans(exp(-dist_sub / two_bandwidth_sq))))
    }, numeric(nrow(dist_sub)))
    expect_equal(observed, expected, tolerance = 1e-12)
  }, seconds = 2)
})

test_that("radial Laplace and Student KDE NLL match normalized densities", {
  with_test_timeout({
    dist_sq <- matrix(c(
      0.1, 0.4, 0.8,
      0.3, 0.2, 0.6
    ), nrow = 2L, byrow = TRUE)
    bandwidth_vec <- c(0.4, 1.1)
    p <- 3L

    laplace <- relieverChangepoint:::.kde_nll_loss(
      dist_sq, bandwidth_vec, var_dim = p, kernel = "laplace"
    )
    expected_laplace <- vapply(bandwidth_vec, function(bandwidth) {
      log_normalizer <-
        lgamma(p / 2) - log(2) - (p / 2) * log(pi) -
        lgamma(p) - p * log(bandwidth)
      -(log_normalizer +
          log(rowMeans(exp(-sqrt(dist_sq) / bandwidth))))
    }, numeric(nrow(dist_sq)))
    expect_equal(laplace, expected_laplace, tolerance = 1e-12)

    df <- 7
    student <- relieverChangepoint:::.kde_nll_loss(
      dist_sq, bandwidth_vec, var_dim = p, kernel = "student",
      kernel_args = list(df = df)
    )
    expected_student <- vapply(bandwidth_vec, function(bandwidth) {
      log_normalizer <-
        lgamma((df + p) / 2) - lgamma(df / 2) -
        (p / 2) * log(df * pi) - p * log(bandwidth)
      log_shape <-
        -((df + p) / 2) * log1p(
          dist_sq / (df * bandwidth^2)
        )
      -(log_normalizer + log(rowMeans(exp(log_shape))))
    }, numeric(nrow(dist_sq)))
    expect_equal(student, expected_student, tolerance = 1e-12)
  }, seconds = 2)
})

test_that("KDE NLL remains stable for extremely small bandwidths", {
  with_test_timeout({
    gaussian <- relieverChangepoint:::.kde_nll_loss(
      matrix(c(1e-320, 2e-320), nrow = 1L),
      bandwidth_vec = 1e-300,
      var_dim = 1L,
      kernel = "gaussian"
    )
    expect_true(is.finite(gaussian[1L, 1L]))
    expect_equal(
      gaussian[1L, 1L], 4.999944335913415e279, tolerance = 1e-12
    )

    dist_sq <- matrix(c(1, 4), nrow = 1L)
    bandwidth <- 1e-300
    var_dim <- 3L
    df <- 5
    student <- relieverChangepoint:::.kde_nll_loss(
      dist_sq,
      bandwidth_vec = bandwidth,
      var_dim = var_dim,
      kernel = "student",
      kernel_args = list(df = df)
    )
    log_ratio <- log(dist_sq) - log(df) - 2 * log(bandwidth)
    softplus <- pmax(log_ratio, 0) + log1p(exp(-abs(log_ratio)))
    log_shape <- -((df + var_dim) / 2) * softplus
    log_shape_max <- max(log_shape)
    log_kernel_mean <-
      log_shape_max + log(mean(exp(log_shape - log_shape_max)))
    log_normalizer <-
      lgamma((df + var_dim) / 2) - lgamma(df / 2) -
      (var_dim / 2) * (log(df) + log(pi)) -
      var_dim * log(bandwidth)
    expected_student <- -(log_normalizer + log_kernel_mean)

    expect_true(is.finite(student[1L, 1L]))
    expect_equal(student[1L, 1L], expected_student, tolerance = 1e-12)
  }, seconds = 2)
})

test_that("non-Gaussian KDE defaults use variance-calibrated bandwidths", {
  with_test_timeout({
    set.seed(20260719)
    x <- matrix(stats::rnorm(30L), nrow = 15L)
    dist_sq <- as.matrix(stats::dist(x))^2
    gaussian <- relieverChangepoint:::.kde_resolve_bandwidth_vec(
      NULL, var_dim = 2L, dist_sq = dist_sq
    )
    laplace <- relieverChangepoint:::.kde_resolve_bandwidth_vec(
      NULL, var_dim = 2L, dist_sq = dist_sq, kernel = "laplace"
    )
    student <- relieverChangepoint:::.kde_resolve_bandwidth_vec(
      NULL, var_dim = 2L, dist_sq = dist_sq, kernel = "student",
      kernel_args = list(df = 5)
    )

    expect_equal(laplace, gaussian / sqrt(3), tolerance = 1e-12)
    expect_equal(student, gaussian * sqrt(3 / 5), tolerance = 1e-12)
    expect_error(
      relieverChangepoint:::.kde_resolve_bandwidth_vec(
        NULL, var_dim = 2L, dist_sq = dist_sq, kernel = "student",
        kernel_args = list(df = 2)
      ),
      "bandwidth_vec must be supplied"
    )
  }, seconds = 2)
})

test_that("fixed L2 kernels construct the expected Gram matrices", {
  with_test_timeout({
    x <- matrix(c(0, 1, 2, -1, 1, 3), nrow = 3L, byrow = TRUE)
    dist_sq <- relieverChangepoint:::.kernel_dist_sq(x)
    inner_product <- tcrossprod(x)

    expect_equal(
      relieverChangepoint:::.kernel_l2_matrix(
        x, kernel = "gaussian", bandwidth = 2
      ),
      exp(-dist_sq / 8),
      tolerance = 1e-12
    )
    expect_equal(
      relieverChangepoint:::.kernel_l2_matrix(
        x, kernel = "laplace", bandwidth = 2
      ),
      exp(-sqrt(dist_sq) / 2),
      tolerance = 1e-12
    )
    expect_equal(
      relieverChangepoint:::.kernel_l2_matrix(x, kernel = "linear"),
      inner_product,
      tolerance = 1e-12
    )
    expect_equal(
      relieverChangepoint:::.kernel_l2_matrix(
        x, kernel = "polynomial",
        kernel_args = list(degree = 3, scale = 0.5, offset = 2)
      ),
      (0.5 * inner_product + 2)^3,
      tolerance = 1e-12
    )
    expect_equal(
      relieverChangepoint:::.kernel_l2_matrix(
        x, kernel = "rational_quadratic", bandwidth = 2,
        kernel_args = list(alpha = 1.5)
      ),
      (1 + dist_sq / 12)^(-1.5),
      tolerance = 1e-12
    )
    scaled_matern32 <- sqrt(3 * dist_sq) / 2
    expect_equal(
      relieverChangepoint:::.kernel_l2_matrix(
        x, kernel = "matern32", bandwidth = 2
      ),
      (1 + scaled_matern32) * exp(-scaled_matern32),
      tolerance = 1e-12
    )
    scaled_matern52 <- sqrt(5 * dist_sq) / 2
    expect_equal(
      relieverChangepoint:::.kernel_l2_matrix(
        x, kernel = "matern52", bandwidth = 2
      ),
      (
        1 + scaled_matern52 + 5 * dist_sq / 12
      ) * exp(-scaled_matern52),
      tolerance = 1e-12
    )

    for (kernel in c("matern32", "matern52")) {
      tiny_bandwidth <- relieverChangepoint:::.kernel_l2_matrix(
        x, kernel = kernel, bandwidth = 1e-300
      )
      expect_true(all(is.finite(tiny_bandwidth)))
      expect_equal(diag(tiny_bandwidth), rep(1, nrow(x)))
      expect_equal(
        tiny_bandwidth[row(tiny_bandwidth) != col(tiny_bandwidth)],
        rep(0, nrow(x) * (nrow(x) - 1L))
      )
    }
  }, seconds = 2)
})

test_that("pairwise squared distances are stable under a large translation", {
  x <- rbind(
    c(0, 0),
    c(1, 2),
    c(3, 4)
  )
  y <- rbind(
    c(-1, 1),
    c(2, 5)
  )
  offset <- c(1e12, -1e12)
  shifted_x <- sweep(x, 2L, offset, "+")
  shifted_y <- sweep(y, 2L, offset, "+")

  expected_within <- unname(as.matrix(stats::dist(shifted_x))^2)
  expected_cross <- outer(
    seq_len(nrow(x)),
    seq_len(nrow(y)),
    Vectorize(function(i, j) sum((x[i, ] - y[j, ])^2))
  )

  expect_equal(
    relieverChangepoint:::.kernel_dist_sq(shifted_x),
    expected_within,
    tolerance = 1e-12
  )
  expect_equal(
    relieverChangepoint:::.kernel_dist_sq(shifted_x, shifted_y),
    expected_cross,
    tolerance = 1e-12
  )
})

test_that("pairwise squared distances remain stable with an extreme outlier", {
  x <- rbind(
    c(1e12, -1e12),
    c(0, 0),
    c(1, 2),
    c(3, 4)
  )
  y <- rbind(c(2, 1), c(-1, 3))

  expect_equal(
    relieverChangepoint:::.kernel_dist_sq(x),
    unname(as.matrix(stats::dist(x))^2),
    tolerance = 1e-12
  )
  expect_equal(
    relieverChangepoint:::.kernel_dist_sq(x, y),
    vapply(seq_len(nrow(y)), function(j) {
      rowSums((x - rep(y[j, ], each = nrow(x)))^2)
    }, numeric(nrow(x))),
    tolerance = 1e-12
  )
})

test_that("KDE solpath reg_fun reproduces non-ReCV kernel losses", {
  with_test_timeout({
    set.seed(20260711)
    x <- matrix(stats::rnorm(18L * 2L), nrow = 18L)
    dist_sq <- as.matrix(stats::dist(x))^2
    bandwidth_vec <- c(0.8, 1.3)
    l <- 5L
    r <- 14L
    l_end <- 4L
    r_end <- 16L

    virtual <- relieverChangepoint::reg_fun_kde_nll_solpath(
      dist_sq, 1L, 1L, is_virtual_run = TRUE,
      bandwidth_vec = bandwidth_vec, var_dim = 2L
    )
    expect_equal(virtual$n_loss_outputs, length(bandwidth_vec))
    expect_equal(
      virtual$loss_output_meta$row_type,
      rep("kde_nll", length(bandwidth_vec))
    )
    expect_equal(
      unlist(virtual$loss_output_meta$hyper_value), bandwidth_vec
    )
    expect_equal(
      virtual$loss_output_meta$hyper_name,
      rep("bandwidth", length(bandwidth_vec))
    )
    expect_equal(
      virtual$loss_output_meta$loss_kind,
      rep("negative_log_likelihood", length(bandwidth_vec))
    )
    expect_false("pelt_pruning" %in% names(virtual$loss_output_meta))

    nll_path <- relieverChangepoint::reg_fun_kde_nll_solpath(
      dist_sq, l, r, l_end, r_end,
      bandwidth_vec = bandwidth_vec, var_dim = 2L
    )
    expect_equal(dim(nll_path$loss), c(r_end - l_end + 1L, 2L))
    expect_true(all(is.finite(nll_path$loss)))
  })
})

test_that("direct KDE NLL validation is limited to the used distance block", {
  x <- matrix(seq_len(36L), nrow = 18L, ncol = 2L)
  dist_sq <- as.matrix(stats::dist(x))^2
  dist_sq[18L, 18L] <- -1

  direct <- relieverChangepoint::reg_fun_kde_nll_solpath(
    data = dist_sq, l = 1L, r = 8L, l_end = 1L, r_end = 8L,
    bandwidth_vec = 1, var_dim = 2L
  )
  expect_true(all(is.finite(direct$loss)))

  make_loss_fun <- getFromNamespace(
    ".reliever_make_individual_loss_fun", "relieverChangepoint"
  )
  expect_error(
    make_loss_fun(
      data = dist_sq,
      reg_fun = relieverChangepoint::reg_fun_kde_nll_solpath,
      para_list = list(bandwidth_vec = 1, var_dim = 2L),
      n_loss_outputs = 1L
    ),
    "finite non-negative"
  )
})

test_that("reg_fun_kde_l2 matches native mean loss on a precomputed kernel matrix", {
  with_test_timeout({
    set.seed(20260711)
    x <- matrix(stats::rnorm(20L * 3L), nrow = 20L)
    dist_sq <- as.matrix(stats::dist(x))^2
    kernel_mat <- exp(-0.7 * dist_sq)

    out <- relieverChangepoint::reg_fun_kde_l2(
      kernel_mat, l = 4L, r = 13L, l_end = 3L, r_end = 15L,
      save_model = TRUE
    )
    virtual <- relieverChangepoint::reg_fun_kde_l2(
      kernel_mat, 1L, 1L, is_virtual_run = TRUE
    )
    center <- colMeans(kernel_mat[4:13, , drop = FALSE])
    expected <- rowMeans(
      sweep(kernel_mat[3:15, , drop = FALSE], 2L, center, "-")^2
    )

    expect_equal(dim(out$loss), c(13L, 1L))
    expect_equal(unname(out$loss[, 1L]), unname(expected), tolerance = 1e-12)
    expect_equal(out$model$center, center, tolerance = 1e-12)
    expect_equal(virtual$n_loss_outputs, 1L)
    expect_identical(virtual$loss_output_meta$row_type, "kde_l2")
    expect_identical(virtual$loss_output_meta$loss_kind, "rss")
    expect_false("pelt_pruning" %in% names(virtual$loss_output_meta))

    by_reg_fun <- relieverChangepoint::reliever_generic(
      kernel_mat, cpn_max = 1L, dm = 4L, cov_rate = 0.8,
      reg_fun = relieverChangepoint::reg_fun_kde_l2, method = "SN", echo = FALSE
    )
    by_native <- relieverChangepoint::reliever_mean(
      kernel_mat, cpn_max = 1L, dm = 4L, cov_rate = 0.8,
      method = "SN", echo = FALSE
    )

    expect_same_cpd_path(by_reg_fun, by_native, tolerance = 1e-10)
  })
})

test_that("reg_fun_crossfit_template supports simple custom loss functions", {
  with_test_timeout({
    data <- matrix(seq_len(18L), ncol = 2L)
    loss_fun <- function(data, train_id, eval_id, hyper_set, ...) {
      center <- colMeans(data[train_id, , drop = FALSE])
      base_loss <- rowMeans(
        sweep(data[eval_id, , drop = FALSE], 2L, center, "-")^2
      )
      outer(base_loss, hyper_set$scale)
    }

    out <- relieverChangepoint::reg_fun_crossfit_template(
      data,
      l = 3L,
      r = 8L,
      l_end = 2L,
      r_end = 9L,
      loss_fun = loss_fun,
      hyper_set = data.frame(scale = c(1, 2)),
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      nfolds = 2L,
      fold_type = "blk"
    )
    virtual <- relieverChangepoint::reg_fun_crossfit_template(
      data,
      l = 1L,
      r = 1L,
      loss_fun = loss_fun,
      hyper_set = data.frame(scale = c(1, 2)),
      hyper_name = "scale",
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      is_virtual_run = TRUE
    )

    expect_equal(dim(out$loss), c(8L, 4L))
    expect_true(all(is.finite(out$loss)))
    expect_equal(
      virtual$loss_output_meta$hyper_name,
      c(NA, NA, rep("scale", 2L))
    )
    expect_equal(
      virtual$loss_output_meta$loss_kind,
      c(
        "crossfit_loss", NA,
        rep("crossfit_loss", 2L)
      )
    )
    expect_false("pelt_pruning" %in% names(virtual$loss_output_meta))
    expect_error(
      relieverChangepoint::reg_fun_crossfit_template(
        data,
        l = 3L,
        r = 4L,
        loss_fun = loss_fun,
        hyper_set = data.frame(scale = 1),
        nfolds = 5L,
        fold_type = "blk"
      ),
      "length 2 is smaller than nfolds = 5"
    )
  })
})

test_that("reg_fun_crossfit_template works with reliever dc_grid", {
  with_test_timeout({
    data <- matrix(seq_len(24L), ncol = 1L)
    dc_grid <- seq(2L, 24L, by = 2L)
    loss_fun <- function(data, train_id, eval_id, hyper_set, ...) {
      center <- mean(data[train_id, 1L])
      base_loss <- (data[eval_id, 1L] - center)^2
      outer(base_loss, hyper_set)
    }
    reg_fun <- function(data, l, r, l_end = l, r_end = r,
                        save_model = FALSE, is_virtual_run = FALSE) {
      relieverChangepoint::reg_fun_crossfit_template(
        data = data,
        l = l,
        r = r,
        l_end = l_end,
        r_end = r_end,
        loss_fun = loss_fun,
        hyper_set = c(1, 2),
        nfolds = 2L,
        fold_type = "blk",
        save_model = save_model,
        is_virtual_run = is_virtual_run
      )
    }

    res <- relieverChangepoint::reliever_generic(
      data = data,
      reg_fun = reg_fun,
      cpn_max = 1L,
      dm = 2L,
      cov_rate = 0.6,
      method = "SN",
      echo = FALSE,
      dc_grid = dc_grid
    )

    expect_equal(nrow(res$run_meta), 1L)
    expect_true(all(unlist(res$cpd_path$candidates$cpd) %in% dc_grid))
    expect_true(all(is.finite(res$cpd_path$candidates$loss)))
  })
})

test_that("reg_fun_clf_crossfit_template supports custom probability functions", {
  with_test_timeout({
    data <- matrix(seq_len(20L), ncol = 2L)
    prob_fun <- function(data, y, train_id, eval_id, hyper,
                         probability_scale = 1, ...) {
      rep(
        probability_scale * (mean(y[train_id]) + hyper$offset),
        length(eval_id)
      )
    }

    out <- relieverChangepoint::reg_fun_clf_crossfit_template(
      data,
      l = 4L,
      r = 7L,
      l_end = 3L,
      r_end = 8L,
      prob_fun = prob_fun,
      hyper_set = data.frame(offset = c(0, 0.05)),
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      nfolds = 2L,
      fold_type = "op"
    )
    for (unsupported_fold in c("blk", "stable_blk")) {
      expect_error(
        relieverChangepoint::reg_fun_clf_crossfit_template(
          data,
          l = 1L,
          r = 1L,
          prob_fun = prob_fun,
          hyper_set = data.frame(offset = 0),
          nfolds = 2L,
          fold_type = unsupported_fold,
          is_virtual_run = TRUE
        ),
        "only supports fold_type = \"op\"",
        info = unsupported_fold
      )
    }
    expect_error(
      relieverChangepoint::reg_fun_clf_crossfit_template(
        data,
        l = 4L,
        r = 7L,
        l_end = 3L,
        r_end = 8L,
        prob_fun = prob_fun,
        hyper_set = data.frame(offset = c(0, 0.05)),
        nfolds = 2L,
        fold_type = "op",
        loss = "binary_nll"
      ),
      "density_ratio_nll"
    )
    expect_error(
      relieverChangepoint::reg_fun_clf_crossfit_template(
        data,
        l = 4L,
        r = 7L,
        l_end = 3L,
        r_end = 8L,
        prob_fun = prob_fun,
        hyper_set = data.frame(offset = 0),
        nfolds = 2L,
        fold_type = "op",
        eps = 1e-20
      ),
      "Machine.*double.eps"
    )
    expect_error(
      relieverChangepoint::reg_fun_clf_crossfit_template(
        data,
        l = 4L,
        r = 7L,
        l_end = 3L,
        r_end = 8L,
        prob_fun = prob_fun,
        hyper_set = data.frame(offset = 0),
        nfolds = 2L,
        fold_type = "op",
        op_size = 10L
      ),
      "Reduce op_size"
    )
    boundary_prob <- function(data, y, train_id, eval_id, hyper, ...) {
      rep(1, length(eval_id))
    }
    boundary <- relieverChangepoint::reg_fun_clf_crossfit_template(
      data,
      l = 4L,
      r = 7L,
      l_end = 3L,
      r_end = 8L,
      prob_fun = boundary_prob,
      hyper_set = data.frame(offset = 0),
      nfolds = 2L,
      fold_type = "op",
      eps = .Machine$double.eps
    )
    expect_true(all(is.finite(boundary$loss)))
    virtual <- relieverChangepoint::reg_fun_clf_crossfit_template(
      data,
      l = 1L,
      r = 1L,
      prob_fun = prob_fun,
      hyper_set = data.frame(offset = c(0, 0.05)),
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      is_virtual_run = TRUE
    )

    expect_equal(dim(out$loss), c(6L, 4L))
    expect_true(all(is.finite(out$loss)))
    expect_equal(
      virtual$loss_output_meta$loss_kind,
      c(
        "crossfit_loss", "negative_log_likelihood",
        rep("crossfit_loss", 2L)
      )
    )
    expect_false("pelt_pruning" %in% names(virtual$loss_output_meta))
  })
})

test_that("classifier crossfit excludes complete full-data folds", {
  with_test_timeout({
    data <- matrix(seq_len(12L), ncol = 1L)
    train_ids <- list()
    prob_fun <- function(data, y, train_id, eval_id, hyper, ...) {
      train_ids[[length(train_ids) + 1L]] <<- train_id
      rep(mean(y[train_id]), length(eval_id))
    }

    out <- relieverChangepoint::reg_fun_clf_crossfit_template(
      data,
      l = 5L,
      r = 8L,
      l_end = 5L,
      r_end = 8L,
      prob_fun = prob_fun,
      hyper_set = 1,
      nfolds = 2L,
      fold_type = "op"
    )

    expect_equal(train_ids, list(seq(2L, 12L, by = 2L),
                                 seq(1L, 12L, by = 2L)))
    expect_equal(unname(out$loss), matrix(0, 4L, 1L), tolerance = 1e-12)
  })
})

test_that("classifier density ratios correct each fold's class prior", {
  with_test_timeout({
    set.seed(20260724)
    data <- matrix(stats::rnorm(100L), ncol = 1L)
    null_prob <- function(data, y, train_id, eval_id, hyper, ...) {
      rep(mean(y[train_id]), length(eval_id))
    }
    out <- relieverChangepoint::reg_fun_clf_crossfit_template(
      data,
      l = 21L,
      r = 40L,
      l_end = 21L,
      r_end = 40L,
      prob_fun = null_prob,
      hyper_set = 1,
      nfolds = 2L,
      fold_type = "op"
    )
    whole_interval <- relieverChangepoint::reg_fun_clf_crossfit_template(
      data,
      l = 1L,
      r = nrow(data),
      l_end = 1L,
      r_end = nrow(data),
      prob_fun = null_prob,
      hyper_set = 1,
      nfolds = 2L,
      fold_type = "op"
    )

    expect_equal(unname(out$loss), matrix(0, 20L, 1L), tolerance = 1e-12)
    expect_equal(
      unname(whole_interval$loss),
      matrix(0, nrow(data), 1L),
      tolerance = 1e-12
    )
  })
})

test_that("built-in classifier hyperparameter settings are explicit", {
  with_test_timeout({
    data <- matrix(stats::rnorm(40L), nrow = 10L)

    expect_error(
      relieverChangepoint::reg_fun_ranger_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE,
        hyper_set = c(5, 10)
      ),
      "data frame or matrix with named columns"
    )
    expect_error(
      relieverChangepoint::reg_fun_mlp_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE,
        hyper_set = c(2, 4)
      ),
      "data frame or matrix with named columns"
    )
    expect_error(
      relieverChangepoint::reg_fun_ranger_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE,
        hyper_set = matrix(c(5, 10), ncol = 1L)
      ),
      "uniquely named columns"
    )
    expect_error(
      relieverChangepoint::reg_fun_mlp_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE,
        hyper_set = list(list(2))
      ),
      "named argument lists"
    )
    expect_error(
      relieverChangepoint::reg_fun_ranger_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE,
        hyper_set = data.frame(num.trees = 5L),
        ranger_args = list(num.trees = 10L)
      ),
      "supplied in only one"
    )
    expect_error(
      relieverChangepoint::reg_fun_mlp_crossfit(
        data, 1L, 1L, is_virtual_run = TRUE,
        nnet_args = list(20L)
      ),
      "nnet_args must be a named list"
    )

    virtual <- relieverChangepoint::reg_fun_ranger_crossfit(
      data, 1L, 1L, is_virtual_run = TRUE,
      hyper_set = list(list(), list(num.trees = 5L)),
      loss_output_types = c("recv", "crossfit_homo_hyper")
    )
    expect_equal(virtual$n_loss_outputs, 3L)

    named_matrix <- matrix(
      c(5L, 10L), ncol = 1L,
      dimnames = list(NULL, "num.trees")
    )
    matrix_virtual <- relieverChangepoint::reg_fun_ranger_crossfit(
      data, 1L, 1L, is_virtual_run = TRUE,
      hyper_set = named_matrix,
      loss_output_types = c("recv", "crossfit_homo_hyper")
    )
    expect_equal(matrix_virtual$n_loss_outputs, 3L)
  })
})

test_that("reg_fun_ranger_crossfit exposes classifier ReCV rows", {
  with_test_timeout({
    testthat::skip_if_not_installed("ranger")

    set.seed(20260613)
    data <- matrix(stats::rnorm(28L * 4L), nrow = 28L)
    data[10:18, 1:2] <- data[10:18, 1:2] + 1.5

    ranger_out <- relieverChangepoint::reg_fun_ranger_crossfit(
      data, l = 10L, r = 18L, l_end = 8L, r_end = 21L,
      nfolds = 2L,
      hyper_set = data.frame(num.trees = c(5L, 8L), min.node.size = 2L),
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      ranger_args = list(seed = 20260613L, num.threads = 1L),
      fold_type = "op"
    )
    expect_equal(dim(ranger_out$loss), c(14L, 4L))
    expect_true(all(is.finite(ranger_out$loss)))
    expect_error(
      relieverChangepoint::reg_fun_ranger_crossfit(
        data, l = 10L, r = 18L, l_end = 8L, r_end = 21L,
        nfolds = 2L,
        hyper_set = data.frame(num.trees = 5L, min.node.size = 2L),
        ranger_args = list(class.weights = c(1, 2)),
        fold_type = "op"
      ),
      "empirical training class prior"
    )
    expect_error(
      relieverChangepoint::reg_fun_ranger_crossfit(
        data, l = 10L, r = 18L, l_end = 8L, r_end = 21L,
        nfolds = 2L,
        hyper_set = list(list(data = data)),
        fold_type = "op"
      ),
      "managed by the crossfit wrapper"
    )
  })
})

test_that("reg_fun_mlp_crossfit exposes classifier ReCV rows", {
  with_test_timeout({
    testthat::skip_if_not_installed("nnet")

    set.seed(20260613)
    data <- matrix(stats::rnorm(28L * 4L), nrow = 28L)
    data[10:18, 1:2] <- data[10:18, 1:2] + 1.5

    mlp_out <- relieverChangepoint::reg_fun_mlp_crossfit(
      data, l = 10L, r = 18L, l_end = 8L, r_end = 21L,
      nfolds = 2L,
      hyper_set = data.frame(size = c(2L, 3L), maxit = 20L),
      loss_output_types = c(
        "recv", "incv", "crossfit_homo_hyper"
      ),
      fold_type = "op"
    )
    expect_equal(dim(mlp_out$loss), c(14L, 4L))
    expect_true(all(is.finite(mlp_out$loss)))
    expect_error(
      relieverChangepoint::reg_fun_mlp_crossfit(
        data, l = 10L, r = 18L, l_end = 8L, r_end = 21L,
        nfolds = 2L,
        hyper_set = data.frame(size = 2L, maxit = 20L),
        nnet_args = list(weights = rep(1, nrow(data))),
        fold_type = "op"
      ),
      "empirical training class prior"
    )
    expect_error(
      relieverChangepoint::reg_fun_mlp_crossfit(
        data, l = 10L, r = 18L, l_end = 8L, r_end = 21L,
        nfolds = 2L,
        hyper_set = list(list(x = data)),
        fold_type = "op"
      ),
      "managed by the crossfit wrapper"
    )
  })
})

test_that("the public reg_fun extension API is exported", {
  with_test_timeout({
    exports <- getNamespaceExports("relieverChangepoint")
    expect_true("reg_fun_mean" %in% exports)
    expect_true("reg_fun_lasso_solpath" %in% exports)
    expect_true("reg_fun_lasso_crossfit" %in% exports)
    expect_true("reg_fun_kde_nll_solpath" %in% exports)
    expect_true("reg_fun_kde_l2" %in% exports)
    expect_true("reg_fun_nmcd" %in% exports)
    expect_true("reg_fun_crossfit_template" %in% exports)
    expect_true("reg_fun_clf_crossfit_template" %in% exports)
  })
})

test_that("reg_fun_nmcd wraps the C++ univariate nonparametric loss", {
  with_test_timeout({
    x <- c(1, 3, 2, 5, 4, 7, 6)

    virtual <- relieverChangepoint::reg_fun_nmcd(x, 1L, 1L, is_virtual_run = TRUE)
    out <- relieverChangepoint::reg_fun_nmcd(x, l = 2L, r = 5L, l_end = 1L, r_end = 7L)

    expect_equal(virtual$n_loss_outputs, 1L)
    expect_identical(virtual$loss_output_meta$row_type, "nmcd")
    expect_identical(
      virtual$loss_output_meta$loss_kind,
      "negative_log_likelihood"
    )
    expect_false("pelt_pruning" %in% names(virtual$loss_output_meta))
    expect_equal(dim(out$loss), c(7L, 1L))
    expect_true(all(is.finite(out$loss)))
    expect_error(
      relieverChangepoint::reg_fun_nmcd(
        x, l = 2L, r = 5L, l_end = 1L, r_end = 7L,
        w_trunc = "bad"
      ),
      "w_trunc"
    )
    expect_error(
      relieverChangepoint::reg_fun_nmcd(
        x, l = 2L, r = 5L, l_end = 1L, r_end = 7L,
        sort_X = rev(sort(x))
      ),
      "in nondecreasing order"
    )

    stacked <- c(x, x + 0.25)
    external <- relieverChangepoint::reg_fun_nmcd(
      stacked, l = 1L, r = length(x),
      l_end = length(x) + 1L, r_end = length(stacked),
      sort_X = sort(x)
    )
    expect_equal(dim(external$loss), c(length(x), 1L))
    expect_true(all(is.finite(external$loss)))

    reference_loss <- function(x, l, r, l_end, r_end,
                               sorted_reference, tail_truncation) {
      reference_n <- length(sorted_reference)
      loss <- numeric(r_end - l_end + 1L)
      for (cutpoint in seq.int(
        tail_truncation + 1L,
        reference_n - 1L - tail_truncation
      )) {
        probability <- mean(x[l:r] <= sorted_reference[cutpoint])
        if (abs(probability) <= 1e-15) {
          probability <- 1 / (2 * (r - l + 1L))
        }
        if (abs(probability - 1) <= 1e-15) {
          probability <- 1 - 1 / (2 * (r - l + 1L))
        }
        indicator <- x[l_end:r_end] <= sorted_reference[cutpoint]
        weight <- reference_n / (cutpoint * (reference_n - cutpoint))
        loss <- loss - (
          indicator * log(probability) +
            (1 - indicator) * log1p(-probability)
        ) * weight
      }
      loss
    }

    x_tied <- c(0, 0, 1, 2, 2, 4, 5, 5)
    sorted_reference <- c(-1, 0, 0, 1, 2, 2, 4, 5, 5, 8)
    observed <- relieverChangepoint::reg_fun_nmcd(
      x_tied, l = 2L, r = 7L, l_end = 1L, r_end = 8L,
      w_trunc = 0.1, sort_X = sorted_reference
    )$loss[, 1L]
    expected <- reference_loss(
      x_tied, 2L, 7L, 1L, 8L, sorted_reference, tail_truncation = 1L
    )
    expect_equal(observed, expected, tolerance = 1e-12)
    expect_error(
      relieverChangepoint::reg_fun_nmcd(
        x_tied, l = 0L, r = 7L, l_end = 1L, r_end = 8L,
        sort_X = sorted_reference
      ),
      "out of range"
    )
  })
})

test_that("reg_fun_nmcd excludes appended evaluation rows from its reference", {
  with_test_timeout({
    train <- c(-2, -1, 0, 1, 2, 3)
    eval <- c(50, 60, 70)
    stacked <- c(train, eval)
    attr(stacked, ".reliever_external_train_n") <- length(train)

    marked <- relieverChangepoint::reg_fun_nmcd(
      stacked,
      l = 1L,
      r = length(train),
      l_end = length(train) + 1L,
      r_end = length(stacked)
    )
    explicit <- relieverChangepoint::reg_fun_nmcd(
      stacked,
      l = 1L,
      r = length(train),
      l_end = length(train) + 1L,
      r_end = length(stacked),
      sort_X = sort(train)
    )
    leaky <- relieverChangepoint::reg_fun_nmcd(
      stacked,
      l = 1L,
      r = length(train),
      l_end = length(train) + 1L,
      r_end = length(stacked),
      sort_X = sort(stacked)
    )

    expect_equal(marked$loss, explicit$loss)
    expect_false(isTRUE(all.equal(marked$loss, leaky$loss)))
  })
})
