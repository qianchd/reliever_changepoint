make_native_mean_fixture <- function(seed = 20240527L, n = 28L, p = 4L) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p)
  list(x = x)
}

mean_loss_reg <- function(data, l, r, l_end = l, r_end = r,
                          is_virtual_run = FALSE, ...) {
  if (is_virtual_run) {
    return(1L)
  }
  center <- colMeans(data[l:r, , drop = FALSE])
  centered <- sweep(data[l_end:r_end, , drop = FALSE], 2L, center)
  list(loss = matrix(rowMeans(centered^2), ncol = 1L))
}

expect_cpp_result_matches_reliever <- function(cpp, reliever_res,
                                               tolerance = 1e-10) {
  expect_equal(
    as.numeric(cpp_loss(cpp)),
    as.numeric(cpd_path_loss(reliever_res)),
    tolerance = tolerance
  )
  expect_equal(
    as.numeric(cpp_cps_num(cpp)),
    as.numeric(cpd_path_K(reliever_res))
  )
  expect_equal(
    cpp_cpd_list(cpp),
    unclass(cpd_path_run(reliever_res)$cpd)
  )
}

make_wbs_family_intervals <- function(method, n, d, M = 8L, wbs_seed = 77L) {
  switch(
    method,
    WBS = relieverChangepoint::create_wbs_itv(n, 2L * d, M, wbs_seed),
    SeedBS = relieverChangepoint::create_seed_itv(n, 2L * d),
    BS = relieverChangepoint::create_wbs_itv(n, 2L * d, 0L, wbs_seed)
  )
}

decode_cost_col <- function(col, n) {
  interval_index <- col - 1L
  for (j in seq_len(n)) {
    base <- j * (j - 1L) / 2L
    i <- interval_index - base
    if (i >= 1L && i <= j) {
      return(c(i, j))
    }
  }
  stop("Invalid cost matrix column.", call. = FALSE)
}

exact_mean_segment_cost <- function(data, l, r) {
  center <- colMeans(data[l:r, , drop = FALSE])
  centered <- sweep(data[l:r, , drop = FALSE], 2L, center)
  sum(rowMeans(centered^2))
}

expect_finite_cost_mat_matches_exact_mean <- function(cost_mat, data,
                                                      tolerance = 1e-10) {
  finite_cols <- which(is.finite(cost_mat[1L, ]))
  expect_gt(length(finite_cols), 0L)
  actual <- as.numeric(cost_mat[1L, finite_cols])
  expected <- vapply(finite_cols, function(col) {
    interval <- decode_cost_col(col, nrow(data))
    exact_mean_segment_cost(data, interval[1L], interval[2L])
  }, numeric(1))
  expect_equal(actual, expected, tolerance = tolerance)
}

expect_elementwise_relative_loss_error_lte <- function(actual, expected,
                                                       tolerance = 1e-7) {
  expect_equal(length(actual), length(expected))
  relative_error <- abs(actual - expected) / pmax(1, abs(expected))
  expect_lte(max(relative_error), tolerance)
}

test_that("native C++ SN mean helper matches reliever baseline", {
  with_test_timeout({
    n <- 24L
    L <- 2L
    d <- 5L
    cov_rate <- 0.8
    ratio <- 0.9
    fixture <- make_native_mean_fixture(n = n)
    int_set <- relieverChangepoint::create_relief_itv(n, cov_rate, d)

    mean_cost_mat <- matrix(Inf, 1L, 1L + n * (n + 1L) / 2L)
    mean_cpp <- relieverChangepoint:::cpd_mean_by_cost_mat_cpp(
      "SN", n, L, d, matrix(integer(), 0L, 2L), numeric(), 0,
      fixture$x,
      mean_cost_mat,
      ratio,
      int_set$miss_cover_len,
      int_set$int_len,
      int_set$layer_point,
      int_set$int_eps,
      FALSE
    )
    mean_r <- relieverChangepoint::reliever_generic(
      data = fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = cov_rate,
      reg_fun = mean_loss_reg,
      method = "SN",
      run_cpd_ids = 1L,
      echo = FALSE
    )
    expect_cpp_result_matches_reliever(mean_cpp, mean_r)
  })
})

test_that("native mean cost engines store exact full-search segment losses", {
  with_test_timeout({
    fixture <- make_native_mean_fixture(seed = 20240601L, n = 18L, p = 3L)
    n <- nrow(fixture$x)
    L <- 2L
    d <- 4L
    ratio <- 0.9
    empty_int_eps <- matrix(integer(), 0L, 2L)

    sn_cost_mat <- matrix(Inf, 1L, 1L + n * (n + 1L) / 2L)
    relieverChangepoint:::cpd_mean_by_cost_mat_cpp(
      "SN", n, L, d, matrix(integer(), 0L, 2L), numeric(), 0,
      fixture$x,
      sn_cost_mat,
      ratio,
      integer(),
      integer(),
      integer(),
      empty_int_eps,
      TRUE
    )
    expect_finite_cost_mat_matches_exact_mean(sn_cost_mat, fixture$x)

    seed_intervals <- relieverChangepoint::create_seed_itv(n, 2L * d)
    seedbs_cost_mat <- matrix(Inf, 1L, 1L + n * (n + 1L) / 2L)
    relieverChangepoint:::cpd_mean_by_cost_mat_cpp(
      "SeedBS", n, L, d, seed_intervals, numeric(), 0,
      fixture$x,
      seedbs_cost_mat,
      ratio,
      integer(),
      integer(),
      integer(),
      empty_int_eps,
      TRUE
    )
    expect_finite_cost_mat_matches_exact_mean(seedbs_cost_mat, fixture$x)
  })
})

test_that("native mean loss is stable under large common offsets", {
  with_test_timeout({
    set.seed(1)
    centered <- matrix(stats::rnorm(30L * 2L), nrow = 30L)
    shifted <- centered + 1e8
    common <- list(
      cpn_max = 3L,
      dm = 3L,
      cov_rate = 1,
      method = "SN",
      cache_backend = "by_cost_mat"
    )
    native <- do.call(
      relieverChangepoint::reliever_mean,
      c(list(data = shifted), common)
    )
    generic <- do.call(
      relieverChangepoint::reliever_generic,
      c(
        list(
          data = shifted,
          reg_fun = relieverChangepoint::reg_fun_mean
        ),
        common
      )
    )
    unshifted <- do.call(
      relieverChangepoint::reliever_mean,
      c(list(data = centered), common)
    )

    expect_elementwise_relative_loss_error_lte(
      native$cpd_path$candidates$loss,
      generic$cpd_path$candidates$loss,
      tolerance = 1e-7
    )
    expect_elementwise_relative_loss_error_lte(
      native$cpd_path$candidates$loss,
      unshifted$cpd_path$candidates$loss,
      tolerance = 1e-7
    )
    expect_true(all(native$cpd_path$candidates$loss >= 0))
    expect_silent(
      relieverChangepoint::select_by_run(native, cpn_crit = "rss_sic")
    )

    set.seed(1)
    piecewise_shift <- c(
      stats::rnorm(20L),
      stats::rnorm(20L, 1e8, 1),
      stats::rnorm(20L, -1e8, 1)
    )
    piecewise_common <- list(
      cpn_max = 5L,
      dm = 5L,
      cov_rate = 1,
      method = "SN",
      cache_backend = "by_cost_mat"
    )
    piecewise_native <- do.call(
      relieverChangepoint::reliever_mean,
      c(list(data = piecewise_shift), piecewise_common)
    )
    piecewise_generic <- do.call(
      relieverChangepoint::reliever_generic,
      c(
        list(
          data = piecewise_shift,
          reg_fun = relieverChangepoint::reg_fun_mean
        ),
        piecewise_common
      )
    )
    expect_elementwise_relative_loss_error_lte(
      piecewise_native$cpd_path$candidates$loss,
      piecewise_generic$cpd_path$candidates$loss,
      tolerance = 1e-7
    )
    expect_true(all(piecewise_native$cpd_path$candidates$loss >= 0))

    set.seed(28)
    stress_levels <- c(-1e12, -1e12, -1e12, -1e12, 0, 0)
    stress_data <- unlist(lapply(
      stress_levels,
      function(level) stats::rnorm(6L, level, 1)
    ))
    stress_common <- list(
      cpn_max = 4L,
      dm = 3L,
      cov_rate = 0.8,
      method = "SN",
      cache_backend = "by_loss_block",
      ratio = 0.9
    )
    stress_native <- do.call(
      relieverChangepoint::reliever_mean,
      c(list(data = stress_data), stress_common)
    )
    stress_common$ratio <- NULL
    stress_generic <- do.call(
      relieverChangepoint::reliever_generic,
      c(
        list(
          data = stress_data,
          reg_fun = relieverChangepoint::reg_fun_mean
        ),
        stress_common
      )
    )
    expect_elementwise_relative_loss_error_lte(
      stress_native$cpd_path$candidates$loss,
      stress_generic$cpd_path$candidates$loss,
      tolerance = 1e-7
    )

    set.seed(2)
    mixed_levels <- sample(c(-1e12, 0, 1e12), 6L, replace = TRUE)
    mixed_data <- unlist(lapply(
      mixed_levels,
      function(level) stats::rnorm(6L, level, 1)
    ))
    mixed_common <- list(
      cpn_max = 4L,
      dm = 3L,
      cov_rate = 0.8,
      method = "SN",
      cache_backend = "by_loss_block"
    )
    mixed_native <- do.call(
      relieverChangepoint::reliever_mean,
      c(list(data = mixed_data, ratio = 0.2), mixed_common)
    )
    mixed_generic <- do.call(
      relieverChangepoint::reliever_generic,
      c(
        list(data = mixed_data, reg_fun = relieverChangepoint::reg_fun_mean),
        mixed_common
      )
    )
    expect_elementwise_relative_loss_error_lte(
      mixed_native$cpd_path$candidates$loss,
      mixed_generic$cpd_path$candidates$loss,
      tolerance = 1e-7
    )
  })
})

test_that("native mean backends recover strong piecewise-constant changepoints", {
  with_test_timeout({
    data <- rbind(
      matrix(c(0, 0), 30L, 2L, byrow = TRUE),
      matrix(c(5, -5), 30L, 2L, byrow = TRUE),
      matrix(c(-4, 4), 30L, 2L, byrow = TRUE)
    )
    true_cps <- c(30L, 60L)
    method_specs <- list(
      SN = list(),
      BS = list(M = 0L),
      WBS = list(M = 30L, wbs_seed = 123L),
      SeedBS = list(M = 30L)
    )

    for (method in names(method_specs)) {
      common_args <- c(
        list(
          data = data,
          cpn_max = length(true_cps),
          dm = 10L,
          cov_rate = 0.8,
          method = method,
          ratio = 0.9,
          detail = TRUE,
          echo = FALSE
        ),
        method_specs[[method]]
      )

      by_cost_mat <- do.call(
        relieverChangepoint::reliever_mean,
        c(common_args, list(cache_backend = "by_cost_mat"))
      )
      by_loss_block <- do.call(
        relieverChangepoint::reliever_mean,
        c(common_args, list(cache_backend = "by_loss_block"))
      )

      by_cost_mat_k2 <- cpd_path_run(by_cost_mat)
      by_cost_mat_k2 <- by_cost_mat_k2[by_cost_mat_k2$K == 2L, ,
                                       drop = FALSE]
      by_loss_block_k2 <- cpd_path_run(by_loss_block)
      by_loss_block_k2 <- by_loss_block_k2[
        by_loss_block_k2$K == 2L, , drop = FALSE
      ]

      expect_equal(by_cost_mat_k2$loss, 0, tolerance = 1e-12)
      expect_equal(by_loss_block_k2$loss, 0, tolerance = 1e-12)
      expect_cpd_error_lte(by_cost_mat_k2$cpd[[1L]], true_cps)
      expect_cpd_error_lte(by_loss_block_k2$cpd[[1L]], true_cps)
      expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-12)
    }
  })
})

test_that("native C++ SN mean by_loss_block backend matches by_cost_mat mean baseline", {
  with_test_timeout({
    n <- 30L
    L <- 2L
    d <- 6L
    cov_rate <- 0.8
    ratio <- 0.9
    fixture <- make_native_mean_fixture(seed = 20240528L, n = n)
    int_set <- relieverChangepoint::create_relief_itv(n, cov_rate, d)

    by_cost_mat_cost_mat <- matrix(Inf, 1L, 1L + n * (n + 1L) / 2L)
    by_cost_mat_cpp <- relieverChangepoint:::cpd_mean_by_cost_mat_cpp(
      "SN", n, L, d, matrix(integer(), 0L, 2L), numeric(), 0,
      fixture$x,
      by_cost_mat_cost_mat,
      ratio,
      int_set$miss_cover_len,
      int_set$int_len,
      int_set$layer_point,
      int_set$int_eps,
      FALSE
    )
    by_loss_block_cpp <- relieverChangepoint:::cpd_mean_by_loss_block_cpp(
      "SN", n, L, d, matrix(integer(), 0L, 2L), numeric(), 0,
      fixture$x,
      ratio,
      int_set$miss_cover_len,
      int_set$int_len,
      int_set$layer_point,
      int_set$int_eps,
      return_cache_profile = TRUE,
      use_owner_key = TRUE
    )
    reliever_r <- relieverChangepoint::reliever_generic(
      data = fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = cov_rate,
      reg_fun = mean_loss_reg,
      method = "SN",
      run_cpd_ids = 1L,
      echo = FALSE
    )

    expect_equal(cpp_loss(by_loss_block_cpp), cpp_loss(by_cost_mat_cpp),
                 tolerance = 1e-10)
    expect_equal(cpp_cps_num(by_loss_block_cpp), cpp_cps_num(by_cost_mat_cpp))
    expect_equal(cpp_cpd_list(by_loss_block_cpp),
                 cpp_cpd_list(by_cost_mat_cpp))
    expect_cpp_result_matches_reliever(by_loss_block_cpp, reliever_r)
    expect_true(by_loss_block_cpp$n_model_fit > 0)
  })
})

test_that("reliever_mean exposes by_cost_mat and by_loss_block native mean backends", {
  with_test_timeout({
    n <- 30L
    L <- 2L
    d <- 6L
    cov_rate <- 0.8
    ratio <- 0.9
    fixture <- make_native_mean_fixture(seed = 20240529L, n = n)

    generic <- relieverChangepoint::reliever_generic(
      data = fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = cov_rate,
      reg_fun = mean_loss_reg,
      method = "SN",
      run_cpd_ids = 1L,
      echo = FALSE
    )
    by_cost_mat <- relieverChangepoint::reliever_mean(
      fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = cov_rate,
      method = "SN",
      ratio = ratio,
      cache_backend = "by_cost_mat",
      owner_key = FALSE,
      detail = TRUE,
      echo = FALSE
    )
    by_loss_block <- relieverChangepoint::reliever_mean(
      fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = cov_rate,
      method = "SN",
      ratio = ratio,
      cache_backend = "by_loss_block",
      detail = TRUE,
      echo = FALSE
    )
    no_owner_key <- relieverChangepoint::reliever_mean(
      fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = cov_rate,
      method = "SN",
      ratio = ratio,
      cache_backend = "by_loss_block",
      owner_key = FALSE,
      detail = TRUE,
      echo = FALSE
    )

    expect_same_cpd_path(by_cost_mat, generic, tolerance = 1e-10)
    expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-10)
    expect_same_cpd_path(no_owner_key, by_loss_block, tolerance = 1e-10)
    expect_false("cost_mat" %in% names(by_cost_mat))
    expect_false("int_set" %in% names(by_cost_mat))
    expect_true(is.matrix(by_cost_mat$cache_profile$objects$cost_mat))
    expect_false(is.null(by_loss_block$cache_profile$objects$int_set))
    expect_true(
      "owner_key" %in%
        names(by_loss_block$cache_profile$objects$loss_block_cache)
    )
    expect_false(
      "owner_key" %in%
        names(no_owner_key$cache_profile$objects$loss_block_cache)
    )
    expect_error(
      relieverChangepoint::reliever_mean(
        fixture$x, cpn_max = L, dm = d, cov_rate = cov_rate, owner_key = NA
      ),
      "owner_key"
    )

    switched <- NULL
    expect_warning(
      switched <- relieverChangepoint::reliever_mean(
        fixture$x,
        cpn_max = L,
        dm = d,
        cov_rate = 1,
        method = "SN",
        ratio = ratio,
        cache_backend = "by_loss_block",
        echo = FALSE
      ),
      "using cache_backend = \"by_cost_mat\""
    )
    expect_equal(switched$settings$cache_backend, "by_cost_mat")
    expect_true(all(is.finite(cpd_path_loss(switched))))
  })
})

test_that("native C++ WBS-family mean by_loss_block backend matches by_cost_mat baselines", {
  with_test_timeout({
    n <- 30L
    L <- 2L
    d <- 6L
    cov_rate <- 0.8
    ratio <- 0.9
    fixture <- make_native_mean_fixture(seed = 20240530L, n = n)
    int_set <- relieverChangepoint::create_relief_itv(n, cov_rate, d)

    for (method in c("WBS", "SeedBS", "BS")) {
      lr_m <- make_wbs_family_intervals(method, n, d, M = 8L, wbs_seed = 77L)
      by_cost_mat_cost_mat <- matrix(Inf, 1L, 1L + n * (n + 1L) / 2L)
      by_cost_mat_cpp <- relieverChangepoint:::cpd_mean_by_cost_mat_cpp(
        method, n, L, d, lr_m, numeric(), 0,
        fixture$x,
        by_cost_mat_cost_mat,
        ratio,
        int_set$miss_cover_len,
        int_set$int_len,
        int_set$layer_point,
        int_set$int_eps,
        FALSE
      )
      by_loss_block_cpp <- relieverChangepoint:::cpd_mean_by_loss_block_cpp(
        method, n, L, d, lr_m, numeric(), 0,
        fixture$x,
        ratio,
        int_set$miss_cover_len,
        int_set$int_len,
        int_set$layer_point,
        int_set$int_eps,
        return_cache_profile = TRUE,
        use_owner_key = TRUE
      )

      expect_equal(cpp_loss(by_loss_block_cpp), cpp_loss(by_cost_mat_cpp),
                   tolerance = 1e-10)
      expect_equal(cpp_cps_num(by_loss_block_cpp),
                   cpp_cps_num(by_cost_mat_cpp))
      expect_equal(cpp_cpd_list(by_loss_block_cpp),
                   cpp_cpd_list(by_cost_mat_cpp))
    expect_true(by_loss_block_cpp$n_model_fit > 0)
    }
  })
})

test_that("reliever_mean WBS-family methods match generic by_cost_mat reliever", {
  with_test_timeout({
    n <- 30L
    L <- 2L
    d <- 6L
    cov_rate <- 0.8
    ratio <- 0.9
    fixture <- make_native_mean_fixture(seed = 20240531L, n = n)

    method_specs <- list(
      WBS = list(M = 8L, wbs_seed = 77L),
      SeedBS = list(M = 8L),
      BS = list(M = 0L)
    )

    for (method in names(method_specs)) {
      generic <- do.call(
        relieverChangepoint::reliever_generic,
        c(
          list(
            data = fixture$x,
            cpn_max = L,
            dm = d,
            cov_rate = cov_rate,
            reg_fun = mean_loss_reg,
            method = method,
            run_cpd_ids = 1L,
            echo = FALSE
          ),
          method_specs[[method]]
        )
      )
      by_cost_mat <- do.call(
        relieverChangepoint::reliever_mean,
        c(
          list(
            data = fixture$x,
            cpn_max = L,
            dm = d,
            cov_rate = cov_rate,
            method = method,
            ratio = ratio,
            cache_backend = "by_cost_mat",
            detail = TRUE,
            echo = FALSE
          ),
          method_specs[[method]]
        )
      )
      by_loss_block <- do.call(
        relieverChangepoint::reliever_mean,
        c(
          list(
            data = fixture$x,
            cpn_max = L,
            dm = d,
            cov_rate = cov_rate,
            method = method,
            ratio = ratio,
            cache_backend = "by_loss_block",
            detail = TRUE,
            echo = FALSE
          ),
          method_specs[[method]]
        )
      )

      expect_same_cpd_path(by_cost_mat, generic, tolerance = 1e-10)
      expect_same_cpd_path(by_loss_block, by_cost_mat, tolerance = 1e-10)
    }
  })
})

test_that("native C++ WBS and SeedBS mean dispatch matches full reliever", {
  with_test_timeout({
    n <- 28L
    L <- 2L
    d <- 5L
    ratio <- 0.9
    fixture <- make_native_mean_fixture(n = n)

    wbs_seed <- 77L
    M <- 8L
    wbs_intervals <- relieverChangepoint::create_wbs_itv(n, 2L * d, M, wbs_seed)
    wbs_cost_mat <- matrix(Inf, 1L, 1L + n * (n + 1L) / 2L)
    wbs_cpp <- relieverChangepoint:::cpd_mean_by_cost_mat_cpp(
      "WBS", n, L, d, wbs_intervals, numeric(), 0,
      fixture$x,
      wbs_cost_mat,
      ratio,
      integer(),
      integer(),
      integer(),
      matrix(integer(), 0L, 2L),
      TRUE
    )
    wbs_r <- relieverChangepoint::reliever_generic(
      data = fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = 1,
      reg_fun = mean_loss_reg,
      method = "WBS",
      M = M,
      wbs_seed = wbs_seed,
      run_cpd_ids = 1L,
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    expect_cpp_result_matches_reliever(wbs_cpp, wbs_r)

    seed_intervals <- relieverChangepoint::create_seed_itv(n, 2L * d)
    sbs_cost_mat <- matrix(Inf, 1L, 1L + n * (n + 1L) / 2L)
    sbs_cpp <- relieverChangepoint:::cpd_mean_by_cost_mat_cpp(
      "SeedBS", n, L, d, seed_intervals, numeric(), 0,
      fixture$x,
      sbs_cost_mat,
      ratio,
      integer(),
      integer(),
      integer(),
      matrix(integer(), 0L, 2L),
      TRUE
    )
    seedbs_r <- relieverChangepoint::reliever_generic(
      data = fixture$x,
      cpn_max = L,
      dm = d,
      cov_rate = 1,
      reg_fun = mean_loss_reg,
      method = "SeedBS",
      run_cpd_ids = 1L,
      cache_backend = "by_cost_mat",
      echo = FALSE
    )
    expect_cpp_result_matches_reliever(sbs_cpp, seedbs_r)
  })
})
