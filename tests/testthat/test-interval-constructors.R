# Relief and search interval constructors ------------------------------------

test_that("create_relief_itv matches a small golden output", {
  with_test_timeout({
    int_set <- relieverChangepoint::create_relief_itv(20L, 0.8, 5L)

    expect_equal(
      int_set$int_eps[seq_len(8L), ],
      structure(
        c(
          1L, 3L, 5L, 7L, 9L, 11L, 13L, 15L,
          5L, 7L, 9L, 11L, 13L, 15L, 17L, 19L
        ),
        .Dim = c(8L, 2L),
        .Dimnames = list(NULL, c("l", "r"))
      )
    )
    expect_equal(
      int_set$layer_point,
      c(8L, 16L, 23L, 30L, 36L, 40L, 45L, 48L, 50L, 51L)
    )
    expect_equal(
      int_set$int_len,
      c(4L, 5L, 6L, 7L, 8L, 9L, 11L, 12L, 14L, 16L)
    )
    expect_equal(
      int_set$miss_cover_len,
      c(4L, 5L, 6L, 7L, 8L, 10L, 11L, 13L, 15L, 17L)
    )
  })
})

test_that("create_relief_itv full-search warning gives a partial-search action", {
  expect_warning(
    relieverChangepoint::create_relief_itv(20L, 1, 5L),
    "use cov_rate < .* when a partial Reliever search is intended"
  )
})

test_that("WBS and SeedBS interval constructors are deterministic", {
  with_test_timeout({
    expect_equal(
      relieverChangepoint::create_wbs_itv(20L, 5L, 6L, wbs_seed = 99L),
      structure(
        c(0L, 2L, 1L, 5L, 3L, 3L, 15L, 11L, 9L, 12L, 19L, 8L),
        .Dim = c(6L, 2L),
        .Dimnames = list(NULL, c("l", "r"))
      )
    )
    expect_equal(
      relieverChangepoint::create_seed_itv(20L, 5L),
      structure(
        c(0L, 5L, 10L, 0L, 5L, 0L, 10L, 15L, 20L, 15L, 20L, 20L),
        .Dim = c(6L, 2L),
        .Dimnames = list(NULL, c("l", "r"))
      )
    )

    set.seed(1L)
    seed_one <- relieverChangepoint::create_seed_itv(61L, 10L)
    set.seed(2L)
    seed_two <- relieverChangepoint::create_seed_itv(61L, 10L)
    expect_equal(seed_one, seed_two)

    set.seed(2026L)
    expected <- stats::runif(3L)
    set.seed(2026L)
    relieverChangepoint::create_seed_itv(61L, 10L)
    expect_equal(stats::runif(3L), expected)

    set.seed(2026L)
    expected <- stats::runif(3L)
    set.seed(2026L)
    relieverChangepoint::create_wbs_itv(20L, 5L, 6L, wbs_seed = 99L)
    expect_equal(stats::runif(3L), expected)
  })
})

test_that("create_seed_itv keeps the full interval for coarse grids", {
  with_test_timeout({
    expect_equal(
      relieverChangepoint::create_seed_itv(10L, 6L),
      structure(
        c(0L, 10L), .Dim = c(1L, 2L),
        .Dimnames = list(NULL, c("l", "r"))
      )
    )

    res <- relieverChangepoint::reliever_generic(
      matrix(seq_len(10L), ncol = 1L),
      cpn_max = 1L,
      dm = 3L,
      reg_fun = reg_null,
      method = "SeedBS",
      echo = FALSE
    )
    expect_equal(nrow(res$cpd_path$candidates), 2L)
    expect_true(all(is.finite(res$cpd_path$candidates$loss)))
  })
})
