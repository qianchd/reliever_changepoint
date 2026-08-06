# Relief interval ownership algebra ------------------------------------------
#
# These exhaustive tests belong to the exact-to-relief and relief-to-exact
# mapping layer. They are intentionally excluded from routine quick/full runs.
# Run them after changing create_relief_itv(), either mapping routine, or cache
# ownership logic:
#   Sys.setenv(RELIEVER_TEST_MODE = "interval")
#   devtools::test(filter = "interval-routines")

exact2relief_itv_test <- function(l, r, int_set) {
  relieverChangepoint:::exact2relief_itv_routine_c(
    l, r, int_set$miss_cover_len, int_set$int_len,
    int_set$layer_point, int_set$int_eps
  )
}

relief2exact_itv_test <- function(id, int_set) {
  relieverChangepoint:::relief2exact_itv_routine_c(
    id, int_set$int_len, int_set$layer_point, int_set$int_eps, int_set$n
  )
}

expand_relief_cells <- function(id, int_set) {
  saw <- relief2exact_itv_test(id, int_set)
  if (is.null(dim(saw))) {
    saw <- matrix(saw, ncol = 2L)
  }

  l <- int_set$int_eps[id, 1L] + 1L
  r <- int_set$int_eps[id, 2L]
  chunks <- vector("list", nrow(saw))
  for (s in seq_len(nrow(saw))) {
    ll_to <- max(1L, saw[s, 1L])
    rr_from <- if (s == 1L) r else saw[s - 1L, 2L] + 1L
    rr_to <- min(int_set$n, saw[s, 2L])
    if (ll_to <= l && rr_from <= rr_to) {
      chunks[[s]] <- as.matrix(expand.grid(
        ll = seq(l, ll_to, by = -1L),
        rr = rr_from:rr_to
      ))
    }
  }
  chunks <- Filter(Negate(is.null), chunks)
  if (length(chunks) == 0L) {
    return(matrix(integer(), 0L, 2L))
  }
  do.call(rbind, chunks)
}

relief_partition_diagnostics <- function(n, cov_rate, dm) {
  int_set <- suppressWarnings(
    relieverChangepoint::create_relief_itv(n, cov_rate, dm)
  )
  owner <- matrix(NA_integer_, n, n)
  inverse_errors <- matrix(integer(), 0L, 4L)
  duplicate_cells <- matrix(integer(), 0L, 4L)

  for (id in seq_len(nrow(int_set$int_eps))) {
    cells <- expand_relief_cells(id, int_set)
    if (nrow(cells) == 0L) {
      next
    }
    cells <- cells[cells[, 2L] - cells[, 1L] + 1L >= dm, , drop = FALSE]
    for (z in seq_len(nrow(cells))) {
      ll <- cells[z, 1L]
      rr <- cells[z, 2L]
      inverse_id <- as.integer(exact2relief_itv_test(
        ll - 1L, rr, int_set
      )[1L])
      if (inverse_id != id) {
        inverse_errors <- rbind(
          inverse_errors, c(ll, rr, expected = id, actual = inverse_id)
        )
      }
      if (!is.na(owner[ll, rr])) {
        duplicate_cells <- rbind(
          duplicate_cells, c(ll, rr, first = owner[ll, rr], second = id)
        )
      } else {
        owner[ll, rr] <- id
      }
    }
  }

  eligible <- row(owner) <= col(owner) &
    col(owner) - row(owner) + 1L >= dm
  list(
    inverse_errors = inverse_errors,
    duplicate_cells = duplicate_cells,
    missing_cells = which(eligible & is.na(owner), arr.ind = TRUE)
  )
}

expect_relief_partition <- function(n, cov_rate, dm) {
  diagnostics <- relief_partition_diagnostics(n, cov_rate, dm)
  info <- sprintf("n=%d, cov_rate=%g, dm=%d", n, cov_rate, dm)
  expect_equal(nrow(diagnostics$inverse_errors), 0L, info = info)
  expect_equal(nrow(diagnostics$duplicate_cells), 0L, info = info)
  expect_equal(nrow(diagnostics$missing_cells), 0L, info = info)
}

test_that("relief interval routines match small golden mappings", {
  skip_if_not_interval_tests()
  with_test_timeout({
    int_set <- relieverChangepoint::create_relief_itv(20L, 0.8, 5L)
    expect_equal(
      exact2relief_itv_test(0L, 20L, int_set),
      structure(c(51L, 2L, 18L), .Dim = c(3L, 1L))
    )
    expect_equal(
      exact2relief_itv_test(3L, 15L, int_set),
      structure(c(43L, 4L, 15L), .Dim = c(3L, 1L))
    )
    expect_equal(
      relief2exact_itv_test(1L, int_set),
      structure(c(2L, 6L), .Dim = c(1L, 2L))
    )
    expect_equal(
      relief2exact_itv_test(10L, int_set),
      structure(c(3L, 8L), .Dim = c(1L, 2L))
    )

    special_set <- relieverChangepoint::create_relief_itv(120L, 0.9, 15L)
    relief <- exact2relief_itv_test(
      0L, 113L, special_set
    )
    expect_equal(relief, structure(c(880L, 4L, 109L), .Dim = c(3L, 1L)))
    special_cells <- expand_relief_cells(relief[1L], special_set)
    expect_true(any(special_cells[, 1L] == 1L & special_cells[, 2L] == 113L))
  })
})

test_that("cov_rate 1 enumerates each eligible interval exactly once", {
  skip_if_not_interval_tests()
  with_test_timeout({
    for (dm in c(2L, 5L)) {
      n <- 20L
      int_set <- suppressWarnings(
        relieverChangepoint::create_relief_itv(n, 1, dm)
      )
      expected_eps <- do.call(
        rbind,
        lapply(dm:n, function(len) cbind(0L:(n - len), len:n))
      )
      mode(expected_eps) <- "integer"
      expect_equal(unname(int_set$int_eps), unname(expected_eps))
      expect_false(any(duplicated(data.frame(int_set$int_eps))))
      expect_relief_partition(n, 1, dm)
    }
  })
})

test_that("relief interval routines partition representative designs", {
  skip_if_not_interval_tests()
  with_test_timeout({
    specifications <- list(
      c(n = 18L, cov_rate = 0.6, dm = 3L),
      c(n = 30L, cov_rate = 1.0, dm = 1L),
      c(n = 30L, cov_rate = 1.0, dm = 3L),
      c(n = 30L, cov_rate = 1.0, dm = 10L),
      c(n = 20L, cov_rate = 0.5, dm = 1L),
      c(n = 30L, cov_rate = 0.6, dm = 1L),
      c(n = 30L, cov_rate = 0.9, dm = 1L),
      c(n = 30L, cov_rate = 0.6, dm = 5L)
    )
    for (spec in specifications) {
      expect_relief_partition(
        as.integer(spec[["n"]]),
        spec[["cov_rate"]],
        as.integer(spec[["dm"]])
      )
    }
  })
})
