reliever_test_mode <- function() {
  mode <- tolower(Sys.getenv("RELIEVER_TEST_MODE", "quick"))
  valid_modes <- c("quick", "full", "interval", "all")
  if (!mode %in% valid_modes) {
    stop(
      "RELIEVER_TEST_MODE must be 'quick', 'full', 'interval', or 'all'.",
      call. = FALSE
    )
  }
  mode
}

skip_if_not_full_tests <- function(reason = "full test suite only") {
  testthat::skip_if(
    !reliever_test_mode() %in% c("full", "all"),
    paste0(reason, "; set RELIEVER_TEST_MODE=full or all to run")
  )
}

skip_if_not_interval_tests <- function(reason = "relief interval suite only") {
  testthat::skip_if(
    !reliever_test_mode() %in% c("interval", "all"),
    paste0(reason, "; set RELIEVER_TEST_MODE=interval or all to run")
  )
}
