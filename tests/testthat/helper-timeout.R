with_test_timeout <- function(code, seconds = 120) {
  setTimeLimit(cpu = seconds, elapsed = seconds, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  force(code)
}
