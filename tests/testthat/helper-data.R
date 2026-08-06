test_scale_snr <- function(b, snr, sig = 1) {
  b <- as.vector(b)
  sqrt(snr / sum(b^2)) * sig * b
}
