reg_null <- function(data, l, r, l_end = l, r_end = r,
                     is_virtual_run = FALSE, ...) {
  if (is_virtual_run) {
    return(1L)
  }
  list(loss = matrix(0, r_end - l_end + 1L, 1L))
}
