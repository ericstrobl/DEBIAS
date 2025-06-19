cor_to_pval <- function(r, n) {
  if (abs(r) == 1){
    return(0)
  }
  t_stat <- r * sqrt((n - 2) / (1 - r^2))
  pval <- 2 * (1 - pt(abs(t_stat), df = n - 2))
  return(pval)
}
