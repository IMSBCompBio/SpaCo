compute_C <- function(x, neighbourindexmatrix) {
  n <- nrow(neighbourindexmatrix)
  W <- sum(neighbourindexmatrix)
  active_Indices <- which(neighbourindexmatrix != 0)
  i <- active_Indices %% n
  i[which(i == 0)] <- n
  j <- ceiling(active_Indices / n)
  C <- (((n - 1) / (2 * n * W)) *
          sum(neighbourindexmatrix[active_Indices] * (x[i] - x[j])^2)) / var(x)
  return(C)
}
