normalizeA <- function(x, A)
{
  return((x / rep(normA(x, A), length(x))))
}
