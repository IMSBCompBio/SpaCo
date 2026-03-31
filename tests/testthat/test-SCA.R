test_that("spaco works", {
  n_row <- 5
  n_col <- 5
  n <- n_row * n_col
  p <- 10
  n_total <- n * p
  genes <- matrix(rnorm(n_total), ncol = p, nrow = n)
  coords <- expand.grid(1:n_row, 1:n_col)
  spotnames <- paste0("spot_", 1:n)
  genenames <- paste0("gene_", 1:p)
  rownames(coords) <- spotnames
  colnames(coords) <- c("x", "y")
  rownames(genes) <- spotnames
  colnames(genes) <- genenames
  d_row <- abs(outer(coords[, "x"], coords[, "x"], "-"))
  d_col <- abs(outer(coords[, "y"], coords[, "y"], "-"))
  W <- (d_row <= 1 & d_col <= 1 & (d_row + d_col) > 0) * 1
  object <- SpaCoObject(W, genes, coords)
  object <- SpaCo(object)
  W_Matrix <- Matrix::Matrix(W)
  object_Matrix <- SpaCoObject(W_Matrix, genes, coords)
  object_Matrix <- SpaCo(object_Matrix)
  expect_equal(object@Lambdas, object_Matrix@Lambdas)
  expect_error(SpaCoObject(cbind(W, 1), genes, coords))
  expect_error(SpaCoObject(rbind(cbind(W, 1), 1), genes, coords))
  test <- SVGTest(object)
  test_Matrix <- SVGTest(object_Matrix)
  expect_equal(test$score, test_Matrix$score)
  new_genes <- genes <- matrix(rnorm(n_total), ncol = p, nrow = n)
  colnames(new_genes) <- paste0("gene_", 1:p)
  object_more <- set_projection(object, new_genes)
  object_more_Matrix <- set_projection(object_Matrix, new_genes)
  test_more <- SVGTest(object_more)
  test_more_Matrix <- SVGTest(object_more_Matrix)
  expect_equal(test_more$score, test_more_Matrix$score)
  Spaco_plot(object, 1)
  Spaco_plot(object_Matrix, 1)
  object <- denoise_profiles(object)
  object_Matrix <- denoise_profiles(object_Matrix)
  feature_plot(object, "gene_1")
  feature_plot(object_Matrix, "gene_1")
  denoised_projection_plot(object, "gene_1")
  denoised_projection_plot(object_Matrix, "gene_1")
})
test_that("rcpp function", {
  A <- matrix(1:6, 2, 3)
  B <- matrix(1:4, 2, 2)
  expect_error(eigenMapMatMult(A, B))
})
