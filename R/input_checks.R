#' @keywords internal
# ---- Helper function for input checks ----
check_RunSCA_args <- function(args) {
  # SpaCoObject
  if (!inherits(args$SpaCoObject, "SpaCoObject")) {
    stop("`SpaCoObject` must be an object of class 'SpaCoObject'.")
  }

  # PC_criterion
  if (!args$PC_criterion %in% c("percent", "number")) {
    stop("`PC_criterion` must be either 'percent' or 'number'.")
  }

  # PC_value
  if (args$PC_criterion == "percent") {
    if (!is.numeric(args$PC_value) || length(args$PC_value) != 1 ||
        args$PC_value <= 0 || args$PC_value >= 1) {
      stop("When `PC_criterion` is 'percent', `PC_value` must be a numeric scalar > 0 and < 1.")
    }
  } else if (args$PC_criterion == "number") {
    if (!is.numeric(args$PC_value) ||
        length(args$PC_value) != 1 ||
        args$PC_value != as.integer(args$PC_value)) {
      stop("When `PC_criterion` is 'number', `PC_value` must be an integer scalar.")
    }
  }

  # set_nspacs
  if (!is.null(args$set_nspacs)) {
    if (!is.numeric(args$set_nspacs) ||
        length(args$set_nspacs) != 1 ||
        args$set_nspacs != as.integer(args$set_nspacs)) {
      stop("`set_nspacs` must be either NULL or a single integer value.")
    }
  }

  # nSim
  if (!is.numeric(args$nSim) ||
      length(args$nSim) != 1 ||
      args$nSim != as.integer(args$nSim)) {
    stop("`nSim` must be a single integer value.")
  }

  # nSpacQuantile
  if (!is.numeric(args$nSpacQuantile) ||
      length(args$nSpacQuantile) != 1 ||
      args$nSpacQuantile <= 0 || args$nSpacQuantile > 1) {
    stop("`nSpacQuantile` must be a numeric scalar > 0 and <= 1.")
  }

  # reducedSpots
  if (!is.logical(args$reducedSpots) ||
      length(args$reducedSpots) != 1) {
    stop("`reducedSpots` must be a single logical value (TRUE or FALSE).")
  }

  # nReduce
  if (!is.numeric(args$nReduce) ||
      length(args$nReduce) != 1 ||
      args$nReduce != as.integer(args$nReduce)) {
    stop("`nReduce` must be a single integer value.")
  }

  invisible(TRUE)
}
.check_SpaCo_inputs <- function(args) {
  # Extract expected inputs from the list
  neighbours  <- args$neighbours
  data        <- args$data
  coordinates <- args$coordinates

  ## --- Check neighbours ---
  if (nrow(neighbours) != ncol(neighbours))
    stop("'neighbours' must be a square matrix.")

  if (!Matrix::isSymmetric(neighbours))
    stop("'neighbours' must be symmetric.")

  n <- ncol(neighbours)

  ## --- Check coordinates ---
  if (!is.matrix(coordinates) && !is.data.frame(coordinates))
    stop("'coordinates' must be a matrix or data.frame.")

  if (ncol(coordinates) != 2)
    stop("'coordinates' must have exactly 2 columns.")

  if (nrow(coordinates) != n)
    stop("'coordinates' must have the same number of rows as neighbours has columns.")

  ## --- Check data ---
  if (!is.matrix(data))
    stop("'data' must be a matrix.")

  if (nrow(data) != n)
    stop("'data' must have the same number of rows as neighbours.")

  if (any(diag(as.matrix(neighbours)) != 0))
    warning("Diagonal entries of 'neighbours' are not zero.")

  invisible(TRUE)
}
