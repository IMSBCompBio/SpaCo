.SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
.rescale_to_range <- function(x) {
  2*rescale(x)
}
.rescale_to_spac <- function(x) {
  4*(rescale(x) - 0.5)
}

#' Plot SPaCo meta genes.
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param spac component to plot
#'
#' @return returns a ggplot object with the meta gene expression.
#' @export
#' @import ggplot2
#' @import ggforce
#' @import tidyr
#' @import rcartocolor
#' @import patchwork
#' @import dplyr
#' @import scales
Spaco_plot <- function(SpaCoObject, spac = 1, ncol = NULL, combine = TRUE, dataset = "dataset1") {
  require(scales)
  require(patchwork) # Ensure patchwork is also loaded for combining plots

  # Initialize the plots list
  plots <- vector("list", length = length(spac))

  # Determine the class of the object
  object_class <- class(SpaCoObject)[1] # Assuming the first class name is the most specific one

  for (i in seq_along(spac)) {
    # Choose the function based on object class
    if (object_class == "mergedSpaCoObject") {
      # Call .mergedsinglespacplot for mergedSpaCoObject
      plots[[i]] <- suppressWarnings(.mergedsinglespacplot(SpaCoObject, i = i, spac = spac, dataset = dataset))
    } else if (object_class == "SpaCoObject") {
      # Call .singlespacplot for SpaCoObject
      plots[[i]] <- suppressWarnings(.singlespacplot(SpaCoObject, i = i, spac = spac))
    } else {
      stop("Unsupported object class: ", object_class)
    }
  }

  # Combine the plots if requested
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol, guides = "auto")
  }

  return(plots)
}


.singlespacplot <- function(SpaCoObject, i = i, spac = spac) {
  name_arg <- paste0("spac_", spac[i])
  rescale_spac <- SpaCoObject@projection[, spac[i], drop = FALSE]
  rescale_spac[,1] <- .rescale_to_spac(rescale_spac[,1])
  singleplot <- ggplot(data = tibble(
  tidyr::as_tibble(SpaCoObject@pixel_positions_list, rownames = "BC"),
  assign(paste0("spac_", i), tibble::as_tibble(rescale_spac[, 1, drop = FALSE], rownames = NA))))  +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

  if (any((SpaCoObject@pixel_positions_list$imagerow[1] %% 1) > 0)) {
    singleplot <- singleplot + ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                                                       sides = 4, r = 3.5, angle = pi / 4, fill = !!as.name(paste0("spac_",  spac[i]))))+
      scale_fill_gradientn(colours = .SpatialColors(n=100),
                           limits = c(-2, 2))+scale_x_continuous(name = NULL, breaks = NULL) +
      scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    singleplot <- singleplot + geom_tile(aes(x = imagecol, y = imagerow, fill = !!as.name(paste0("spac_",  spac[i]))))+
      scale_fill_gradient(low = "white", high = "black") + coord_flip()+ scale_x_reverse(name = NULL, breaks = NULL)

  }

  return(singleplot)
}

.mergedsinglespacplot <- function(SpaCoObject, i = 1, spac = spac, dataset = "dataset2") {
  pattern <- paste0("(.*)", dataset, "$")

  # Extract all rownames from the projections list for the specified dataset
  all_rownames <- rownames(SpaCoObject@projection[["merged"]])

  # Use grep to find rownames that match the pattern
  matching_rownames <- grep(pattern, all_rownames, value = TRUE)

  name_arg <- paste0("spac_", spac[i])
  rescale_spac <- SpaCoObject@projection$merged[matching_rownames, spac[i], drop = FALSE]
  rescale_spac[,1] <- .rescale_to_spac(rescale_spac[,1])
  mergedsingleplot <- ggplot(data = tibble(
    tidyr::as_tibble(SpaCoObject@pixel_positions_list[[dataset]], rownames = "BC"),
    assign(paste0("spac_", i), tibble::as_tibble(rescale_spac[, 1, drop = FALSE], rownames = NA))))  +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

  if (any((SpaCoObject@pixel_positions_list[[dataset]]$imagerow[1] %% 1) > 0)) {
    mergedsingleplot <- mergedsingleplot + ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                                                       sides = 4, r = 3.5, angle = pi / 4, fill = !!as.name(paste0("spac_",  spac[i]))))+
      scale_fill_gradientn(colours = .SpatialColors(n=100),
                           limits = c(-2, 2))+scale_x_continuous(name = NULL, breaks = NULL) +
      scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    mergedsingleplot <- mergedsingleplot + geom_tile(aes(x = imagecol, y = imagerow, fill = !!as.name(paste0("spac_",  spac[i]))))+
      scale_fill_gradient(low = "white", high = "black") + coord_flip()+ scale_x_reverse(name = NULL, breaks = NULL)

  }

  return(mergedsingleplot)
}






#' Plot denoised gene expression
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param spac component to plot
#'
#' @return returns a ggplot object with the denoised gene expression.
#' @export
#'
#' @import ggplot2
#' @import ggforce
#' @import tidyr
#' @import rcartocolor
#' @import patchwork
#' @import scales
denoised_projection_plot <- function(SpaCoObject, features = NULL, ncol = NULL, combine = TRUE)
{
  plots <- vector(
    mode = "list",
    length = length(features))
  for (i in 1:length(features)) {
    plots[[i]] <- suppressWarnings(.singledenoisedprojectionplot(SpaCoObject, i = i, features))
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol, guides = "auto")
  }
  return(plots)
}


.singledenoisedprojectionplot <- function(SpaCoObject, i = i, features) {
  require(scales)
  name_arg <- features[i]
  rescaled_denoised <- SpaCoObject@denoised[ ,features[i]  , drop = FALSE]
  rescaled_denoised[,1] <- .rescale_to_range(rescaled_denoised[,1])
  singleplot <- ggplot(data = tibble(
    tidyr::as_tibble(SpaCoObject@pixel_positions_list, rownames = "BC"),
    as_tibble(rescaled_denoised[ ,features[i]  , drop = FALSE],rownames=NA))) +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

  if (any((SpaCoObject@pixel_positions_list$imagerow[1] %% 1) > 0)) {
    singleplot <- singleplot + ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                                                       sides = 4, r = 3.5, angle = pi / 4, fill = !!as.symbol(paste0(features[i]))))+
      scale_fill_gradientn(colours = .SpatialColors(n=100),
                           limits = c(0, 2)) +scale_x_continuous(name = NULL, breaks = NULL) +
      scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    singleplot <- singleplot + geom_tile(aes(x = imagecol, y = imagerow,fill = !!as.symbol(paste0(features[i])))) +
      scale_fill_gradient(low = "white", high = "black") + coord_flip()+ scale_x_reverse(name = NULL, breaks = NULL)+scale_y_continuous(name = NULL, breaks = NULL)

  }

  return(singleplot)
}






#' Plot gene expression
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param features gene to plot
#'
#' @return returns a ggplot object with gene expression.
#' @export
#'
#' @import ggplot2
#' @import ggforce
#' @import tidyr
#' @import rcartocolor
#' @import patchwork
#' @import scales
feature_plot <- function(SpaCoObject, features = NULL, ncol = NULL, combine = TRUE)
{
  plots <- vector(
    mode = "list",
    length = length(features))
  for (i in 1:length(features)) {
    plots[[i]] <- suppressWarnings(.singledataplot(SpaCoObject, i = i, features))
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol, guides = "auto")
  }
  return(plots)
}


.singledataplot <- function(SpaCoObject, i = i, features) {
  require(scales)
  name_arg <- features[i]
  rescaled_data <- SpaCoObject@data[ ,features[i]  , drop = FALSE]
  rescaled_data[,1] <- .rescale_to_range(rescaled_data[,1])
  singleplot <- ggplot(data = tibble(
    tidyr::as_tibble(SpaCoObject@pixel_positions_list, rownames = "BC"),
    as_tibble(rescaled_data[ ,features[i]  , drop = FALSE],rownames=NA))) +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

  if (any((SpaCoObject@pixel_positions_list$imagerow[1] %% 1) > 0)) {
    singleplot <- singleplot + ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                                                       sides = 4, r = 3.5, angle = pi / 4, fill = !!as.symbol(paste0(features[i]))))+
      scale_fill_gradientn(colours = .SpatialColors(n=100),
                           limits = c(0, 2)) +scale_x_continuous(name = NULL, breaks = NULL) +
      scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    singleplot <- singleplot + geom_tile(aes(x = imagecol, y = imagerow,fill = !!as.symbol(paste0(features[i])))) +
      scale_fill_gradient(low = "white", high = "black") + coord_flip()+ scale_x_reverse(name = NULL, breaks = NULL)+scale_y_continuous(name = NULL, breaks = NULL)

  }

  return(singleplot)
}














