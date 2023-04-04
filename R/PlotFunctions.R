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
Spaco_plot <- function(SpaCoObject, spac = 1, ncol = NULL, combine = TRUE)
{
 plots <- vector(
    mode = "list",
    length = length(spac))
 for (i in spac) {
  plots[[i]] <- suppressWarnings(.singlespacplot(SpaCoObject, i = i))
       }
 if (combine) {
  plots <- patchwork::wrap_plots(plots, ncol = ncol, guides = "auto")
      }
  return(plots)
}

.SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))

.singlespacplot <- function(SpaCoObject, i = i) {
  name_arg <- paste0("spac_", i)
  singleplot <- ggplot(data = tibble(
  tidyr::as_tibble(SpaCoObject@pixel_positions_list, rownames = "BC"),
  assign(paste0("spac_", i), tibble::as_tibble(SpaCoObject@projection[, i, drop = FALSE], rownames = NA))))  +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

  if (any((SpaCoObject@pixel_positions_list$imagerow[1] %% 1) > 0)) {
    singleplot <- singleplot + ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                                                       sides = 4, r = 3.5, angle = pi / 4, fill = !!as.name(paste0("spac_", i))))+
      scale_fill_gradientn(name = name_arg,colours = .SpatialColors(n=100))+scale_x_continuous(name = NULL, breaks = NULL) +
      scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    singleplot <- singleplot + geom_tile(aes(x = imagecol, y = imagerow, fill = !!as.name(paste0("spac_", i))))+
      scale_fill_gradient(low = "white", high = "black") + coord_flip()+ scale_x_reverse(name = NULL, breaks = NULL)

  }

  return(singleplot)
}


#' Plot smoothed gene expression
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param spac component to plot
#'
#' @return returns a ggplot object with the smoothed gene expression.
#' @export
#'
#' @import ggplot2
#' @import ggforce
#' @import tidyr
#' @import rcartocolor
#' @import patchwork
smoothed_projection_plot <- function(SpaCoObject, features = NULL, ncol = NULL, combine = TRUE)
{
  plots <- vector(
    mode = "list",
    length = length(features))
  for (i in 1:length(features)) {
    plots[[i]] <- suppressWarnings(.singlesmoothedprojectionplot(SpaCoObject, i = i, features))
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol, guides = "auto")
  }
  return(plots)
}


.singlesmoothedprojectionplot <- function(SpaCoObject, i = i, features) {
  name_arg <- features[i]
  singleplot <- ggplot(data = tibble(
    tidyr::as_tibble(SpaCoObject@pixel_positions_list, rownames = "BC"),
    as_tibble(SpaCoObject@smoothed[ ,features[i]  , drop = FALSE],rownames=NA))) +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

  if (any((SpaCoObject@pixel_positions_list$imagerow[1] %% 1) > 0)) {
    singleplot <- singleplot + ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                                                       sides = 4, r = 3.5, angle = pi / 4, fill = !!as.symbol(paste0(features[i]))))+
      scale_fill_gradientn(name = name_arg,colours = .SpatialColors(n=100)) +scale_x_continuous(name = NULL, breaks = NULL) +
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
feature_plot <- function(SpaCoObject, features = NULL, ncol = NULL, combine = TRUE)
{
  plots <- vector(
    mode = "list",
    length = length(features))
  for (i in 1:length(features)) {
    plots[[i]] <- suppressWarnings(.singlesmoothedprojectionplot(SpaCoObject, i = i, features))
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol, guides = "auto")
  }
  return(plots)
}


.singlesmoothedprojectionplot <- function(SpaCoObject, i = i, features) {
  name_arg <- features[i]
  singleplot <- ggplot(data = tibble(
    tidyr::as_tibble(SpaCoObject@pixel_positions_list, rownames = "BC"),
    as_tibble(SpaCoObject@data[ ,features[i]  , drop = FALSE],rownames=NA))) +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")

  if (any((SpaCoObject@pixel_positions_list$imagerow[1] %% 1) > 0)) {
    singleplot <- singleplot + ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                                                       sides = 4, r = 3.5, angle = pi / 4, fill = !!as.symbol(paste0(features[i]))))+
      scale_fill_gradientn(name = name_arg,colours = .SpatialColors(n=100)) +scale_x_continuous(name = NULL, breaks = NULL) +
      scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    singleplot <- singleplot + geom_tile(aes(x = imagecol, y = imagerow,fill = !!as.symbol(paste0(features[i])))) +
      scale_fill_gradient(low = "white", high = "black") + coord_flip()+ scale_x_reverse(name = NULL, breaks = NULL)+scale_y_continuous(name = NULL, breaks = NULL)

  }

  return(singleplot)
}














