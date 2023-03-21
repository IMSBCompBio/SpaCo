plotGeneFunction <- function(data, geneName, Coordinate_Data = coordinate_Data,
                             pointSize = 5, customScale = TRUE)
{
  GenePlotData <- Coordinate_Data
  GenePlotData$Locus <- rownames(Coordinate_Data)
  GeneData <- data.frame(Expression = data[,which(colnames(data) == geneName)],
                         Locus = rownames(data))
  Spaco1PlotData <- left_join(GenePlotData, GeneData, by = "Locus")

  Spaco1Plot <- ggplot(Spaco1PlotData, aes(x = col, y = row, color = Expression)) +
    geom_point(size = pointSize, pch = 18) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("") +
    ggtitle(paste("Expression pattern for gene ", geneName, sep = ""))
  if(customScale)
  {
    Spaco1Plot <- Spaco1Plot + scale_colour_gradient2(
      low = "blue",
      mid = "lightgray",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "colour"
    )
  }
  return(Spaco1Plot)
}
plotSimGeneFunction <- function(gene, data = data_clean,
                                Coordinate_Data = coordinate_Data, pointSize = 5,
                                customScale = TRUE)
{
  GenePlotData <- Coordinate_Data
  GenePlotData$Locus <- rownames(Coordinate_Data)
  GeneData <- data.frame(Expression = gene,
                         Locus = rownames(data))
  Spaco1PlotData <- left_join(GenePlotData, GeneData, by = "Locus")

  Spaco1Plot <- ggplot(Spaco1PlotData, aes(x = col, y = row, color = Expression)) +
    geom_point(size = pointSize, pch = 18) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("")
  if(customScale)
  {
    Spaco1Plot <- Spaco1Plot +
      scale_colour_gradient2(
        low = "blue",
        mid = "lightgray",
        high = "red",
        midpoint = 0,
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
  }
  return(Spaco1Plot)
}
plotSpacoProjectionFunction <- function(Spaco, data = data_clean,
                                Coordinate_Data = coordinate_Data, pointSize = 5,
                                customScale = TRUE)
{
  data_centered <- t(scale(data, scale = FALSE))
  GenePlotData <- Coordinate_Data
  GenePlotData$Locus <- rownames(Coordinate_Data)
  GeneData <- data.frame(Expression = t(Spaco %*% data_centered),
                         Locus = rownames(data))
  Spaco1PlotData <- left_join(GenePlotData, GeneData, by = "Locus")

  Spaco1Plot <- ggplot(Spaco1PlotData, aes(x = col, y = row, color = Expression)) +
    geom_point(size = pointSize, pch = 18) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("")
  if(customScale)
  {
    Spaco1Plot <- Spaco1Plot +
      scale_colour_gradient2(
        low = "blue",
        mid = "lightgray",
        high = "red",
        midpoint = 0,
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
  }
  return(Spaco1Plot)
}

plotSpacoProjectionFunction_objectified <- function(SpaCoObject, pointSize = 5,
                                        customScale = TRUE,spac=1)
{

  data = SpaCoObject@data
  Coordinate_Data = SpaCoObject@coordinates
  GenePlotData <- Coordinate_Data
  GenePlotData$Locus <- rownames(Coordinate_Data)
  GeneData <- data.frame(Expression = SpaCoObject@projection[,spac],
                         Locus = rownames(data))
  Spaco1PlotData <- left_join(GenePlotData, GeneData, by = "Locus")

  Spaco1Plot <- ggplot(Spaco1PlotData, aes(x = col, y = row, color = Expression)) +
    geom_point(size = pointSize, pch = 15) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("")
  if(customScale)
  {
    Spaco1Plot <- Spaco1Plot +
      scale_colour_gradient2(
        low = "blue",
        mid = "lightgray",
        high = "red",
        midpoint = 0,
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      )
  }
  return(Spaco1Plot)
}


plotSpacoProjectionFunction_objectified_raster <- function(SpaCoObject, pointSize = 5,
                                                    customScale = TRUE,spac=1)
{

  data = SpaCoObject@data
  Coordinate_Data = SpaCoObject@coordinates
  GenePlotData <- Coordinate_Data
  GenePlotData$Locus <- rownames(Coordinate_Data)
  GeneData <- data.frame(Expression = SpaCoObject@projection[,spac],
                         Locus = rownames(data))
  Spaco1PlotData <- left_join(GenePlotData, GeneData, by = "Locus")

  Spaco1Plot <- ggplot(Spaco1PlotData, aes(x = col, y = row, color = Expression)) +
    geom_raster(aes(fill = Expression)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("")
  if(customScale)
  {
    Spaco1Plot <- Spaco1Plot +
      scale_fill_gradientn(colours = viridis::inferno(n=20))
  }
  return(Spaco1Plot)
}

#' Title
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param spac Spatial component to plot
#'
#' @return
#' @export
#'
#' @examples
#' @import ggplot2
#' @import ggforce
#' @import tidyr
#' @import rcartocolor
#' @import patchwork
Spaco_plot <- function(SpaCoObject,spac = 1, ncol = NULL, combine = TRUE)
{
 plots <- vector(
    mode = "list",
    length = length(spac))
 for (i in spac) {
  plots[[i]] <- .singlespacplot(SpaCoObject, i = i)
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
  ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                          sides = 4, r = 3.5, angle = pi / 4, fill = !!as.name(paste0("spac_", i)))) +
  scale_x_continuous(name = NULL, breaks = NULL) +
  scale_y_reverse(name = NULL, breaks = NULL) +
  scale_fill_gradientn(name = name_arg,colours = .SpatialColors(n=100)) +
  #rcartocolor::scale_fill_carto_c(name = name_arg,
                               #   type = "diverging", palette = "TealRose") +
  coord_fixed() +
  theme_linedraw(base_size = 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")
  return(singleplot)

}

#' Title
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param spac Spatial component to plot
#'
#' @return
#' @export
#'
#' @examples
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
    plots[[i]] <- .singlesmoothedprojectionplot(SpaCoObject, i = i, features)
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
    as_tibble(SpaCoObject@smoothed[ ,features[i]  , drop = FALSE],rownames=NA)))  +
    ggforce::geom_regon(aes(x0 = imagecol, y0 = imagerow,
                            sides = 4,r =1, angle = pi / 4, fill = !!as.symbol(paste0(features[i])))) +
    scale_x_continuous(name = NULL, breaks = NULL) +
    scale_y_reverse(name = NULL, breaks = NULL) +
    scale_fill_gradientn(name = name_arg,colours = .SpatialColors(n=100)) +
    #rcartocolor::scale_fill_carto_c(name = name_arg,
    #   type = "diverging", palette = "TealRose") +
    coord_fixed() +
    theme_linedraw(base_size = 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top")
  return(singleplot)

}



