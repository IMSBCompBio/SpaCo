library(ggplot2)
.SpatialColors <- colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
brainCoordinateData <- readRDS("C:/LocalStuff/data_for_david/data_for_david/brain_coordinates.RDS")
spots_intersect <- readRDS("X:/SPACO Projekt/Revision Benchmarking/SVG/data/intersect.RDS")
brainCoordinateData <- brainCoordinateData[spots_intersect,]
plotGeneHelper <- function(gene, pointSize = 1.5)
{
  pltData <- cbind(brainCoordinateData, t(gene))
  colnames(pltData) <- c("row", "col", "Expression")
  Spaco1Plot <- ggplot(pltData, aes(x = col, y = row, color = Expression)) +
    geom_point(size = pointSize, pch = 18) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    xlab("") + ylab("") +
      scale_color_gradientn(colours = .SpatialColors(n=100))
  return(Spaco1Plot)
}
