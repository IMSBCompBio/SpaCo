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

  Spaco1Plot <- ggplot(Spaco1PlotData, aes(x = array_col, y = array_row, color = Expression)) +
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
