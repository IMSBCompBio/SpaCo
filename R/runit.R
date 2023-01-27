data_dir="~/Liver_Organoid_project/scRNA_spatial_tissue/VX03_A006200136/cellranger/141895/outs"
slice = "~/Liver_Organoid_project/scRNA_spatial_tissue/VX03_A006200136/cellranger/141895/outs/spatial"
filename ="filtered_feature_bc_matrix.h5"

data_dir="~/Spca_example_data"
slice="~/Spca_example_data/spatial"
filename ="V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5"



SpaCoObject <-read_10x_for_spaco(data_dir = data_dir,slice = slice,filename =filename,variable_features_n =3333 )

SpaCoObject <- SCA_function(SpaCoObject, PC_criterion = "number",
                         PC_value = 50, orthogonalResult = FALSE)



apply(ONB_OriginalBasis, 2, function(x) {
  Expression <- t(ONB_OriginalBasis[,x]) %*% data_centered})


Expression_Locus <- apply(ONB_OriginalBasis, 2, function(x) {
  Expression <- t(ONB_OriginalBasis[,x]) %*% data_centered
  Locus <- rownames(data)
  return(data.frame(Expression, Locus))
})

ONB_OriginalBasis_small <- ONB_OriginalBasis[,1:4]

m <- apply(ONB_OriginalBasis_small, 2, "*", data_centered)



GeneData <- data.frame(Expression = t(Spaco %*% data_centered_t),
                       Locus = rownames(data))


test <- apply(SpaCoObject@spacs, 2, function(x) t(x %*% t(scale(SpaCoObject@data, scale = FALSE))))


Spaco <- ONB_OriginalBasis[,1,drop=F]
data_centered <- t(scale(data, scale = FALSE))
t(Spaco %*% data_centered)
dim(data_centered)

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

aaaaaa <- plotSpacoProjectionFunction(Spaco=SpaCoObject@spacs[,1], data = SpaCoObject@data,
                                             Coordinate_Data = SpaCoObject@coordinates, pointSize = 5,
                                             customScale = TRUE)

aaaaaa





result <- matrix(nrow = 1595,ncol = ncol(SpaCoObject@spacs))
test_loop <- for(i in 1:ncol(SpaCoObject@spacs)){
  result[,i] <- t(SpaCoObject@spacs[,i] %*% t(scale(SpaCoObject@data, scale = FALSE)))
  return(result)
}



all.equal(test,result)
