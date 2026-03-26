utils::globalVariables(c("x", "y"))
.SpatialColors <-
  grDevices::colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
.rescale_to_range <- function(x) {
  2 * scales::rescale(x)
}
.rescale_to_spac <- function(x) {
  4 * (scales::rescale(x) - 0.5)
}

#' Plot SPaCo meta genes.
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param spac component to plot
#' @return returns a ggplot object with the meta gene expression.
#' @export
Spaco_plot <- function(SpaCoObject,
                       spac = 1) {
  suppressWarnings(.singlespacplot(SpaCoObject, spac))
}


.singlespacplot <- function(SpaCoObject,
                            spac) {
  name_arg <- paste0("spac_", spac)
  rescale_spac <- SpaCoObject@projection[, spac, drop = FALSE]
  rescale_spac[, 1] <- .rescale_to_spac(rescale_spac[, 1])
  singleplot <- ggplot2::ggplot(data = tidyr::tibble(
    tidyr::as_tibble(SpaCoObject@coordinates, rownames = "BC"),
    assign(
      paste0("spac_", spac),
      tibble::as_tibble(rescale_spac[, 1, drop = FALSE], rownames = NA)
    )
  ))  +
    ggplot2::coord_fixed() +
    ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )

  if (any((SpaCoObject@pixel_positions_list$row[1] %% 1) > 0)) {
    singleplot <-
      singleplot + ggforce::geom_regon(ggplot2::aes(
        x0 = col,
        y0 = row,
        sides = 4,
        r = 3.5,
        angle = pi / 4,
        fill = !!as.name(paste0("spac_",  spac))
      )) +
      ggplot2::scale_fill_gradientn(colours = .SpatialColors(n = 100),
                                    limits = c(-2, 2)) +
      ggplot2::scale_x_continuous(name = NULL, breaks = NULL) +
      ggplot2::scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    singleplot <-
      singleplot + ggplot2::geom_tile(ggplot2::aes(
        x = x,
        y = y,
        fill = !!as.name(paste0("spac_",  spac))
      )) +
      ggplot2::scale_fill_gradientn(colours = .SpatialColors(n = 100),
                                    limits = c(-2, 2)) +
      ggplot2::scale_x_continuous(name = NULL, breaks = NULL) +
      ggplot2::scale_y_reverse(name = NULL, breaks = NULL)
  }

  return(singleplot)
}
#' Plot denoised gene expression
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param feature Gene for which to plot denoised projection
#' @return returns a ggplot object with the denoised gene expression.
#' @export
denoised_projection_plot <- function(SpaCoObject,
                                     feature = NULL)
{
  suppressWarnings(.singledenoisedprojectionplot(SpaCoObject, feature))
}
.singledenoisedprojectionplot <-
  function(SpaCoObject, feature) {
    name_arg <- feature
    rescaled_denoised <-
      SpaCoObject@denoised[, feature, drop = FALSE]
    rescaled_denoised[, 1] <-
      .rescale_to_range(rescaled_denoised[, 1])
    singleplot <- ggplot2::ggplot(data = tidyr::tibble(
      tidyr::as_tibble(SpaCoObject@coordinates, rownames = "BC"),
      tidyr::as_tibble(rescaled_denoised[, feature  , drop = FALSE],
                       rownames = NA)
    )) +
      ggplot2::coord_fixed() +
      ggplot2::theme_linedraw(base_size = 10) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top"
      )

    if (any((SpaCoObject@pixel_positions_list$row[1] %% 1) > 0)) {
      singleplot <-
        singleplot + ggforce::geom_regon(ggplot2::aes(
          x0 = x,
          y0 = y,
          sides = 4,
          r = 3.5,
          angle = pi / 4,
          fill = !!as.symbol(paste0(feature))
        )) +
        ggplot2::scale_fill_gradientn(colours = .SpatialColors(n = 100),
                                      limits = c(0, 2)) +
        ggplot2::scale_x_continuous(name = NULL, breaks = NULL) +
        ggplot2::scale_y_reverse(name = NULL, breaks = NULL)
    } else {
      singleplot <-
        singleplot + ggplot2::geom_tile(ggplot2::aes(
          x = x,
          y = y,
          fill = !!as.symbol(paste0(feature))
        )) +
        ggplot2::scale_fill_gradientn(colours = .SpatialColors(n = 100),
                                      limits = c(0, 2)) +
        ggplot2::scale_x_reverse(name = NULL, breaks = NULL) +
        ggplot2::scale_y_continuous(name = NULL, breaks = NULL)
    }
    return(singleplot)
  }
#' Plot gene expression
#'
#' @param SpaCoObject SpacoObject with computed projections
#' @param feature Gene to plot
#' @return returns a ggplot object with gene expression.
#' @export
feature_plot <-
  function(SpaCoObject,
           feature)
  {
    suppressWarnings(.singledataplot(SpaCoObject, feature))
  }


.singledataplot <- function(SpaCoObject, feature) {
  rescaled_data <- SpaCoObject@data[, feature,
                                    drop = FALSE]
  rescaled_data[, 1] <- .rescale_to_range(rescaled_data[, 1])
  singleplot <-
    ggplot2::ggplot(data = tidyr::tibble(
      tidyr::as_tibble(SpaCoObject@coordinates,
                       rownames = "BC"),
      tidyr::as_tibble(rescaled_data[, feature,
                                     drop = FALSE],
                       rownames = NA)
    )) +
    ggplot2::coord_fixed() +
    ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )
  if (any((SpaCoObject@pixel_positions_list$row[1] %% 1) > 0)) {
    singleplot <-
      singleplot + ggforce::geom_regon(ggplot2::aes(
        x0 = x,
        y0 = y,
        sides = 4,
        r = 3.5,
        angle = pi / 4,
        fill = !!as.symbol(paste0(feature))
      )) +
      ggplot2::scale_fill_gradientn(colours = .SpatialColors(n = 100),
                                    limits = c(0, 2)) +
      ggplot2::scale_x_continuous(name = NULL, breaks = NULL) +
      ggplot2::scale_y_reverse(name = NULL, breaks = NULL)
  } else {
    singleplot <-
      singleplot + ggplot2::geom_tile(ggplot2::aes(
        x = x,
        y = y,
        fill = !!as.symbol(paste0(feature))
      )) +
      ggplot2::scale_fill_gradientn(colours = .SpatialColors(n = 100),
                                    limits = c(0, 2)) +
      ggplot2::scale_x_reverse(name = NULL, breaks = NULL) +
      ggplot2::scale_y_continuous(name = NULL, breaks = NULL)

  }
  return(singleplot)
}
