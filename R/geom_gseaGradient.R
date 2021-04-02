#' GSEA color gradient for ggplot2
#'
#' Imports:
#' ggplot2
#'
#' @inheritParams geom_gsea
#' @inheritParams ggplot2::geom_rect
#' @import ggplot2
#'
#' @return a ggplot layer that can be added to a ggplot object
#' @details uses the ggplot2::layer function, with geom="rect"
#' @export
#' @examples
#' library(gggsea)
#' library(ggplot)
#'
#' curve <- gseaCurve()
#'
geom_gseaGradient <- function(df, ...){

  gseaGradient <- ggplot2::layer(
    data = df,
    mapping = ggplot2::aes(xmin=xGradientStart, ymin=y1gradient, xmax=x, ymax=y2gradient),
    geom = "rect",
    stat = "identity",
    position = "identity",
    show.legend = FALSE,
    inherit.aes = TRUE,
    params=list(...)
  )
  userInput <- names(as.list(match.call())) #get all parameters set by the user
  if(!"size" %in% userInput) gseaGradient$aes_params[["size"]] <- 0.2
  if(!"fill" %in% userInput) gseaGradient$aes_params[["fill"]] <- df$color

  return(gseaGradient)
}

