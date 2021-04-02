#' GSEA ticks for ggplot2
#'
#' Imports:
#' ggplot2
#'
#' @inheritParams geom_gsea
#' @inheritParams ggplot2::geom_segment
#' @import ggplot2
#'
#' @return a ggplot layer that can be added to a ggplot object
#' @details uses the ggplot2::layer function, with geom="segment"
#' @export
#' @examples
#' library(gggsea)
#' library(ggplot)
#'
#' curve <- gseaCurve()
#'
geom_gseaTicks <- function(df, ...){

  gseaTicks <- ggplot2::layer(
    data = df,
    mapping = ggplot2::aes(x=x, y=y1ticks, xend=x, yend=y2ticks),
    geom = "segment",
    stat = "identity",
    position = "identity",
    show.legend = FALSE,
    inherit.aes = TRUE,
    params=list(...)
  )
  userInput <- names(as.list(match.call())) #get all parameters set by the user
  if(!"size" %in% userInput) gseaTicks$aes_params[["size"]] <- 0.2

  return(gseaTicks)
}
