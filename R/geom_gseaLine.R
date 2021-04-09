#' GSEA line for ggplot2
#'
#' Imports:
#' ggplot2
#'
#' @inheritParams geom_gsea
#' @inheritParams ggplot2::geom_path
#' @import ggplot2
#'
#' @return a ggplot layer that can be added to a ggplot object
#' @details uses the ggplot2::layer function, with geom="path"
#' @export
#' @examples
#' library(gggsea)
#' library(ggplot)
#'
#' curve <- gseaCurve(myRankedlist, mySetlist)
#' ggplot() + geom_gseaLine(curve) + theme_gsea()
#'
geom_gseaLine <- function(df, ...){

  gseaLine <- ggplot2::layer(
    data = df,
    mapping = ggplot2::aes(x=x, y=y),
    geom = "path",
    stat = "identity",
    position = "identity",
    show.legend = FALSE,
    inherit.aes = TRUE,
    params=list(lineend="round", ...)
  )
  userInput <- names(as.list(match.call())) #get all parameters set by the user
  if(!any(c("color","colour") %in% userInput)) gseaLine$aes_params[["colour"]] <- "green"
  if(!"size" %in% userInput) gseaLine$aes_params[["size"]] <- 1

  return(gseaLine)
}
