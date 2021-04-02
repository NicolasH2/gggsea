# devtools::document()

#' create a venn diagram in ggplot2
#'
#' Imports:
#' ggplot2
#'
#' @inheritParams ggplot2::geom_path
#' @inheritParams ggplot2::geom_segment
#' @inheritParams ggplot2::geom_rect
#' @import ggplot2
#'
#' @param df data.frame, calculated by the function gseaCurve
#' @param withTicks boolean, should there be ticks?
#' @param withGradient boolean, should ther be a color gradient?
#' @param prettyGSEA boolean, should some aesthetics be automatically added? Adds a 0-line, regulates y-breaks and labels the axes.
#' @return a list of ggplot layers, ready to be added to a ggplot object
#' @details uses the functions geom_gseaLine, geom_gseaTicks and geom_gseaGradient, from this package
#' @export
#' @examples
#' library(gggsea)
#' library(ggplot)
#'
#' curve <- gseaCurve()
#'
geom_gsea <- function(df, withTicks=T, withGradient=T, prettyGSEA=T, ...){

  gseaLine <- geom_gseaLine(df, ...)
  main <- list(gseaLine)

  if(withTicks){
    ticks <- geom_gseaTicks(df, ...)
    main <- c(main, ticks)
  }

  if(withGradient){
    gradient <- geom_gseaTicks(df, ...)
    main <- c(main, gradient)
  }

  if(prettyGSEA){
    maxES <- max(df$y)
    minES <- min(df$y)
    sizeFactor <- abs(maxES - minES)

    ystepsize <- signif(sizeFactor/3, digits=1) #round the sizefactor fraction to the nearest digit
    nDown <- round(minES / ystepsize) #number of ticks below 0
    nUp <- round(maxES / ystepsize) # number of ticks above 0
    breaks <- c(
      seq(from = 0, to = -nDown * ystepsize, by = ystepsize),
      seq(from = 0, to = nUp * ystepsize,    by = ystepsize)
    )

    main <- c(list(
      geom_hline(yintercept=0),
      labs(x="rank", y="enrichment score"),
      scale_y_continuous(breaks=breaks),
      main
    ))
  }

  return( main )

}






