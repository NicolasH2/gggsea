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
geom_gsea <- function(df, prettyGSEA=T, ...){

  gseaLine <- geom_gseaLine(df, ...)
  ticks <- geom_gseaTicks(df[!is.na(df$y1ticks), ], ...)
  gradient <- geom_gseaGradient(df[!is.na(df$y1gradient), ], ...)

  main <- list(gseaLine, ticks, gradient)

  # beautify the graph
  if(prettyGSEA){

    break_fun <- function(y){
      maxY <- max(y) #important: y is not the y in the aes, but rather the maximum and minimum y values on the entire graph
      minY <- min(y)
      sizeFactor <- abs(maxY-minY)
      minY <- minY + sizeFactor * 0.22
      ystepsize <- signif(sizeFactor/4, digits=1) #round the sizefactor fraction to the nearest digit

      nDown <- ceiling(minY / ystepsize) #number of ticks below 0
      nUp <- round(maxY / ystepsize) # number of ticks above 0
      breaks <- c(
        seq(from = nDown * ystepsize, to = 0, by = ystepsize),
        seq(from = 0, to = nUp * ystepsize,    by = ystepsize)
      )

      return(breaks)
    }

    main <- c(list(
      geom_hline( mapping = aes(yintercept=bottomline), data = df[!duplicated(df$set),] ),
      geom_hline(yintercept=0),
      labs(x="rank", y="enrichment score"),
      scale_y_continuous(breaks=break_fun),
      main
    ))
  }

  return( main )

}






