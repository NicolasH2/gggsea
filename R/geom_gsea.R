# devtools::document()

#' create a venn diagram in ggplot2
#'
#' Imports:
#' ggplot2
#'
#' @inheritParams ggplot2::geom_path
#' @inheritParams ggplot2::geom_segment
#' @inheritParams ggplot2::geom_rect
#' @inheritParams ggplot2::facet_wrap
#' @import ggplot2
#'
#' @param df data.frame, calculated by the function gseaCurve
#' @param labelsize number, font size of the statistics (if gsea is provided)
#' @param prettyGSEA boolean, should some aesthetics be automatically added? Adds a 0-line, regulates y-breaks and labels the axes.
#' @param linecolor string, color of the main gsea line
#' @param linesize number, thickness of the main gsea line
#' @param tickcolor string, color of the ticks representing hit genes
#' @param ticksize number, thickness of the ticks representing hit genes
#' @return a list of ggplot layers, ready to be added to a ggplot object
#' @details uses the functions geom_gseaLine, geom_gseaTicks and geom_gseaGradient, from this package
#' @export
#' @examples
#' library(gggsea)
#' library(ggplot)
#'
#' curve <- gseaCurve(myRankedlist, mySetlist)
#' ggplot() + geom_gsea(curve) + theme_gsea()
#'
geom_gsea <- function(
  df, labelsize=4, prettyGSEA=T, ncol=NULL, nrow=NULL,
  linecolor="green", linesize=1,
  tickcolor="black", ticksize=0.5,
  ...
){
  # in case the user provides several colors/sizes, each set will get a different color
  nsets <- length(unique(df$set)) #number of sets. The reason it is calculated like this is only in case the user provides more colors (or sizes) than there are sets
  linecolor <- rep(rep(linecolor, nsets)[1:nsets], each=nrow(df)/nsets) #the linecolor vector is repeated n times (in case only 1 color was provided) and then only the first n are taken (in case too many colors were provided), which are repeated as many times as there are rows in df, so every line part gets a color
  linesize <- rep(rep(linesize, nsets)[1:nsets], each=nrow(df)/nsets)

  nticks <- unlist( lapply(unique(df$set), function(x) nrow(df[df$set %in% x & !is.na(df$y1ticks),])) )#number of ticks for each set (other than lines, there is not one tick for every row in the data.frame and the number of ticks is not the same for each set)
  tickcolor <- rep(rep(tickcolor, nsets)[1:nsets], nticks) #each color is repeated a different number of times
  ticksize <- rep(rep(ticksize, nsets)[1:nsets], nticks)

  # plot the main parts
  gseaLine <- geom_gseaLine(df,                     color=linecolor, size=linesize, ...)
  ticks <- geom_gseaTicks(df[!is.na(df$y1ticks), ], color=tickcolor, size=ticksize, ...)
  gradient <- geom_gseaGradient(df[!is.na(df$y1gradient), ], ...)
  statistics <- ggplot2::geom_label(data=df[!is.na(df$stattext),],
                                    mapping = aes(x, ystat, label=stattext ),
                                    size=labelsize, hjust=0, vjust=0 , parse=T, alpha=.7, fill="white")

  # combine all parts and add a facet_wrap
  main <- list(gseaLine, ticks, gradient, statistics,
               ggplot2::facet_wrap(~set, scale="free_y", ncol=ncol, nrow=nrow)
               )

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

