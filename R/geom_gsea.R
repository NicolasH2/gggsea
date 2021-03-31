# devtools::document()

geom_gsea <- function(rl, set, weight=1, showTicks=T, showGradient=T, prettyGSEA=T, ...){

  df <- .gseaCurve(rl, set, weight)

  main <- list( ggplot2::geom_path(ggplot2::aes(x=x, y=y), data=df, color="green") )

  if(showTicks){
    ticksdf <- .presenceTicks(rl, df)
    ticks <- ggplot2::geom_segment(data=ticksdf, mapping=ggplot2::aes(x=x, y=y1, xend=x, yend=y2), size=0.2)

    main <- c(main, ticks)
  }

  if(showGradient){
    gradientdf <- .colorRibbon(rl, df)
    gradient <- ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), data=gradientdf, fill=gradientdf$color)

    main <- c(main, gradient)
  }

  if(prettyGSEA){
    main <- c(
      geom_hline(yintercept=0),
    main)
  }

  return( main )

}


theme_gsea <- function(){
  mytheme <- theme(

  )
}

.gseaCurve <- function(rl, set, weight=1){
  # 0) reduce the set so that only genes are left that come up in the ranked list
  set <- set[set %in% names(rl)]

  # 1) a vector that has a number for each gene in the ranked list: 0 if not in the set, and its metric (adjusted by a defined weight) if it is in the set.
  presence <- ifelse(names(rl) %in% set, rl^weight, 0)

  # 2) a similar vector, except it has a 1 for every gene NOT in the list, and a 0 for everything else.
  absence <- ifelse(presence>0,0,1)

  # 3) calculate the relative cumulative increase for presence and absence
  cumPresence <- cumsum(presence)
  relPresence <- sapply(cumPresence, function(x) x/max(cumPresence)) #stepwise increase for presence
  relAbsence <- cumsum(absence) / (length(rl) - length(set)) #stepwise increase for absence

  # 4) subtract the cumulative absence from the cumulative presence to get the enrichment score
  es <- resPresence-resAbsence #enrichment score
  xcoord <- seq_along(es)

  df <- data.frame(x = xcoord, y = es)
  return(df)
}

.presenceTicks <- function(rl, dfLine){
  minES <- min(dfLine$y)
  maxES <- max(dfLine$y)
  sizeFactor <- maxES - minES

  ticks <- data.frame(x = which(names(rl) %in% set),
                      y1 = minES-(sizeFactor*0.025),
                      y2 = minES-(sizeFactor*0.125))
  return(ticks)
}

.colorRibbon <- function(rl, dfLine, lowcol="blue", midcol="white", highcol="red", resolution=100){

  # 1) create a data.frame that will eventually hold the plotting values. Start with a sequence from -max to +max of the ranked list's metric
  limit <- max(abs(rl))
  gradient <- data.frame( valueMax = seq(-limit, limit, length.out=resolution) )[-1,,drop=FALSE]

  # 2) add color values to the table that correspond to the metric
  colfunc1 <- grDevices::colorRampPalette(c(highcol, midcol)) #functions for getting a color ramp
  colfunc2 <- grDevices::colorRampPalette(c(midcol, lowcol))
  gradient$color <- c(colfunc1(resolution/2), colfunc2(resolution/2)[-1])

  # 3) the x columns will contain the number of genes that contain a value <= than the one in the current row. Know we know where each color starts and ends
  #(the x axis will be as long as the number of genes in the ranked list, therefore the number of genes with a metric smaller than the one that stands for a color defines how long-stretched this color will be)
  gradient$x1 <- sapply( gradient$valueMax, function(x) sum(rl <= x) )
  gradient <- gradient[!duplicated(gradient$x1),]
  gradient$x2 <- c( 1, gradient$x1[-nrow(gradient)] )
  gradient <- gradient[-nrow(gradient),]

  # 4) add y column, which will be the same for all. The y position will only be influenced by the ES values (i.e. where the curve is)
  minES <- min(dfLine$y)
  maxES <- max(dfLine$y)
  sizeFactor <- maxES - minES
  gradient$y1 <- minES - sizeFactor*0.13
  gradient$y2 <- minES - sizeFactor*0.23

  return(gradient)
}

#
#   ##=GSEA=plot=##=======================================================================
#   plot <- ggplot2::ggplot(toPlot, ggplot2::aes(x=x, y=y)) +
#     #==main=line==#
#     ggplot2::geom_point(color=colors[3], size=0.1) +
#     ggplot2::geom_line(color=colors[3], size=lineSize) +
#
#     #==hit=segments==#
#     ggplot2::geom_segment(data=data.frame(x=pathway), mapping=ggplot2::aes(x=x, y=yTicks[1], xend=x, yend=yTicks[2]), size=ticksSize) +
#
#     #==colorRibbon=left=(red)=and=right=(blue)==#
#     ggplot2::geom_rect(data=colorRibbon1, ggplot2::aes(xmin=x1, xmax=x2, ymin=yRibbon[1], ymax=yRibbon[2]), fill=colorRibbon1$barcolor) +
#     ggplot2::geom_rect(data=colorRibbon2, ggplot2::aes(xmin=x1, xmax=x2, ymin=yRibbon[1], ymax=yRibbon[2]), fill=colorRibbon2$barcolor) +
#     ggplot2::annotate(geom="text", size=annoSize, y=yLabel, x=(length(stats)/400), label=annotations[1], color=colors[1], hjust=0) +
#     ggplot2::annotate(geom="text", size=annoSize, y=yLabel, x=length(rankedList)-(length(stats)/400), label=annotations[2], color=colors[2], hjust=1) +
#
#     #==theme=settings==#
#     ggplot2::geom_hline(yintercept=0, colour="gray") +
#     ggplot2::theme_bw() +
#     ggplot2::theme(panel.border=ggplot2::element_rect(fill=NA, size=0.5),
#                    panel.grid=ggplot2::element_blank(), legend.position="none",
#                    axis.text.x=ggplot2::element_text(angle=45, hjust=1),
#                    axis.text= ggplot2::element_text(size=textSize),
#                    axis.title=ggplot2::element_text(size=textSize),
#                    strip.text=ggplot2::element_text(size=textSize)) +
#
#     #==axis=titles=and=labels==#
#     ggplot2::scale_y_continuous(name="enrichment score", labels=axisXtext, breaks=axisXtext) +
#     ggplot2::labs(x="rank") +
#     #==plot=title==#
#     ggplot2::facet_grid(. ~ title)
#
#   #==additional=lines==#
#   if(vlines){plot <- plot + ggplot2::geom_vline(xintercept = vline, linetype="dotted", color="gray20")}
#   if(hlines){plot <- plot + ggplot2::geom_hline(yintercept=c(top, bottom), colour="red", linetype="dotted")}
#
#   #==statistics===#
#   if(is.data.frame(gseaTable)){
#     pval <- round(gseaTable[which(gseaTable[,"GO_ID"] %in% setID),"pval"], 4)
#     NES <- round(gseaTable[which(gseaTable[,"GO_ID"] %in% setID),"NES"], 2)
#     plot <- plot +
#       ggplot2::annotate(geom="text",size=annoSize, x=(length(stats)/10), y=yStats[1], label=substitute(paste(italic("NES"), " = ",NES)),hjust=0,vjust=0) +
#       ggplot2::annotate(geom="text",size=annoSize, x=(length(stats)/10), y=yStats[2], label=substitute(paste(italic("p"), " = ",pval)),hjust=0,vjust=0)
#   }
#
#   #==make=axis=text=and=title=blank=(user=choice)==#
#   if(!yTitle){plot <- plot + ggplot2::theme(axis.title.y=element_blank())}
#   if(!xTitle){plot <- plot + ggplot2::theme(axis.title.x=element_blank())}
#   if(!xText){plot <- plot + ggplot2::theme(axis.text.x=element_blank())}
#
#   #==output==#
#   plot
# }
