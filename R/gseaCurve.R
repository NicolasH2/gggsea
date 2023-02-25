#' calculate coordinates for a GSEA plot
#'
#' Imports:
#' grDevices
#' dplyr
#'
#' @param rl named(!), sorted(!) vector. This ranked list's Values are the ranking metric (e.g. log2FC), names are the genes IDs. Gene IDs have to be of the same type as the ones in setList.
#' @param setlist named(!) list of character vectors. Each vector is a gene signature, each item in that vector is a gene ID (same type as the ones in rl!)
#' @param gsea data.frame with certain columns: pathway, pval, NES. The latter two will be printed on the GSEA plot.
#' @param weight number, the higher the more important are the changes at the extremes. 0: no weight, i.e. each found gene counts the same. 1: each gene counts according to its metric. 2: genes counts according to their squared matric, etc.
#' @return a data.frame with coordinates for a GSEA plot. When given as an input, geom_gsea will automatically take care. Otherwise: x and y plot the regular curve (geom_path); x, y1ticks and y2ticks plot the ticks (use geom_segment); color, x, xGradientStart, y1gradient and y2gradient for color bar (use geom_rect)
#' @details calculating the enrichment score at any given point follows standard rules. See for example https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/
#' @export
#' @examples
#' library(gggsea)
#'
#' curve <- gseaCurve(myRankedlist, mySetlist)
#'
gseaCurve <- function(rl, setlist, gsea=NULL, weight=1){

  dfList <- mapply(function(set, setname){
    if( sum(set %in% names(rl))==0 ) stop("None of the genes in the ranked list are present in the set.")

    # 0) reduce the set so that only genes are left that come up in the ranked list
    set <- set[set %in% names(rl)]

    # 1) a vector that has a number for each gene in the ranked list: 0 if not in the set, and its metric (adjusted by a defined weight) if it is in the set.
    presence <- rep(0,length(rl))
    positions <- which(names(rl) %in% set)
    presence[positions] <- abs(rl[positions])^weight

    # 2) a similar vector, except it has a 1 for every gene NOT in the list, and a 0 for everything else.
    absence <- rep(1,length(rl))
    absence[positions] <- 0

    # 3) calculate the relative cumulative increase for presence and absence
    cumPresence <- cumsum(presence)
    relPresence <- sapply(cumPresence, function(x) x/max(cumPresence)) #stepwise increase for presence
    relAbsence <- cumsum(absence) / (length(rl) - length(set)) #stepwise increase for absence

    # 4) subtract the cumulative absence from the cumulative presence to get the enrichment score
    es <- relPresence-relAbsence #enrichment score
    xcoord <- seq_along(es)

    ## the data.frame df will contain everything for the GSEAplot. For now it just contains the actual curve
    df <- data.frame(x = c(0,xcoord), y = c(0,es), set = setname, gene=c(0, names(rl)) )

    maxES <- max(df$y)
    minES <- min(df$y)
    sizeFactor <- abs(maxES - minES)
    lowestPoint <- minES - sizeFactor / 30
    df$bottomline <- lowestPoint
    df$zeroline <- median(which(rl==sort(abs(rl))[1]))

    #=======================================================
    # add statistics =======================================================
    # label will initially be empty and only be filled if gsea was provided
    statdf <- data.frame(x = 0,
                         ystat = lowestPoint+sizeFactor*.02,
                         stattext = NA )

    if(!is.null(gsea)){
      subgsea <- gsea[gsea$pathway %in% setname,]
      statdf$stattext = paste0("atop(italic(NES)==",as.character(round(subgsea$NES, 2)),
                               ",italic(p)==",      as.character(round(subgsea$pval,4)),")")
    }
    df <- merge(df, statdf, by="x", all=T) #merge the dataframe with the statistics (will add statistics coordinates and label only to the first row (x=0))

    #=======================================================
    # add ticks =======================================================
    df <- merge(df, .presenceTicks(rl, set, lowestPoint, sizeFactor), by="x", all=TRUE)
    lowestPoint <- min(df$y2ticks, na.rm=TRUE) # lowest point is changed for the color gradient

    #=======================================================
    # add color gradient =======================================================
    df <- merge(df, .colorGradient(rl, lowestPoint, sizeFactor), by="x", all=TRUE)

    return(df)

  }, set=setlist, setname=names(setlist), SIMPLIFY=FALSE)

  df <- do.call(rbind, dfList) # combine all df's (were calculated separately for each set)

  return(df)
}

#========================================
# calculate ticks
.presenceTicks <- function(rl, set, lowestPoint, sizeFactor){

  ticks <- data.frame(x = which(names(rl) %in% set),
                      y1ticks = lowestPoint - sizeFactor / 40,
                      y2ticks = lowestPoint - sizeFactor / 8,
                      hitgene = names(rl[names(rl) %in% set]))

  return(ticks)
}

#========================================
# calculate the color gradient
.colorGradient <- function(rl, lowestPoint, sizeFactor, lowcol="blue", midcol="white", highcol="red", resolution=20){

  # 1) create a data.frame that will eventually hold the plotting values. Start with a sequence from -max to +max of the ranked list's metric
  gradient <- unlist(lapply(seq(0,1,length.out=resolution/2+1), function(x) dplyr::nth(sort(abs(rl)), as.integer(length(rl)*x)) ))
  gradient <- gradient[-1]
  gradient <- sort(c(gradient,-gradient))
  gradient <- data.frame(valueMax=gradient[-1])

  # 2) add color values to the table that correspond to the metric
  colfunc1 <- grDevices::colorRampPalette(c(highcol, midcol)) #functions for getting a color ramp
  colfunc2 <- grDevices::colorRampPalette(c(midcol, lowcol))
  gradient$color <- c(colfunc1(resolution/2), colfunc2(resolution/2)[-1])

  # 3) the x columns will contain the number of genes that contain a value <= than the one in the current row. Know we know where each color starts and ends
  #(the x axis will be as long as the number of genes in the ranked list, therefore the number of genes with a metric smaller than the one that stands for a color defines how long-stretched this color will be)
  gradient$x <- sapply( gradient$valueMax, function(x) sum(rl <= x) )
  gradient <- gradient[!duplicated(gradient$x),]
  gradient$xGradientStart <- c( 1, gradient$x[-nrow(gradient)] )
  #gradient <- gradient[-nrow(gradient),]

  # 4) add y column, which will be the same for all. The y position will only be influenced by the ES values (i.e. where the curve is)
  gradient$y1gradient <- lowestPoint #multiplying by the sizefactor is necessary to keep the gradient height and position the same in every graph
  gradient$y2gradient <- lowestPoint - sizeFactor / 8

  return(gradient)
}
