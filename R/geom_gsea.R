# devtools::document()

# seekGSEplot <- function(rankedList, setID, setList="BP", Species=NA, gseaTable=NA,
#                         annotations=c("up","down"), plotTitle=NA,
#                         cutC=log2(1.5), gseaParam=1,
#                         vlines=T, hlines=T, lineSize=1.5, ticksSize=0.2,
#                         textSize=30, annoSize=10, colors=c("red","blue","green"),
#                         yTitle=T, xTitle=T, xText=T, ggarrangeRows=1) {
#
#   #==check=input==#
#   if(is.na(Species)){stop("!!!Error, no Species selected!!!")}
#   stats <- rankedList
#
#   ##=pathway=genes=##===================================================================
#   #==for=GO=(if)=or=customSets=(else)=#
#   if(is.character(setList)){
#     setName <- AnnotationDbi::Term(GO.db::GOTERM[AnnotationDbi::Ontology(GO.db::GOTERM)==setList]) #vector of GO-set-names, with GOIDs as vector names
#     setName <- as.character(setName[names(setName) %in% setID]) #string of GO-set-name corresponding to GO-ID given by the "setID" parameter
#     pathway <- switch(Species,
#                       "human"=as.list(org.Hs.eg.db::org.Hs.egGO2EG)[[setID]], #all gene IDs for the specific GO-ID
#                       "mouse"=as.list(org.Mm.eg.db::org.Mm.egGO2EG)[[setID]])
#   }else{
#     setName <- setID
#     pathway <- setList[[setID]]
#   }
#   print(setName)
#   setName <- gsub("^REACTOME_","",gsub("^KEGG_","",setName))
#   wordcut <- floor(nchar(setName)/2)-3
#   getspace <- FALSE
#
#   while(!getspace){
#     getspace <- {substr(setName, wordcut,wordcut) %in% c("_",",",";"," ","-",".")}
#     if(!getspace){wordcut <- wordcut+1}
#     if(wordcut==nchar(setName)){getspace <- TRUE}
#   }
#   setName <- paste0(substr(setName, 1, wordcut), "\n", substr(setName, wordcut+1, nchar(setName)))
#
#   #=====================================================================================
#   #=====================================================================================
#
#   ##=fgsea=calculation=##===============================================================
#   #==rankedList=preparation==#
#   rnk <- rank(-stats) # gives a vector: 1,2,3,4,... ranking the rankedList (actually unnecessary)
#   ord <- order(rnk) # orders the previous vector (also actually unnecessary)
#   statsAdj <- stats#[ord] # puts the rankedList into correct order (also actually unnecessary)
#   statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam) # increases the rl values by gseaParam without losing the - or + before the value
#   statsAdj <- statsAdj / max(abs(statsAdj)) # normalizing the rl values by the max
#
#   #==find=genes=from=rankedList=in=setID==#
#   pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj))))) # vector with positions genes from geneset in the ranked list
#   pathway <- sort(pathway) # vector is sorted by positions
#
#   #==fgsea=function==#
#   gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats=pathway, returnAllExtremes=T)
#
#   #=====================================================================================
#   #=====================================================================================
#
#   ##=calculations=##====================================================================
#   #==calculating=x=and=y=coordinates==#
#   bottom <- min(gseaRes$bottoms, na.rm=T)
#   top <- max(gseaRes$tops, na.rm=T)
#   xs <- as.vector(rbind(pathway-1, pathway)) #combines the vectors in an alternating fashion
#   ys <- as.vector(rbind(gseaRes$bottoms, gseaRes$tops))
#   n <- length(statsAdj)
#
#   #==gsea=line=(most=important=part)==#
#   toPlot <- data.frame(x=c(0, xs, n+1), y=c(0, ys, 0))
#   #==axis=labels==#
#   axisXtext <- c(ceiling(top*10)/10, ceiling(bottom*10)/10)
#   axisXspace <- abs(round(2*(axisXtext[1]-axisXtext[2]))/10)*2
#   axisXtext <- round(c(seq(axisXtext[2], axisXtext[1], by=ifelse(axisXspace==0, 0.1 ,axisXspace)), 0),1)
#
#   #==Title==#
#   if(!is.na(plotTitle)){ # if a plotTitle is specified, use it
#     toPlot$title = plotTitle
#   }else{
#     toPlot$title = setName # else, use the default name
#   }
#
#   #==Positions==#
#   x=y=NULL
#   sizeFactor <- ((top - bottom))
#   yStats <- c(bottom+(sizeFactor*0.01), bottom+sqrt(ggarrangeRows)*(sizeFactor*0.2*annoSize/10))
#   yTicks <- c(bottom-(sizeFactor*0.025), bottom-(sizeFactor*0.125))
#   yRibbon <- c(bottom-(sizeFactor*0.13), bottom-(sizeFactor*0.23))
#   yLabel <- bottom-(sizeFactor*0.28)
#   vline <- c(length(rankedList[rankedList>=cutC]),
#              length(rankedList)-length(rankedList[rankedList<=-cutC]),
#              length(rankedList[rankedList>=0]))
#
#   #==colorRibbon==#
#   # bar lengths for each threshold
#   barlength <- unique((cutC*c(1/4,1/2,1,2,4,4,8,8,8,8,rep(16,8))[1:(3+ceiling(max(abs(rankedList))))])) # vector with lengths for log2FC thresholds
#   barlength <- na.omit(barlength)
#   barlength <- c(rev(barlength)[2:length(barlength)],0, -barlength) # vector as above but mirrored (highest, ..., 0, ..., -highest)
#   barlength <- sapply(barlength, function(x) sum(rankedList>=x)) # vector with number of genes that have log2FC above the thresholds
#   # bar colors
#   # colfunc1 <- colorRampPalette(c(colors[1], "white")) #functions for getting a color ramp
#   # colfunc2 <- colorRampPalette(c("white", colors[2]))
#   # barcolor <- c(colfunc1(length(barlength)/2), colfunc2(length(barlength)/2))
#   barcolor <- .mycolorgradient(length(barlength), highcol=colors[2], midcol="white", lowcol=colors[1])
#   # data.frame: each row has a bar color intensity, bar length, bar start (x1) and end (x2) position and an Order-Number
#   colorRibbon <- data.frame(barcolor=barcolor,
#                             barlength=barlength,
#                             x1=c(0,barlength[1:(length(barlength)-1)]),
#                             x2=barlength)
#   # divide data frame into red and blue ribbon
#   colorRibbon1 <- colorRibbon[1:(length(barcolor)/2),]
#   colorRibbon2 <- colorRibbon[((length(barcolor)/2)+1):length(barcolor),]
#   #=====================================================================================
#   #=====================================================================================
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
