#' applies a GSEA-friendly theme to a ggplot object
#'
#' Imports:
#' ggplot2
#'
#' @import ggplot2
#'
#' @param textsize number, sets the size for every text in the plot, except those that are defined extra (e.g. annotate)
#' @return a ggplot layer that can be added to a ggplot object
#' @details sets options for the panel.grid, panel.background, panel.border, strip.background, strip.text, text, axis.text, axis.title and plot.title
#' @export
#' @examples
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_hm()
#'
theme_gsea <- function(textsize = 14){
  mytheme <- ggplot2::theme(
    panel.grid =       element_blank(),
    panel.background = element_blank(),
    panel.border =     element_rect(size=.5, fill=NA),
    strip.background = element_rect(colour = "black"),
    strip.text =       element_text(size=textsize, color="black"),
    text =             element_text(size=textsize, color="black"),
    axis.text.y =      element_text(size=textsize, color="black"),
    axis.text.x =      element_text(size=textsize, color="black", angle=45, hjust=1),
    axis.title =       element_text(size=textsize, color="black"),
    plot.title =       element_text(size=textsize, color="black"),

  )
}
