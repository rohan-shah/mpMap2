#' @import qtl
#' @importFrom igraph plot.igraph V E "E<-" graph.empty vertices edges layout.reingold.tilford
#' @import methods
#' @import utils
#' @import ggplot2
#' @import graphics
#' @importFrom progress progress_bar
#' @importFrom graphics plot
#' @importFrom stats as.dist cutree pchisq rnorm cor manova na.omit
#' @importFrom utils combn head tail
#' @importFrom pryr address
#' @importFrom fastcluster hclust
#' @importFrom nnls nnls
#' @importFrom ggplot2 ggplot
#' @importFrom Heatplus heatmap_2
#' @importFrom RColorBrewer brewer.pal
#' @exportClass pedigreeGraph
#' @exportClass probabilities
#' @exportMethod subset 
#' @exportClass mpcross
#' @exportClass pedigreeGraph
#' @exportMethod imputationMap
#' @exportMethod plot
#' @exportMethod flatImputationMapNames
#' @importClassesFrom Matrix index dspMatrix dppMatrix
#' @importFrom methods setClass
#' @importFrom jsonlite toJSON
#' @useDynLib mpMap2
NULL
