#' @import qtl
#' @import igraph
#' @import methods
#' @import utils
#' @import stats
#' @import ggplot2
#' @import graphics
#' @importFrom graphics plot
#' @importFrom stats as.dist cutree pchisq rnorm
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
