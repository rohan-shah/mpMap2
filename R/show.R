#' @include mpcross-class.R
#' @include print.R
setMethod(f = "show", signature = "mpcrossRF", definition = function(object) print(object))
setMethod(f = "show", signature = "mpcross", definition = function(object) print(object))
setMethod(f = "show", signature = "mpcrossLG", definition = function(object) print(object))
