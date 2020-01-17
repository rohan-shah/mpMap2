.rcpp_error_recorder <- function(e){  
    invisible( .Call( "rcpp_error_recorder", e ) )
}

.warningsEnv <- new.env()
.warningsEnv$warnings <- character()

.rcpp_warning_recorder <- function(w){
    .warningsEnv$warnings <- append(.warningsEnv$warnings, w$message)
    invokeRestart("muffleWarning")
}

.rcpp_collect_warnings <- function() {
    warnings <- .warningsEnv$warnings
    .warningsEnv$warnings <- character()
    warnings
}


