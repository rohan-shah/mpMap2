#' @title Get or set number of threads for OpenMP
#' @rdname openmp
#' @export
#' @description Get or set the number of threads for OpenMP
#' @details Some functions in mpMap2 are parallelised. Depending on the number of cores available, and the type of workload, it may be advantageous to turn parallelisation on or off, by setting the number of OpenMP threads appropriately. Setting the number of threads to 1 turns parallelisation off
#' 
#' In particular, for small examples on a computer with a large number of threads, parallelisation may result in a huge decrease in performance. 
#'
#' This function returns an error if the package was not compiled with OpenMP. 
#' @param num New number of threads for OpenMP
#' @return None
omp_set_num_threads <- function(num)
{
	.Call("omp_set_num_threads", num, PACKAGE = "mpMap2")
}
#' @export
#' @rdname openmp
#' @return The number of threads for OpenMP
omp_get_num_threads <- function()
{
	.Call("omp_get_num_threads", PACKAGE = "mpMap2")
}
