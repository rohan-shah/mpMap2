combineGenotypes <- function(finals, hetData)
{
	isHetDataArgument(hetData)
	isNumericMatrixArgument(finals)
	if(ncol(finals) != 2*length(hetData))
	{
		stop("Inconsistent dimensions for inputs")
	}
	if(any(colnames(finals) != names(hetData)))
	{
		stop("Inconsistent names for inputs")
	}
	.Call("combineGenotypes", finals, hetData, PACKAGE = "mpMap2")
}