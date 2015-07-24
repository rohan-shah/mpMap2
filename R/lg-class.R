checkLG <- function(object)
{
	errors <- c()
	if(any(is.na(object@groups)))
	{
		errors <- c(errors, "Slot groups cannot contain NA values")
	}
	if(any(is.na(object@allGroups)))
	{
		errors <- c(errors, "Slot allGroups cannot contain NA values")
	}
	if(any(object@allGroups < 0))
	{
		errors <- c(errors, "Slot allGroups cannot contain negative values")
	}
	if(!all(object@groups %in% object@allGroups))
	{
		errors <- c(errors, "Only values in slot allGroups are allowed as values in slot groups")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.lg <- setClass("lg", slots = list(groups = "integer", allGroups = "integer"), validity = checkLG)