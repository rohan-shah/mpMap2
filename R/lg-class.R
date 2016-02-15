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
	if(!is.null(object@imputedTheta))
	{
		if(!is.list(object@imputedTheta))
		{
			errors <- c(errors, "If slot imputedTheta is not null, it must be a list of rawSymmetricMatrix objects")
			return(errors)
		}
		if(length(object@imputedTheta) != length(object@allGroups))
		{
			errors <- c(errors, "Slot imputedTheta had the wrong length")
			return(errors)
		}
		if(any(unlist(lapply(object@imputedTheta, class)) != "rawSymmetricMatrix"))
		{
			errors <- c(errors, "If slot imputedTheta is not null, it must be a list of rawSymmetricMatrix objects")
			return(errors)
		}
		if(!identical(names(object@imputedTheta), as.character(object@allGroups)))
		{
			errors <- c(errors, "If slot imputedTheta is not null, its names must match the values in slot allGroups")
			return(errors)
		}
		groupCounts <- sapply(object@allGroups, function(x) sum(object@groups == x))
		imputedThetaLengths <- unlist(lapply(object@imputedTheta, function(x) length(x@data)))
		if(any(imputedThetaLengths != groupCounts*(groupCounts + 1)/2))
		{
			errors <- c(errors, "Slot imputedTheta contained objects with the wrong length")
			return(errors)
		}
		#Check that markers are correct
		correctMarkers <- sapply(1:length(object@allGroups), function(x)
			{
				group <- object@allGroups[x]
				imputedMarkers <- object@imputedTheta[[x]]@markers
				groupMarkers <- names(which(object@groups == group))
				return(identical(groupMarkers, imputedMarkers))
			})
		if(any(!correctMarkers))
		{
			errors <- c(errors, "Markers in matrices contained in object@imputedTheta were inconsistent with those in slot object@groups")
			return(errors)
		}
		#Check that the levels are all the same
		if(length(object@imputedTheta) > 0)
		{
			levelsFirst <- object@imputedTheta[[1]]@levels
			for(i in 1:length(object@imputedTheta))
			{
				if(!identical(object@imputedTheta[[i]]@levels, levelsFirst))
				{
					errors <- c(errors, "Slot levels must be the same for all objects contained in slot imputedTheta")
					return(errors)
				}
			}
		}
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
setClassUnion("listOrNULL", c("list", "NULL"))
.lg <- setClass("lg", slots = list(groups = "integer", allGroups = "integer", imputedTheta = "listOrNULL"), validity = checkLG)
