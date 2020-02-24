#' @title Add noise to marker positions
#' @description Add noise to marker positions, so that no markers are co-located
#' @details Add noise to marker positions, so that no markers are located at the same position on a single chromosome. This was necessary before there was an error model implemented in the IBD genotype imputation and IBD genotype probabliity code. There is little reason to use this function now. 
#' @param map The map to add noise to. 
#' @return A copy of the input map, with noise added to genetic positions. 
#' @export
jitterMap <- function(map)
{
	eps <- 1e-3
	for(chromosome in 1:length(map))
	{
		mapPart <- map[[chromosome]]
		if(length(mapPart) > 1)
		{
			if(all(abs(mapPart - mapPart[1]) < eps)) stop("Chromosome consists of a single location, which multiple markers are mapped to. This case should be handled manually. ")
			counter <- 1
			while(TRUE)
			{
				currentLocation <- mapPart[counter]
				start <- counter
				while(counter < length(mapPart) && mapPart[counter+1] == currentLocation) counter <- counter + 1
				#So now start:counter is a range of markers which are at the same genetic location
				if(counter != start)
				{
					if(counter == length(mapPart))
					{
						totalDist <- mapPart[counter] - mapPart[start-1]
						mapPart[start:counter] <- seq(mapPart[start-1] + totalDist/(counter - start + 1), mapPart[counter], length.out = counter - start + 1)
					}
					else if(start == 1)
					{
						totalDist <- mapPart[counter+1]
						mapPart[1:counter] <- seq(0, totalDist*(counter-1)/counter, length.out = counter)
					}
					else
					{
						totalDist <- mapPart[counter+1] - mapPart[start-1]
						nMarkers <- counter - start+1
						mapPart[start:counter] <- seq(mapPart[start-1] + totalDist/(nMarkers+1), mapPart[counter+1] - totalDist / (nMarkers + 1), length.out = nMarkers)
					}
				}
				if(counter == length(mapPart)) break
				counter <- counter + 1
			}
		}
		map[[chromosome]] <- mapPart
	}
	return(map)
}
