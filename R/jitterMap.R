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
						mapPart[start:counter] <- seq(mapPart[start-1] + totalDist/nMarkers, mapPart[counter], length.out = nMarkers)
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
