context("estimateRf two marker tests")
getMap <- function(distance)
{
	map <- list("chr1" = c("a" = 0, "b" = distance))
	class(map)<- "map"
	return(map)
}
distances <- c(1, 5, 10, 20, 50)

test_that("Numerically accurate for an F2 design", 
	{
		set.seed(1)
		pedigree <- f2Pedigree(3000)
		for(distance in distances)
		{
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			rf <- estimateRF(cross, recombValues = (1:100)/200, keepLod = TRUE, keepLkhd=TRUE)
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=distance/80)
		}
	})

test_that("Numerically accurate for a RIL design", 
	{
		set.seed(2)
		pedigree <- rilPedigree(3000, 8)
		for(distance in distances)
		{
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			rf <- estimateRF(cross, recombValues = (1:100)/200, keepLod = TRUE, keepLkhd=TRUE)
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=distance/80)
		}
	})