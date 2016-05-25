context("estimateRf 16-parent tests")
getMap <- function(distance)
{
	map <- list("chr1" = c("a" = 0, "b" = distance))
	class(map)<- "map"
	return(map)
}
distances <- c(1, 5, 10, 20, 50)

test_that("Numerically accurate for randomly chosen funnels", 
	{
		for(intercrossingGenerations in 0:1)
		{
			for(selfingGenerations in 0:1)
			{
				pedigree <- sixteenParentPedigreeRandomFunnels(initialPopulationSize=2000, selfingGenerations = selfingGenerations, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
				pedigree@selfing <- "finite"
				for(distance in distances)
				{
					map <- getMap(distance)
					cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
					cross2 <- cross + fixedNumberOfFounderAlleles(8)
					#Ignore the warning about residual hetrozygosity
					capture.output(rf <- estimateRF(cross2, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = FALSE, keepLkhd=FALSE))
					expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=0.05)
					expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
					expect_identical(rf@rf@theta[1,1], 0)
					expect_identical(rf@rf@theta[2,2], 0)
				}
			}
		}
	})
rm(getMap, distances)
