context("estimateRf 4-parent tests")
getMap <- function(distance)
{
	map <- list("chr1" = c("a" = 0, "b" = distance))
	class(map)<- "map"
	return(map)
}
distances <- c(1, 5, 10, 20, 50)

test_that("Numerically accurate for a single marker", 
	{
		testFunction <- function(pedigree)
		{
			for(distance in distances)
			{
				map <- getMap(distance)
				cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
				rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE)
				expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=0.015)
				expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
				expect_identical(rf@rf@theta[1,1], 0)
				expect_identical(rf@rf@theta[2,2], 0)
			}
		}
		for(intercrossingGenerations in 0:1)
		{
			for(selfingGenerations in 0:1)
			{
				pedigreeSingleFunnel <- fourParentPedigreeSingleFunnel(initialPopulationSize=50000, selfingGenerations = selfingGenerations, nSeeds = 1, intercrossingGenerations = intercrossingGenerations)
				pedigreeSingleFunnel@selfing <- "finite"
				testFunction(pedigreeSingleFunnel)

				pedigreeRandomFunnels <- fourParentPedigreeRandomFunnels(initialPopulationSize=50000, selfingGenerations = selfingGenerations, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
				pedigreeRandomFunnels@selfing <- "finite"
				testFunction(pedigreeRandomFunnels)

			}
		}
		
	})
rm(getMap, distances)
