context("estimateRf 8-parent tests")
getMap <- function(distance)
{
	map <- list("chr1" = c("a" = 0, "b" = distance))
	class(map)<- "map"
	return(map)
}
distances <- c(1, 5, 10, 20, 50)

test_that("Numerically accurate with randomly chosen funnels", 
	{
		testFunc <- function(pedigree)
		{
			for(distance in distances)
			{
				map <- getMap(distance)
				cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
				#Ignore the warning about residual heterozygosity
				capture.output(rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
				expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=0.02)
				expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
				expect_identical(rf@rf@theta[1,1], 0)
				expect_identical(rf@rf@theta[2,2], 0)
			}
		}
		for(intercrossingGenerations in 0:1)
		{
			for(selfingGenerations in 0:1)
			{
				pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize=100000, selfingGenerations = selfingGenerations, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
				pedigree@selfing <- "finite"
				testFunc(pedigree)
				pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize=100000, selfingGenerations = selfingGenerations, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
				pedigree@selfing <- "finite"
				testFunc(pedigree)
			}
		}
	})
rm(getMap, distances)
