context("estimateRF vs estimateRFSingle, 8-parent tests")

test_that("Test actual maps",
{
	testFunc <- function(pedigree, map)
	{
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + multiparentSNP(keepHets = TRUE)
		suppressWarnings(capture.output(rf1 <- estimateRF(cross)))
		suppressWarnings(capture.output(rf2 <- mpMap2:::estimateRFSingleDesign(cross)))
		expect_identical(rf1, rf2)
	}
	map1 <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	map2 <- qtl::sim.map(len = c(100, 100), n.mar = 11, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	maps <- list(map1, map2)

	pedigree1 <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
	pedigree1@selfing <- "infinite"
	pedigree2 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 10, nSeeds = 1, intercrossingGenerations = 0)
	pedigree2@selfing <- "infinite"
	pedigrees <- list(pedigree1, pedigree2)
	for(map in maps)
	{
		for(pedigree in pedigrees)
		{
			testFunc(pedigree, map)
		}
	}
})
test_that("Using randomly chosen funnels, with finite generations of selfing", 
	{
		distances <- c(1, 5)
		for(intercrossingGenerations in 0:1)
		{
			for(selfingGenerations in 0:1)
			{
				pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize=100, selfingGenerations = selfingGenerations, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
				pedigree@selfing <- "finite"
				for(distance in distances)
				{
					map <- list("chr1" = c("a" = 0, "b" = distance))
					class(map)<- "map"
					cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
					#Ignore the warning about residual hetrozygosity
					capture.output(rfSingle <- estimateRFSingleDesign(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
					capture.output(rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
					expect_identical(rfSingle, rf)
					expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
					expect_identical(rf@rf@theta[1,1], 0)
					expect_identical(rf@rf@theta[2,2], 0)
				}
			}
		}
	})

test_that("Using randomly chosen funnels, with infinite generations of selfing", 
	{
		distances <- c(1, 5)
		for(intercrossingGenerations in 0:1)
		{
			pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize=100, selfingGenerations = 6, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
			pedigree@selfing <- "infinite"
			for(distance in distances)
			{
				map <- list("chr1" = c("a" = 0, "b" = distance))
				class(map)<- "map"
				cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
				#Ignore the warning about residual hetrozygosity
				capture.output(rfSingle <- estimateRFSingleDesign(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
				capture.output(rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
				expect_identical(rfSingle, rf)
				expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
				expect_identical(rf@rf@theta[1,1], 0)
				expect_identical(rf@rf@theta[2,2], 0)
			}
		}
	})
