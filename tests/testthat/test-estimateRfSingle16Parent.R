context("estimateRF vs estimateRFSingle, 16-parent tests")

test_that("Using randomly chosen funnels", 
	{
		distances <- c(1, 5, 10, 20, 50)
		for(intercrossingGenerations in 0:1)
		{
			pedigree <- sixteenParentPedigreeRandomFunnels(initialPopulationSize=100, selfingGenerations = 3, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
			pedigree@selfing <- "finite"
			for(distance in distances)
			{
				map <- list("chr1" = c("a" = 0, "b" = distance))
				class(map)<- "map"
				cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
				cross2 <- cross + removeHets()
				#Ignore the warning about residual hetrozygosity
				capture.output(rfSingle <- estimateRFSingleDesign(cross2, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
				capture.output(rf <- estimateRF(cross2, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
				expect_identical(rfSingle, rf)
				expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
				expect_identical(rf@rf@theta[1,1], 0)
				expect_identical(rf@rf@theta[2,2], 0)
			}
		}
	})
