context("estimateRf biparental tests")
getMap <- function(distance)
{
	map <- list("chr1" = c("a" = 0, "b" = distance))
	class(map)<- "map"
	return(map)
}
distances <- c(1, 5, 10, 20, 50)
tolerances <- c(0.015, 0.01, 0.01, 0.01, 0.03)

test_that("Numerically accurate for an F2 design, fully informative", 
	{
		pedigree <- f2Pedigree(100000)
		for(index in 1:length(distances))
		{
			distance <- distances[index]
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = FALSE, keepLkhd=FALSE)
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=tolerances[index])
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})
test_that("Numerically accurate for an F2 design, with dominant markers", 
	{
		pedigree <- f2Pedigree(200000)
		distances <- c(7.5, 10, 20, 50)
		tolerances <- c(0.02, 0.02, 0.02, 0.025)
		#The smaller distances don't pass due to the sample size
		for(index in 1:length(distances))
		{
			distance <- distances[index]
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)+biparentalDominant()
			rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = FALSE, keepLkhd=FALSE)
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=tolerances[index])
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})

test_that("Numerically accurate for a RIL design", 
	{
		pedigree <- rilPedigree(100000, 10)
		for(index in 1:length(distances))
		{
			distance <- distances[index]
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			#There is a warning that we're analysing a population with heterozygotes as if heterozygotes were impossible. 
			capture.output(suppressWarnings(rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = FALSE, keepLkhd=FALSE)))
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=tolerances[index])
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})
test_that("Checking that f2 pedigree split into two sets of markers gives the same non-NA rf estimates",
	{
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)

		#In this case we're double-counting the data for markers 4-8, which explains the strange regions for the likelihood and lod. 
		cross1 <- subset(cross, markers = 1:8)
		cross2 <- subset(cross, markers  = 4:11)
		
		#Change dataset2 line names
		cross2@geneticData[[1]]@pedigree@lineNames <- paste0(cross2@geneticData[[1]]@pedigree@lineNames, ",2")
		rownames(cross2@geneticData[[1]]@finals) <- paste0(rownames(cross2@geneticData[[1]]@finals), ",2")
		rownames(cross2@geneticData[[1]]@founders) <- paste0(rownames(cross2@geneticData[[1]]@founders), ",2")

		combined <- cross1 + cross2
		combinedRf <- estimateRF(combined, keepLod = TRUE, keepLkhd = TRUE)
		rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		#Reorganise the markers of combinedRf
		combinedRf <- subset(combinedRf, markers = markers(rf))
		#Everything that isn't NA (0xFF) in combinedRf should be equal to the values in rf
		expect_that(all((combinedRf@rf@theta@data == as.raw(0xff)) | (combinedRf@rf@theta@data == rf@rf@theta@data)), is_true())
		
		combinedRfLod <- as(combinedRf@rf@lod, "matrix")
		rfLod <- as(rf@rf@lod, "matrix")

		combinedRfLkhd <- as(combinedRf@rf@lkhd, "matrix")
		rfLkhd <- as(rf@rf@lkhd, "matrix")
		expect_that(all(abs(combinedRfLod[1:8,1:3] - rfLod[1:8, 1:3]) < 1e-5, na.rm=TRUE), is_true())
		expect_that(all(abs(combinedRfLod[9:11,4:11] - rfLod[9:11, 4:11]) < 1e-5, na.rm=TRUE), is_true())
		expect_that(all(abs(combinedRfLkhd[1:8,1:3] - rfLkhd[1:8, 1:3]) < 1e-5, na.rm=TRUE), is_true())
		expect_that(all(abs(combinedRfLkhd[9:11,4:11] - rfLkhd[9:11,4:11]) < 1e-5, na.rm=TRUE), is_true())
	})
test_that("Checking that f2 pedigree split into two sets of markers gives the same non-NA rf estimates",
	{
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)

		#In this case we avoid double-counting in the centre, but count only half the lines on 1:3 and 9:11. 
		cross1 <- subset(cross, markers = 1:8)
		cross2 <- subset(cross, markers  = 4:11)

		cross1 <- subset(cross1, lines = 1:250)
		cross2 <- subset(cross2, lines = 251:500)

		combined <- cross1 + cross2
		combinedRf <- estimateRF(combined, keepLod = TRUE, keepLkhd = TRUE)
		rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		#Reorganise the markers of combinedRf
		combinedRf <- subset(combinedRf, markers = markers(rf))
		#Everything that isn't NA (0xFF) in combinedRf should be equal to the values in rf
		expect_identical(rf@rf@theta[4:8, 4:8], combinedRf@rf@theta[4:8,4:8])
		
		combinedRfLod <- as(combinedRf@rf@lod, "matrix")
		rfLod <- as(rf@rf@lod, "matrix")

		combinedRfLkhd <- as(combinedRf@rf@lkhd, "matrix")
		rfLkhd <- as(rf@rf@lkhd, "matrix")
		expect_that(all(abs(combinedRfLod[4:8,4:8] - rfLod[4:8, 4:8]) < 1e-5, na.rm=TRUE), is_true())
		expect_that(all(abs(combinedRfLkhd[4:8,4:8] - rfLkhd[4:8, 4:8]) < 1e-5, na.rm=TRUE), is_true())

	})
rm(getMap, distances, tolerances)
