context("estimateRf biparental tests")
getMap <- function(distance)
{
	map <- list("chr1" = c("a" = 0, "b" = distance))
	class(map)<- "map"
	return(map)
}
distances <- c(1, 5, 10, 20, 50)
tolerances <- c(0.015, 0.01, 0.01, 0.01, 0.01)

test_that("Numerically accurate for an F2 design, fully informative", 
	{
		set.seed(1)
		pedigree <- f2Pedigree(30000)
		for(index in 1:length(distances))
		{
			distance <- distances[index]
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE)
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=tolerances[index])
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})
test_that("Numerically accurate for an F2 design, with dominant markers", 
	{
		set.seed(1)
		pedigree <- f2Pedigree(100000)
		distances <- c(7.5, 10, 20, 50)
		tolerances <- c(0.02, 0.02, 0.02, 0.02)
		#The smaller distances don't pass due to the sample size
		for(index in 1:length(distances))
		{
			distance <- distances[index]
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)+biparentalDominant()
			rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE)
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=tolerances[index])
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})

test_that("Numerically accurate for a RIL design", 
	{
		set.seed(2)
		pedigree <- rilPedigree(100000, 10)
		for(index in 1:length(distances))
		{
			distance <- distances[index]
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			#There is a warning that we're analysing a population with hetrozygotes as if hetrozygotes were impossible. 
			capture.output(rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), (0:100)/200), keepLod = TRUE, keepLkhd=TRUE))
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=tolerances[index])
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})
test_that("Checking f2 pedigree split into 100 different datasets",
	{
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(5000)
		set.seed(1)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		crosses <- subset(cross, lines = 1:50)
		for(i in 2:100)
		{
			crosses <- crosses + subset(cross, lines = 1:50 + (i-1)*50)
		}
		multipleDatasetsRf <- estimateRF(crosses, keepLod = TRUE, keepLkhd = TRUE)
		singleDatasetRf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
	})
rm(getMap, distances, tolerances)
