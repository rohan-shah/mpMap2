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
		pedigree <- f2Pedigree(5000)
		for(distance in distances)
		{
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			rf <- estimateRF(cross, recombValues = (0:100)/200, keepLod = TRUE, keepLkhd=TRUE)
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=0.06)
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})

test_that("Numerically accurate for a RIL design", 
	{
		set.seed(2)
		pedigree <- rilPedigree(20000, 8)
		for(distance in distances)
		{
			map <- getMap(distance)
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			#There is a warning that we're analysing a population with hetrozygotes as if hetrozygotes were impossible. 
			capture.output(rf <- estimateRF(cross, recombValues = (0:100)/200, keepLod = TRUE, keepLkhd=TRUE))
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance=0.06)
			expect_identical(rf@rf@theta[1,2], rf@rf@theta[2,1])
			expect_identical(rf@rf@theta[1,1], 0)
			expect_identical(rf@rf@theta[2,2], 0)
		}
	})
