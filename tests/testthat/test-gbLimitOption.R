context("Test option gbLimit of estimateRF")
check <- function(cross)
{
	rf1 <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
	rf2 <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE, gbLimit = 0)
	rf3 <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE, gbLimit = 61*8*1e-9)
	rf4 <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE, gbLimit = 2*61*8*1e-9)
	rf4 <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE, gbLimit = 4*61*8*1e-9)
	rf5 <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE, gbLimit = 8*61*8*1e-9)
	objects <- list(rf1, rf2, rf3, rf4, rf5)
	indices <- list(1:2, 2:3, 3:4, 4:5, c(1,5))
	lapply(indices, function(x)
		{
			y <- x[1]
			z <- x[2]
			expect_identical(objects[[y]]@rf@theta, objects[[z]]@rf@theta)
			expect_identical(objects[[y]]@rf@lod, objects[[z]]@rf@lod)
			expect_identical(objects[[y]]@rf@lkhd, objects[[z]]@rf@lkhd)
		})
}

test_that("Checking that value of gbLimit option doesn't change results for f2",
	{
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)

		check(cross)
	})
test_that("Check that huge number of markers cannot be analysed except using gbLimit option",
	{
		f2Pedigree <- f2Pedigree(10)
		map <- sim.map(len = 100, n.mar = 4000, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)

		expect_that(rf <- estimateRF(cross), throws_error())
		rf <- estimateRF(cross, gbLimit = 1)
	})
rm(check)
