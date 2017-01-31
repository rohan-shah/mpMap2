context("Test backcross pedigree")
test_that("Test that genetic data is as expected",
	{
		map <- sim.map(len = 500, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- backcrossPedigree(1000)
		cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
		proportions <- table(cross@geneticData[[1]]@finals)/ length(cross@geneticData[[1]]@finals)
		#Expect that there are only two genotypes
		expect_identical(length(proportions), 2L)
		#Expect that proportions are close to half. 
		expect_true(all(abs(proportions - 0.5) < 0.02))
	})
