context("Test backcross pedigree")
test_that("Test that genetic data is as expected",
	{
		map <- sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- backcrossPedigree(5000)
		cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
		proportions <- table(cross@geneticData[[1]]@finals)/ length(cross@geneticData[[1]]@finals)
		#Expect that there are only two genotypes
		expect_identical(length(proportions), 2L)
		#Expect that proportions are close to half. 
		expect_true(all(abs(proportions - 0.5) < 0.015))
	})
test_that("Test that estimation of recombination fractions is accurate",
	{
		distances <- c(1, 5, 10, 20, 50)
		pedigree <- backcrossPedigree(100000)
		for(distance in distances)
		{
			map <- sim.map(len = distance, n.mar = 2, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
			cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
			rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), seq(0, 0.5, length.out = 200)))
			expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance = 0.03)
		}
	})
