context("Test recombination fraction estimation for eight parent pedigrees with improper funnels")
test_that("Test that estimation of recombination fractions is accurate with finite generations of selfing and fully informative markers",
{
	distances <- c(1, 5, 10, 20, 50)
	pedigree <- eightParentPedigreeImproperFunnels(initialPopulationSize = 12000, selfingGenerations = 6, nSeeds = 1)
	pedigree@selfing <- "infinite"
	for(distance in distances)
	{
		map <- sim.map(len = distance, n.mar = 2, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
		#Ignore message about memory usage
		capture.output(rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), seq(0, 0.5, length.out = 150))))
		expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance = 0.06)
	}
})
test_that("Test that estimation of recombination fractions is accurate with finite generations of selfing and SNP markers",
{
	distances <- c(1, 5, 10, 20, 50)
	pedigree <- eightParentPedigreeImproperFunnels(initialPopulationSize = 30000, selfingGenerations = 10, nSeeds = 1)
	pedigree@selfing <- "infinite"
	for(distance in distances)
	{
		map <- sim.map(len = distance, n.mar = 2, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane) + multiparentSNP(keepHets = FALSE)
		#Ignore message about memory usage
		capture.output(rf <- estimateRF(cross, recombValues = c(haldaneToRf(distance), seq(0, 0.5, length.out = 150))))
		expect_equal(rfToHaldane(rf@rf@theta[1,2]), distance, tolerance = 0.1)
	}
})
