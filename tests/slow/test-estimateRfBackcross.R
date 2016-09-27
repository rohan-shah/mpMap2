context("Test recombination fraction estimation for backcross pedigree")
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
