context("estimateRf generic tests")
map <- list("chr1" = c("a" = 0, "b" = 100))
class(map)<- "map"

test_that("Checking that we correctly reject a line not in the pedigree",
	{
		pedigree <- f2Pedigree(100)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		rownames(cross@geneticData[[1]]@finals)[1] <- "abcdef"
		expect_that(rf <- estimateRF(cross, recombValues = (0:100)/200, keepLod = TRUE, keepLkhd=TRUE), throws_error())
	})
test_that("Checking that we require 0 and 0.5 in the recombination values", 
	{
		pedigree <- f2Pedigree(100)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_that(rf <- estimateRF(cross, recombValues = c(0, 0.1), keepLod = TRUE, keepLkhd=TRUE), throws_error())
		expect_that(rf <- estimateRF(cross, recombValues = c(0.1, 0.5), keepLod = TRUE, keepLkhd=TRUE), throws_error())
		expect_error(rf <- estimateRF(cross, recombValues = c(0, 0.1, 0.5), keepLod = TRUE, keepLkhd=TRUE), NA)
	})
test_that("Checking that we must have less than 255 possible levels",
	{
		pedigree <- f2Pedigree(100)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		rf <- estimateRF(cross, recombValues = (0:253)/506, keepLod = TRUE, keepLkhd=TRUE)
		expect_that(rf <- estimateRF(cross, recombValues = (0:254)/508, keepLod = TRUE, keepLkhd=TRUE), throws_error())
	})
test_that("Check that a warning is generated if there are hetrozygotes but we're assuming infinite selfing",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 101)
		cross <- simulateMPCross(map = map, pedigree = pedigree, seed = 1, mapFunction = haldane)
		cross@geneticData[[1]]@pedigree@selfing <- "infinite"
		expect_warning(rf <- estimateRF(cross), "Input data had heterozygotes but was analysed assuming infinite selfing")
})
rm(map)
test_that("Check that we can compute recombination fractions if there are no lines",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 11, include.x = FALSE)
		cross <- simulateMPCross(map = map, pedigree = pedigree, seed = 1, mapFunction = haldane)
		cross <- subset(cross, lines = as.character(c()))
		rf <- estimateRF(cross)

		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
		map <- qtl::sim.map(len = 100, n.mar = 11, include.x = FALSE)
		cross <- simulateMPCross(map = map, pedigree = pedigree, seed = 1, mapFunction = haldane)
		cross <- subset(cross, lines = as.character(c()))
		expect_error(rf <- estimateRF(cross), NA)
	})
