context("generateGenotypes")
test_that("Reject out of range recombination fractions",
{
	adjacentRecombination <- seq(0, 0.5, by = 0.01)
	markerNames <- as.character(1:(length(adjacentRecombination)+1))
	pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 10, selfingGenerations = 1, intercrossingGenerations = 0)
	genotypes <- .Call("generateGenotypes", adjacentRecombination, markerNames, pedigree, PACKAGE="mpMap2")

	markerNames <- as.character(1:(length(adjacentRecombination)))
	expect_error(genotypes <- .Call("generateGenotypes", adjacentRecombination, markerNames, pedigree, PACKAGE="mpMap2"))

	adjacentRecombination <- seq(-0.1, 0.5, by = 0.01)
	expect_error(genotypes <- .Call("generateGenotypes", adjacentRecombination, markerNames, pedigree, PACKAGE="mpMap2"), "All recombination fractions must be between 0 and 0.5", class = "std::runtime_error")

	adjacentRecombination <- seq(0, 0.6, by = 0.01)
	expect_error(genotypes <- .Call("generateGenotypes", adjacentRecombination, markerNames, pedigree, PACKAGE="mpMap2"), "All recombination fractions must be between 0 and 0.5", class = "std::runtime_error")
})
