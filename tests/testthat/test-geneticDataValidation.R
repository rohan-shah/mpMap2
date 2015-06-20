context("geneticData validation")

#Set up a basic mpcross object
twoWayGeneticDataSetUp <- function()
{
	env <- sys.frame(sys.parent(1))
	pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
	map <- sim.map(len = rep(100, 1), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
	cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
	return(cross@geneticData[[1]])
}
test_that("geneticData rejects floating point values", 
	{
		geneticData <- twoWayGeneticDataSetUp()
		expect_identical(validObject(geneticData), TRUE)
	})
