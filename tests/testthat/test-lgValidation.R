context("lg validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree@selfing <- "auto"
map <- sim.map(len = rep(100, 1), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
rf <- estimateRF(cross)
grouped <- formGroups(rf, groups=2, clusterBy = "theta")
lg <- grouped@lg

test_that("Simulated lg passes validation",
	{
		expect_identical(validObject(lg, complete=TRUE), TRUE)
	})
test_that("Slot groups cannot contain NA values",
	{
		copied <- lg
		copied@groups[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot allGroups cannot contain NA values",
	{
		copied <- lg
		copied@allGroups[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot allGroups must contain nonnegative values",
	{
		copied <- lg
		copied@groups[1] <- -1L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot groups must contain values that are also in allGroups",
	{
		copied <- lg
		copied@groups[1] <- 99999L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
rm(pedigree, map, cross, rf, grouped, lg)
