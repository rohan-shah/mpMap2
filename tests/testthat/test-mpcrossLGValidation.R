context("mpcrossLG validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree@selfing <- "auto"
map <- sim.map(len = rep(100, 1), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
rf <- estimateRF(cross)
grouped <- formGroups(rf, groups = 2, clusterBy = "theta")

test_that("Simulated mpcrossLG passes validation",
	{
		expect_identical(validObject(grouped, complete=TRUE), TRUE)
	})
test_that("Names of slot lg@groups must be the marker names",
	{
		copied <- grouped
		names(copied@lg@groups)[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- grouped
		names(copied@lg@groups) <- NULL
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- grouped
		names(copied@lg@groups)[1] <- "invalidMarkerName"
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- grouped
		names(copied@lg@groups)[1:2] <- names(copied@lg@groups)[2:1]
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
rm(pedigree, map, cross, rf, grouped)
