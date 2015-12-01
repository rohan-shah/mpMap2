context("mpcrossMapped validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree@selfing <- "auto"
map <- sim.map(len = rep(100, 2), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
rf <- estimateRF(cross)
mapped <- mpcrossMapped(rf, map = map)

test_that("Simulated mpcrossMapped passes validation",
	{
		expect_identical(validObject(mapped, complete=TRUE), TRUE)
	})
test_that("Cross object must have same number of markers as map",
	{
		copied <- mapped
		copied@map[[1]]<- head(copied@map[[1]], -1)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Marker of cross object must agree with markers of map",
	{
		copied <- mapped
		names(copied@map[[1]])[1] <- "invalidMarkerName"
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- mapped
		names(copied@map[[1]])[1:2] <- names(copied@map[[1]])[2:1]
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
rm(pedigree, map, cross, rf, mapped)
