context("mpcrossMapped validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree@selfing <- "finite"
map <- qtl::sim.map(len = rep(100, 2), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
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
test_that("Check of imputed founders matrix works",
	{
		key <- cbind(1:2, 1:2, 1:2)
		copied <- mapped

		imputedData <- matrix(1L, nrow = nLines(copied), ncol = nMarkers(copied))
		colnames(imputedData) <- markers(mapped)

		copied@geneticData[[1]]@imputed <- new("imputed", data = imputedData, key = key, map = map)
		dimnames(copied@geneticData[[1]]@imputed@data) <- dimnames(copied@geneticData[[1]]@finals)
		expect_identical(validObject(copied, complete=TRUE), TRUE)

		copied@geneticData[[1]]@imputed@data[1,1] <- 0L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied@geneticData[[1]]@imputed@data[1,1] <- 3L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied@geneticData[[1]]@imputed@data[1,1] <- 2L
		validObject(copied, complete=TRUE)

		names(copied@geneticData[[1]]@imputed@map[[1]])[1] <- "invalidMarkerName"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
rm(pedigree, map, cross, rf, mapped)
