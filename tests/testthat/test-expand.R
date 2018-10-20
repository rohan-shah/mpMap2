context("expand")
pedigree <- f2Pedigree(10)
map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
rf <- estimateRF(cross, keepLod = TRUE, keepLkhd=TRUE)

test_that("Can expand to same marker set (without a warning)",
	{
		expect_warning(expand(rf, newMarkers = markers(rf)), NA)

		lg <- formGroups(rf, groups = 1)
		expect_warning(expand(lg, newMarkers = markers(lg)), NA)
	})
test_that("Function warns when discarding recombination data",
	{
		newMarkers <- c(markers(rf), "extraMarker")
	
		expect_that(expand(rf, newMarkers = newMarkers), gives_warning("recombination data will be lost"))

		lg <- formGroups(rf, groups = 1)
		expect_that(expand(lg, newMarkers = newMarkers), gives_warning("recombination and linkage group data will be lost"))
	})

test_that("Function warns when discarding linkage group data",
	{
		newMarkers <- c(markers(rf), "extraMarker")
		lg <- formGroups(rf, groups = 1)
		expect_that(expand(lg, newMarkers = newMarkers), gives_warning("linkage group data will be lost"))
	})

test_that("New markers must contain old markers",
	{
		newMarkers <- c("extraMarker", "extraMarker2")
		expect_that(expand(rf, newMarkers = newMarkers), throws_error())

		lg <- formGroups(rf, groups = 1)
		expect_that(expand(lg, newMarkers = newMarkers), throws_error())
	})
rm(pedigree, map, cross, rf)
