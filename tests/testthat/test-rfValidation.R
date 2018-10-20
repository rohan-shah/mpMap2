context("rf validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree@selfing <- "finite"
map <- qtl::sim.map(len = rep(100, 1), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
rf <- estimateRF(cross)

test_that("Simulated rf passes validation",
	{
		expect_identical(validObject(rf, complete=TRUE), TRUE)
	})
test_that("Slot rf@theta must have storage.mode double",
	{
		#Can't be a character
		copy <- rf
		expect_that(copy@rf@theta <- as.character(copy@rf@theta), throws_error())

		#Or an integer
		copy <- rf
		expect_that(copy@rf@theta <- as.integer(copy@rf@theta), throws_error())
	})
test_that("Slot rf@lod can be NULL",
	{
		copy <- rf
		copy@rf@lod <- NULL
		expect_identical(validObject(copy, complete=TRUE), TRUE)
	})
test_that("Slot rf@lkhd can be NULL",
	{
		copy <- rf
		copy@rf@lkhd <- NULL
		expect_identical(validObject(copy, complete=TRUE), TRUE)
	})
test_that("Slot rf@lod must have storage.mode double if it is not NULL",
	{
		copy <- rf
		expect_that(copy@rf@lod <- as.character(copy@rf@lod), throws_error())

		copy <- rf
		expect_that(copy@rf@lod <- as.integer(copy@rf@lod), throws_error())
	})
test_that("Slot rf@lkhd must have storage.mode double if it is not NULL",
	{
		copy <- rf
		expect_that(copy@rf@lkhd <- as.character(copy@rf@lkhd), throws_error())

		copy <- rf
		expect_that(copy@rf@lkhd <- as.integer(copy@rf@lkhd), throws_error())
	})
rm(pedigree, map, cross, rf)
