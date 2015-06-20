context("map validation")

map <- sim.map(len = rep(100, 3), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)

test_that("Generated map passes validation", 
	{
		#An actual map passes validation
		expect_identical({map <- new("map", map); validObject(map, complete=T)}, TRUE)
	})
test_that("Map must have chromosomes", 
	{
		#If there are no chromosomes then throw an error
		expect_that({map <- new("map", list())}, throws_error())
	})
test_that("Every chromosome must have markers",
	{
		#If there are no markers in a chromosome then throw an error
		expect_that({map <- new("map", vector(mode = "list", length = 10))}, throws_error())

		map[[2]] <- vector(mode="numeric", 0)
		expect_that(new("map", map), throws_error())
	})
test_that("Chromosomes must be numeric vectors",
	{
		#If a vector is not numeric then throw an error
		previousNames <- names(map[[1]])
		map[[1]] <- as.logical(map[[1]])
		names(map[[1]]) <- previousNames
		expect_that(new("map", map), throws_error())
	})
test_that("Chromosome names must be unique",
	{
		#If the chromosome names are not unique, throw an error
		names(map)[1:2] <- "repeated"
		expect_that(new("map", map), throws_error())
	})
test_that("Marker names must be unique",
	{
		#If the marker names are not unique, throw an error
		names(map[[1]]) <- names(map[[2]])
		expect_that(new("map", map), throws_error())
	})
