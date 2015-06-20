context("map validation")

test_that("map validation code functions", 
	{
		#If there are no chromosomes then throw an error
		expect_that({map <- new("map", list())}, throws_error())

		#If there are no markers in a chromosome then throw an error
		expect_that({map <- new("map", vector(mode = "list", length = 10))}, throws_error())
		
		map <- sim.map(len = rep(100, 3), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
		map[[2]] <- vector(mode="numeric", 0)
		expect_that(new("map", map), throws_error())

		#If a vector is not numeric then throw an error
		map <- sim.map(len = rep(100, 3), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
		previousNames <- names(map[[1]])
		map[[1]] <- as.logical(map[[1]])
		names(map[[1]]) <- previousNames
		expect_that(new("map", map), throws_error())

		#If the chromosome names are not unique, throw an error
		map <- sim.map(len = rep(100, 3), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
		names(map)[1:2] <- "repeated"
		expect_that(new("map", map), throws_error())

		#If the marker names are not unique, throw an error
		map <- sim.map(len = rep(100, 3), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
		names(map[[1]]) <- names(map[[2]])
		expect_that(new("map", map), throws_error())

		#An actual map passes validation
		map <- sim.map(len = rep(100, 3), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
		expect_identical({map <- new("map", map); validObject(map, complete=T)}, TRUE)
	})
