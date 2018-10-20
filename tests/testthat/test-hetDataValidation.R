context("hetData validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
map <- qtl::sim.map(len = rep(100, 1), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
hetData <- cross@geneticData[[1]]@hetData

test_that("Simulated hetData passes validation",
	{
		expect_identical(validObject(hetData, complete=TRUE), TRUE)
	})

test_that("Every entry must have a valid name",
	{
		copied <- hetData
		names(copied) <- NULL
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- hetData
		names(copied)[1] <- ""
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- hetData
		names(copied)[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Names must be unique",
	{
		copied <- hetData
		names(copied)[1] <- names(copied[2])
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
#Now the stuff checked by the C code
test_that("hetData entries must be numeric matrices",
	{
		#Numeric is OK
		copied <- hetData
		storage.mode(copied[[1]]) <- "numeric"
		expect_identical(validObject(copied, complete=TRUE), TRUE)

		#Integer is OK
		copied <- hetData
		storage.mode(copied[[1]]) <- "integer"
		expect_identical(validObject(copied, complete=TRUE), TRUE)

		#Character isn't
		copied <- hetData
		storage.mode(copied[[1]]) <- "logical"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("hetData entries must have dimension attribute",
	{
		copied <- hetData
		dim(copied[[1]]) <- NULL
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("hetData entries must be two dimensional",
	{
		copied <- hetData
		copied[[1]] <- c(copied[[1]], copied[[1]])
		dim(copied[[1]]) <- c(dim(hetData[[1]]), 2)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("hetData entries must have three columns",
	{
		copied <- hetData
		copied[[1]] <- cbind(copied[[1]], 1:nrow(copied[[1]]))
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("hetData encodings must be symmetric",
	{
		copied <- hetData
		index1 <- which(copied[[1]][,1] == 1 & copied[[1]][,2] == 2)
		index2 <- which(copied[[1]][,1] == 2 & copied[[1]][,2] == 1)
		copied[[1]][index2, 3] <- copied[[1]][index1, 3] + 1
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Duplicate rows not allowed in hetData",
	{
		copied <- hetData
		copied[[1]] <- rbind(copied[[1]], copied[[1]][1,])
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Contradictory rows not allowed in hetData",
	{
		copied <- hetData
		#Duplicate the first row, but give it a different encoding for the finals
		extraRow <- copied[[1]]
		extraRow[3] <- extraRow[3] + 1
		copied[[1]] <- rbind(copied[[1]], extraRow)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
rm(pedigree, map, cross, hetData)
