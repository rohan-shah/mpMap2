context("pedigree validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree <- as(pedigree, "pedigree")

test_that("two-parent pedigree passes validation",
	{
		expect_identical(validObject(pedigree, complete=TRUE), TRUE)
	})
test_that("pedigree slots must have correct length", 
	{
		#Wrong length lineNames gives an error
		copied <- pedigree
		copied@lineNames <- head(copied@lineNames, -1)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Wrong length mother gives an error
		copied <- pedigree
		copied@mother <- head(copied@mother, -1)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Wrong length father given an error
		copied <- pedigree
		copied@father <- head(copied@father, -1)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("pedigree slots cannot contain NA", 
	{
		pedigreeLength <- length(pedigree@father)

		#NA in lineNames gives an error
		copied <- pedigree
		copied@lineNames[pedigreeLength] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#NA in mother gives an error
		copied <- pedigree
		copied@mother[pedigreeLength] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#NA in father gives an error
		copied <- pedigree
		copied@lineNames[pedigreeLength] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Pedigree slots mother and father must be 0 or row indices",
	{
		copied <- pedigree
		copied@mother[1] <- -1L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- pedigree
		copied@mother[1] <- length(copied@mother)+1L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- pedigree
		copied@father[1] <- -1L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- pedigree
		copied@father[1] <- length(copied@mother)+1L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("lineNames must be unique",
	{	
		pedigree@lineNames[3] <- pedigree@lineNames[4]
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Father and mother must come before offspring in pedigree",
	{
		copied <- pedigree
		copied@mother[10] <- 11L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- pedigree
		copied@father[10] <- 11L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})