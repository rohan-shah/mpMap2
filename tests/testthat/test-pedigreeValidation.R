context("pedigree validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree <- as(pedigree, "pedigree")

test_that("two-parent pedigree passes validation",
	{
		expect_identical(validObject(pedigree, complete=TRUE), TRUE)
	})
test_that("pedigree slots must have correct length", 
	{
		#Shorter lineNames gives an error
		copied <- pedigree
		copied@lineNames <- head(copied@lineNames, -1L)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Longer lineNames gives an error
		copied <- pedigree
		copied@lineNames <- c(copied@lineNames, "aaaaaaa")
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Shorter mother gives an error
		copied <- pedigree
		copied@mother <- head(copied@mother, -1L)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Longer mother gives an error
		copied <- pedigree
		copied@mother <- c(copied@mother, 1L)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Shorter father given an error
		copied <- pedigree
		copied@father <- head(copied@father, -1L)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Longer father gives an error
		copied <- pedigree
		copied@father <- c(copied@father, 1L)
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
test_that("Father and mother must come before offspring in pedigree",
	{
		copied <- pedigree
		copied@mother[10] <- 11L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- pedigree
		copied@father[10] <- 11L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Attribute selfing must be either \"infinite\" or \"finite\"",
	{
		copied <- pedigree
		copied@selfing <- "other"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("lineNames must be unique",
	{	
		copied <- pedigree
		copied@lineNames[3] <- copied@lineNames[4]
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Number of founders must be at least 2",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 1000, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		pedigree@selfing <- "finite"
		pedigree@initial <- 1L
		pedigree@mother <- pedigree@father <- 0:1005
		expect_that(validObject(pedigree, complete=TRUE), throws_error())
	})
test_that("Replacement function lineNames checks for duplicates and correct length",
	{
		copied <- pedigree
		expect_that(lineNames(copied)[3] <- lineNames(copied)[4], throws_error("duplicates"))
		expect_that(lineNames(copied) <- lineNames(copied)[-1], throws_error("wrong length"))
	})
rm(pedigree)
