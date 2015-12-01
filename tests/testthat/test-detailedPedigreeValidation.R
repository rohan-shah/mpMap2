context("detailedPedigree validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)

test_that("Generated detailedPedigree passes validation", 
	{
		#An actual map passes validation
		expect_identical(validObject(pedigree, complete=TRUE), TRUE)
	})
test_that("detailedPedigree slot observed must have the correct length", 
	{
		#Longer observed slot throws an error
		copied <- pedigree
		copied@observed <- c(copied@observed, FALSE)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Shorter observed slot throws an error
		copied <- pedigree
		copied@observed <- head(copied@observed, -1)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slots observed and initial cannot contain NA values",
	{
		copied <- pedigree
		copied@observed[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- pedigree
		copied@initial[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})

test_that("Slot initial cannot have length 0",
	{
		copied <- pedigree
		copied@initial <- vector(mode="integer", length=0)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})

test_that("Slot initial must contain row indices",
	{
		#Cannot contain 0
		copied <- pedigree
		copied@initial[1] <- 0L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Cannot contain length + 1
		copied <- pedigree
		copied@initial[1] <- length(copied@observed)+1L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Initial lines cannot have mother or father",
	{
		copied <- pedigree
		copied@mother[copied@initial[2]] <- 1L
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- pedigree
		copied@father[copied@initial[2]] <- 1L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot initial must contain unique values",
	{
		copied <- pedigree
		copied@initial <- c(1L,1L,2L)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Initial lines must be at the start of the pedigree",
	{
		copied <- pedigree
		copied@mother <- c(0L, 0L, copied@mother)
		copied@father <- c(0L, 0L, copied@father)
		copied@lineNames <- c("aa", "bb", copied@lineNames)
		copied@initial <- c(3:4)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
