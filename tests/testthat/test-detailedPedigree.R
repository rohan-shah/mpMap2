context("detailedPedigree")

test_that("Test use of mother / father names as input to detailedPedigree", 
	{
		lineNames <- c("founder1", "founder2", "child1", "child2")
		pedigree <- detailedPedigree(lineNames = c("founder1", "founder2", "child1", "child2"), mother = c("", "", "founder1", "child1"), father = c("", "", "founder2", "child1"), initial = c("founder1", "founder2"), observed = "child2", selfing = "finite")
		expect_identical(pedigree@initial, 1:2)
		expect_identical(pedigree@observed, c(FALSE, FALSE, FALSE, TRUE))
		expect_identical(pedigree@lineNames, lineNames)
		expect_identical(pedigree@mother, c(0L, 0L, 1L, 3L))
		expect_identical(pedigree@father, c(0L, 0L, 2L, 3L))
	})
test_that("Test use of mother / father indices as input to detailedPedigree", 
	{
		lineNames <- c("founder1", "founder2", "child1", "child2")
		pedigree <- detailedPedigree(lineNames = c("founder1", "founder2", "child1", "child2"), mother = c(0, 0, 1, 3), father = c(0, 0, 2, 3), initial = 1:2, observed = 4, selfing = "finite")
		expect_identical(pedigree@initial, 1:2)
		expect_identical(pedigree@observed, c(FALSE, FALSE, FALSE, TRUE))
		expect_identical(pedigree@lineNames, lineNames)
		expect_identical(pedigree@mother, c(0L, 0L, 1L, 3L))
		expect_identical(pedigree@father, c(0L, 0L, 2L, 3L))
	})
