context("singleIndexToPair")
singleIndexToPair <- function(markerRows, markerColumns, index)
{
	.Call("singleIndexToPair", markerRows, markerColumns, index, PACKAGE="mpMap2")
}
parameteriseRegion <- function(markerRows, markerColumns)
{
	nValues <- .Call("countValuesToEstimate", markerRows, markerColumns, PACKAGE="mpMap2")
	return(do.call(rbind, sapply(1:nValues, function(x) singleIndexToPair(markerRows, markerColumns, x), simplify=FALSE)))
}

test_that("Checking that singleIndexToPair works for rectangular regions",
	{
		expect_equal(parameteriseRegion(1, 1:10), cbind(rep(1, 10), 1:10))
		expect_equal(parameteriseRegion(1:2, 2:10), cbind(rep(1:2, times = 9), rep(2:10, each = 2)))

		expect_equal(parameteriseRegion(2, 2:10), cbind(rep(2, 9), 2:10))
		expect_equal(parameteriseRegion(2:3, 3:10), cbind(rep(2:3, times = 8), rep(3:10, each = 2)))
	})
test_that("Checking that singleIndexToPair works for triangular regions",
	{
		expect_equal(parameteriseRegion(1, 1), rbind(c(1,1)))
		expect_equal(parameteriseRegion(1:2, 1:2), rbind(c(1,1), c(1,2), c(2,2)))
		expect_equal(parameteriseRegion(1:3, 1:3), rbind(c(1,1), c(1,2), c(2,2), c(1,3), c(2, 3), c(3, 3)))

		expect_equal(parameteriseRegion(2:3, 2:3), rbind(c(2,2), c(2,3), c(3,3)))
		expect_equal(parameteriseRegion(2:4, 2:4), rbind(c(2,2), c(2,3), c(3,3), c(2,4), c(3,4), c(4,4)))

		expect_equal(parameteriseRegion(1:3, 1:2), rbind(c(1,1), c(1,2), c(2,2)))
		expect_equal(parameteriseRegion(1:9, 1:2), rbind(c(1,1), c(1,2), c(2,2)))
	})
test_that("Checking that singleIndexToPair works for regions that have a rectangular region on the right",
	{
		expect_equal(parameteriseRegion(1:2, 1:3), rbind(c(1,1), c(1,2), c(2,2), c(1,3), c(2,3)))
		expect_equal(parameteriseRegion(2:3, 2:4), rbind(c(2,2), c(2,3), c(3,3), c(2,4), c(3,4)))
		expect_equal(parameteriseRegion(1:2, 1:4), rbind(c(1,1), c(1,2), c(2,2), c(1,3), c(2,3), c(1,4), c(2,4)))
		expect_equal(parameteriseRegion(2:3, 1:4), rbind(c(2,2), c(2,3), c(3,3), c(2,4), c(3,4)))
		expect_equal(parameteriseRegion(2:3, 2:5), rbind(c(2,2), c(2,3), c(3,3), c(2,4), c(3,4), c(2,5), c(3,5)))

	})
test_that("Checking that singleIndexToPair works for regions that start with a rectangular region",
	{
		expect_equal(parameteriseRegion(1:3, 2:3), rbind(c(1,2), c(2,2), c(1,3), c(2,3), c(3,3)))
		expect_equal(parameteriseRegion(1:3, 2:4), rbind(c(1,2), c(2,2), c(1,3), c(2,3), c(3,3), c(1,4), c(2,4), c(3,4)))
		expect_equal(parameteriseRegion(1:3, 2:5), rbind(c(1,2), c(2,2), c(1,3), c(2,3), c(3,3), c(1,4), c(2,4), c(3,4), c(1,5), c(2,5), c(3,5)))
	})
rm(singleIndexToPair, parameteriseRegion)
