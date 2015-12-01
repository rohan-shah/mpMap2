context("countValuesToEstimate")
countValuesToEstimate <- function(markerRows, markerColumns)
{
	.Call("countValuesToEstimate", markerRows, markerColumns, PACKAGE="mpMap2")
}

test_that("Checking that countValuesToEstimate works for rectangular regions",
	{
		expect_equal(countValuesToEstimate(1,1:10), 10)
		expect_equal(countValuesToEstimate(1:2,2:10), 18)
		expect_equal(countValuesToEstimate(1:3,3:10), 24)
	})
test_that("Checking that countValuesToEstimate works for triangular regions",
	{
		expect_equal(countValuesToEstimate(1:10,1:10), 55)
		expect_equal(countValuesToEstimate(2:10,2:10), 45)
		expect_equal(countValuesToEstimate(1:9,1:2), 3)
	})
test_that("Checking that countValuesToEstimate works for regions that start with a rectangular region",
	{
		expect_equal(countValuesToEstimate(1:10,9:10), 19)
		expect_equal(countValuesToEstimate(1:10,8:10), 27)
		expect_equal(countValuesToEstimate(1:3,2:4), 8)
	})
test_that("Checking that countValuesToEstimate works for regions that have a rectangular region on the right",
	{
		expect_equal(countValuesToEstimate(1:2,1:10), 19)
		expect_equal(countValuesToEstimate(1:3,1:10), 27)
	})
rm(countValuesToEstimate)
