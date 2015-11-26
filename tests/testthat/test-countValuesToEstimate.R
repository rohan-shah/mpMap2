context("countValuesToEstimate")
countValuesToEstimate <- function(marker1Start, marker1End, marker2Start, marker2End)
{
	.Call("countValuesToEstimate", marker1Start, marker1End, marker2Start, marker2End, PACKAGE="mpMap2")
}

test_that("Checking that countValuesToEstimate works for rectangular regions",
	{
		test_equal(countValuesToEstimate(1,2,1,11), 10)
		test_equal(countValuesToEstimate(1,3,2,11), 18)
		test_equal(countValuesToEstimate(1,4,3,11), 24)
	})
test_that("Checking that countValuesToEstimate works for triangular regions",
	{
		test_equal(countValuesToEstimate(1,11,1,11), 55)
		test_equal(countValuesToEstimate(2,11,2,11), 45)
	})
test_that("Checking that countValuesToEstimate works for regions that start with a rectangular region",
	{
		test_equal(countValuesToEstimate(1,11,9,11), 19)
		test_equal(countValuesToEstimate(1,11,8,11), 27)
	})
test_that("Checking that countValuesToEstimate works for regions that have a rectangular region on the right",
	{
		test_equal(countValuesToEstimate(1,3,1,11), 19)
		test_equal(countValuesToEstimate(1,4,1,11), 27)
	})
