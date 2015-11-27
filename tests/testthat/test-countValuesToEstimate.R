context("countValuesToEstimate")
countValuesToEstimate <- function(marker1Start, marker1End, marker2Start, marker2End)
{
	.Call("countValuesToEstimate", marker1Start, marker1End, marker2Start, marker2End, PACKAGE="mpMap2")
}

test_that("Checking that countValuesToEstimate works for rectangular regions",
	{
		expect_equal(countValuesToEstimate(1,2,1,11), 10)
		expect_equal(countValuesToEstimate(1,3,2,11), 18)
		expect_equal(countValuesToEstimate(1,4,3,11), 24)
	})
test_that("Checking that countValuesToEstimate works for triangular regions",
	{
		expect_equal(countValuesToEstimate(1,11,1,11), 55)
		expect_equal(countValuesToEstimate(2,11,2,11), 45)
		expect_equal(countValuesToEstimate(1,10,1,3), 3)
	})
test_that("Checking that countValuesToEstimate works for regions that start with a rectangular region",
	{
		expect_equal(countValuesToEstimate(1,11,9,11), 19)
		expect_equal(countValuesToEstimate(1,11,8,11), 27)
		expect_equal(countValuesToEstimate(1,4,2,5), 8)
	})
test_that("Checking that countValuesToEstimate works for regions that have a rectangular region on the right",
	{
		expect_equal(countValuesToEstimate(1,3,1,11), 19)
		expect_equal(countValuesToEstimate(1,4,1,11), 27)
	})
