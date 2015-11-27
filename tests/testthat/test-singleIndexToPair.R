context("singleIndexToPair")
singleIndexToPair <- function(marker1Start, marker1End, marker2Start, marker2End, index)
{
	.Call("singleIndexToPair", marker1Start, marker1End, marker2Start, marker2End, index, PACKAGE="mpMap2")
}
parameteriseRegion <- function(marker1Start, marker1End, marker2Start, marker2End)
{
	nValues <- .Call("countValuesToEstimate", marker1Start, marker1End, marker2Start, marker2End, PACKAGE="mpMap2")
	return(do.call(rbind, sapply(1:nValues, function(x) singleIndexToPair(marker1Start, marker1End, marker2Start, marker2End, x), simplify=FALSE)))
}

test_that("Checking that singleIndexToPair works for rectangular regions",
	{
		expect_equal(parameteriseRegion(1, 2, 1, 11), cbind(1:10, rep(1, 10)))
		expect_equal(parameteriseRegion(1, 3, 2, 11), rbind(cbind(2:10, rep(1, 9)), cbind(2:10, rep(2, 9))))

		expect_equal(parameteriseRegion(2, 3, 2, 11), cbind(2:10, rep(2, 9)))
		expect_equal(parameteriseRegion(2, 4, 3, 11), rbind(cbind(3:10, rep(2, 8)), cbind(3:10, rep(3, 8))))
	})
