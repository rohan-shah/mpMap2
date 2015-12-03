context("rawSymmetricMatrix")
test_that("Check that validation works correctly", 
	{
		#Throws as data is too long
		expect_that(raw <- new("rawSymmetricMatrix", levels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), markers = c("a", "b", "c"), data = as.raw(as.integer(c(0:6)))), throws_error())
		#Throws as data is out of range
		expect_that(raw <- new("rawSymmetricMatrix", levels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), markers = c("a", "b", "c"), data = as.raw(as.integer(c(0:4,6)))), throws_error())
		#Throws as too many levels
		expect_that(raw <- new("rawSymmetricMatrix", levels = c(0:300)/600, markers = c("a", "b", "c"), data = as.raw(as.integer(c(0:4,6)))), throws_error())
		#Throws as values in slot levels are not ordered
		expect_that(raw <- new("rawSymmetricMatrix", levels = c(0:100, 99)/300, markers = c("a", "b", "c"), data = as.raw(as.integer(c(0:4,6)))), throws_error())
		expect_that(raw <- new("rawSymmetricMatrix", levels = c(0:90, 100, 99)/300, markers = c("a", "b", "c"), data = as.raw(as.integer(c(0:4,6)))), throws_error())
	})
test_that("Checking that indexing works correctly",
	{
		raw  <- new("rawSymmetricMatrix", levels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), markers = c("a", "b", "c"), data = as.raw(as.integer(c(0:5))))
		expect_that(raw[1:4, 1:3], throws_error())
		expect_that(raw[1:4, 1:4], throws_error())
		expect_that(raw[4, 1], throws_error())
		expect_that(raw[1, 4], throws_error())
		expect_identical(raw[1,1], 0)
		expect_identical(raw[1,2], 0.1)
		expect_identical(raw[2,1], 0.1)
		expect_identical(raw[3,1], 0.3)
		expect_identical(raw[1,3], 0.3)
		expect_identical(raw[1,3,drop=F], cbind(0.3))
		expect_identical(raw[3,1,drop=F], cbind(0.3))

		expect_identical(raw[1,1:3], c(0, 0.1, 0.3))
		expect_identical(raw[1:3,1], c(0, 0.1, 0.3))
		
		expect_identical(raw[1,1:3,drop=F], rbind(c(0, 0.1, 0.3)))
		expect_identical(raw[1:3,1,drop=F], cbind(c(0, 0.1, 0.3)))

		#Check NA values work
		raw  <- new("rawSymmetricMatrix", levels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), markers = c("a", "b", "c"), data = as.raw(as.integer(c(0:4, 255))))
		expect_identical(raw[3, 3], as.numeric(NA))
	})
test_that("Checking that source and destination cannot be the same",
	{
		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(rep(0, 10)))
		expect_that(.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 1:4, 1:4, raw@data), throws_error())

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(rep(0, 10)))
		expect_that(.Call("assignRawSymmetricMatrixDiagonal", raw, 1:4, 1:4, raw@data), throws_error())
	})
test_that("Checking that subsetting works correctly",
	{
		raw  <- new("rawSymmetricMatrix", levels = c(0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		expect_equal(subset(raw, markers = 1:3)@data, as.raw(as.integer(c(0,1,2,3,4,5))))
		expect_equal(subset(raw, markers = 2:4)@data, as.raw(as.integer(c(2,4,5,7,8,9))))
		expect_equal(subset(raw, markers = c(1,3,4))@data, as.raw(as.integer(c(0,3,5,6,8,9))))
	})
test_that("Checking that rawSymmetricMatrixEstimateRF works correctly, for ordered indices",
	{
		#First with a constant value of 100
		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 1, 1, as.raw(100))
		expect_equal(raw@data[1], as.raw(100))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 1:2, 1:2, as.raw(rep(100, 3)))
		expect_equal(sum(raw@data == as.raw(100)), 3)
		expect_equal(raw@data[1:3], as.raw(c(100,100,100)))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 2:3, 2:3, as.raw(rep(100, 3)))
		expect_equal(sum(raw@data == as.raw(100)), 3)
		expect_equal(raw@data[c(3,5,6)], as.raw(c(100,100,100)))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 2, 2, as.raw(100))
		expect_equal(raw@data[3], as.raw(100))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 1:3, 1:3, as.raw(rep(100, 6)))
		expect_equal(sum(raw@data == as.raw(100)), 6)
		expect_equal(raw@data[1:6], as.raw(rep(100, 6)))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 1:2, 3:4, as.raw(rep(100, 4)))
		expect_equal(sum(raw@data == as.raw(100)), 4)
		expect_equal(raw@data[c(4,5,7,8)], as.raw(rep(100, 4)))

		#Now with different values
		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 1:2, 3:4, as.raw(c(100, 101, 102, 103)))
		expect_equal(sum(raw@data > 99), 4)
		expect_equal(raw@data[c(4,5,7,8)], as.raw(c(100, 101, 102, 103)))
	})
test_that("Checking that rawSymmetricMatrixEstimateRF works correctly, for out-of-order indices",
	{
		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 2:1, 1:2, as.raw(rep(100, 3)))
		expect_equal(sum(raw@data == as.raw(100)), 3)
		expect_equal(raw@data[1:3], as.raw(c(100,100,100)))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 3:2, 3:2, as.raw(rep(100, 3)))
		expect_equal(sum(raw@data == as.raw(100)), 3)
		expect_equal(raw@data[c(3,5,6)], as.raw(c(100,100,100)))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 3:1, 3:1, as.raw(rep(100, 6)))
		expect_equal(sum(raw@data == as.raw(100)), 6)
		expect_equal(raw@data[1:6], as.raw(rep(100, 6)))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 2:1, 4:3, as.raw(rep(100, 4)))
		expect_equal(sum(raw@data == as.raw(100)), 4)
		expect_equal(raw@data[c(4,5,7,8)], as.raw(rep(100, 4)))

		#Now with different values
		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(as.integer(c(0:9))))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 2:1, 4:3, as.raw(c(100, 101, 102, 103)))
		expect_equal(sum(raw@data > 99), 4)
		expect_equal(raw@data[c(4,5,7,8)], as.raw(c(103, 102, 101, 100)))

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(rep(0, 10)))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, 2:1, 4:3, as.raw(c(1,2,3,4)))
		expect_identical(as.vector(raw[1:2,3:4]), ((0:9)/18)[c(4,3,2,1)+1])

		raw  <- new("rawSymmetricMatrix", levels = (0:9)/18, markers = c("a", "b", "c", "d"), data = as.raw(rep(0, 10)))
		.Call("assignRawSymmetricMatrixFromEstimateRF", raw, c(4,1), c(1,4), as.raw(c(0, 0, 9)))
		expect_identical(raw[1:4,1:4], cbind(c(0,0,0,0.5), rep(0, 4), rep(0, 4), c(0.5,0,0,0)))
	})
