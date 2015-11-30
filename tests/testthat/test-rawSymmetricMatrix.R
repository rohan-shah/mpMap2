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
