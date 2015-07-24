context("combineGenotypes")

test_that("Function copes with input matrix having zero columns", 
	{
		#Three lines, no markers
		finals <- matrix(0L, nrow=3, ncol=0)
		hetData <- new("hetData")
		result <- mpMap2:::combineGenotypes(finals, hetData)

		#Returns three lines, no markers
		expect_equal(dim(result), c(3,0))
	})

test_that("Function copes with input matrix having zero rows", 
	{
		#no lines, 2 markers
		finals <- matrix(0L, nrow=0, ncol=4)
		hetData <- new("hetData")
		#Wrong number of hetData entries
		expect_that(mpMap2:::combineGenotypes(finals, hetData), throws_error())

		hetDataEntry <- rbind(c(0,0,0), c(0,1,2), c(1,0,2), c(1,1,1))
		hetData <- new("hetData", list(marker1 = hetDataEntry, marker2 = hetDataEntry))
		#Right number of hetData entries, returns 2 markers, no lines
		result <- mpMap2:::combineGenotypes(finals, hetData)
		expect_equal(dim(result), c(0, 2))
	})
test_that("Copes with final alleles not listed in hetData", 
	{
		finals <- matrix(0L, nrow=1, ncol=2)
		finals[] <- c(3,3)
		hetData <- new("hetData", list(marker1 = rbind(c(0,0,0), c(0,1,2), c(1,0,2), c(1,1,1))))
		expect_that(mpMap2:::combineGenotypes(finals, hetData), throws_error())
	})
test_that("Checking functionality", 
	{
		finals <- matrix(0L, nrow=1, ncol=2)
		hetData <- new("hetData", list(marker1 = rbind(c(0,0,0), c(0,1,2), c(1,0,2), c(1,1,1))))
		#Value of 0,0 -> 0
		finals[] <- c(0,0)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(0)))
		#Value of 1,0 -> 2
		finals[] <- c(1,0)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(2)))
		#Value of 0,1 -> 2
		finals[] <- c(0,1)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(2)))
		#Value of 1,1 -> 1
		finals[] <- c(1,1)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(1)))

		finals <- matrix(0L, nrow=1, ncol=2)
		hetData <- new("hetData", list(marker1 = rbind(c(10,10,10), c(10,20,30), c(20,10,30), c(20,20,20))))
		#Value of 10,10 -> 10
		finals[] <- c(10,10)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(10)))
		#Value of 20,10 -> 30
		finals[] <- c(20,10)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(30)))
		#Value of 10,20 -> 30
		finals[] <- c(10,20)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(30)))
		#Value of 20,20 -> 10
		finals[] <- c(20,20)
		expect_equal(mpMap2:::combineGenotypes(finals, hetData), rbind(c(20)))

		#Invalid values
		finals[] <- c(0,0)
		expect_that(mpMap2:::combineGenotypes(finals, hetData), throws_error())
		finals[] <- c(1,0)
		expect_that(mpMap2:::combineGenotypes(finals, hetData), throws_error())
		finals[] <- c(2,0)
		expect_that(mpMap2:::combineGenotypes(finals, hetData), throws_error())
	})