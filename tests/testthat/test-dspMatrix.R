context("dspMatrix")
test_that("Checking that C assignment call works correctly",
	{
		#First with a constant value of 100
		dsp <- new("dspMatrix", x = as.numeric(1:10), Dim = c(4L,4L))
		.Call("assignDspMatrixFromEstimateRF", dsp, 1, 1, 100)
		expect_equal(dsp@x[1], 100)

		dsp <- new("dspMatrix", x = as.numeric(1:10), Dim = c(4L,4L))
		.Call("assignDspMatrixFromEstimateRF", dsp, 1:2, 1:2, rep(100, 3))
		expect_equal(sum(dsp@x == 100), 3)
		expect_equal(dsp@x[1:3], c(100,100,100))

		dsp <- new("dspMatrix", x = as.numeric(1:10), Dim = c(4L,4L))
		.Call("assignDspMatrixFromEstimateRF", dsp, 2:3, 2:3, rep(100, 3))
		expect_equal(sum(dsp@x == 100), 3)
		expect_equal(dsp@x[c(3,5,6)], c(100,100,100))

		dsp <- new("dspMatrix", x = as.numeric(1:10), Dim = c(4L,4L))
		.Call("assignDspMatrixFromEstimateRF", dsp, 2, 2, 100)
		expect_equal(dsp@x[3], 100)

		dsp <- new("dspMatrix", x = as.numeric(1:10), Dim = c(4L,4L))
		.Call("assignDspMatrixFromEstimateRF", dsp, 1:3, 1:3, rep(100, 6))
		expect_equal(sum(dsp@x == 100), 6)
		expect_equal(dsp@x[1:6], rep(100, 6))

		dsp <- new("dspMatrix", x = as.numeric(1:10), Dim = c(4L,4L))
		.Call("assignDspMatrixFromEstimateRF", dsp, 1:2, 3:4, rep(100, 4))
		expect_equal(sum(dsp@x == 100), 4)
		expect_equal(dsp@x[c(4,5,7,8)], rep(100, 4))

		#Now with different values
		dsp <- new("dspMatrix", x = as.numeric(1:10), Dim = c(4L,4L))
		.Call("assignDspMatrixFromEstimateRF", dsp, 1:2, 3:4, c(100, 101, 102, 103))
		expect_equal(sum(dsp@x > 99), 4)
		expect_equal(dsp@x[c(4,5,7,8)], c(100, 101, 102, 103))
	})
