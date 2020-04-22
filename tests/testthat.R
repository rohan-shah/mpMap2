library(testthat)
library(mpMap2)

try(omp_set_num_threads(1), silent = TRUE)
test_check("mpMap2")
