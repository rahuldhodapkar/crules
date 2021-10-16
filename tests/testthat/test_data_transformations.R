context("Test data transformations")
library(crules)
library(Matrix)

test_that("Test 'GenerateIncidenceMatrix'", {

    x <- sparseMatrix(i=c(), j=c(), dims=c(10,10))
    x <- as(x, 'dgCMatrix')
    x[1,1] <- 5
    x[2,1] <- 1
    x[3,1] <- 4
    y <- GenerateIncidenceMatrix(x)

    expect_equal(y[1,1], 1)
    expect_equal(y[2,1], 1)

    y <- GenerateIncidenceMatrix(x, percentile.cutoff=0.95)
    expect_equal(y[2,1], 0)
})