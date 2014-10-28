context("rankshifts")

test_that("rankshifts loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))
    
    library(codyn)
    
    # Load our example data set
    knz_001d <- read.csv(system.file("data", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))
    
    # Basic test if mean rank produces data frame with right structure and values
    result <- meanrank(knz_001d)
    expect_that(length(names(result)), equals(2))
    expect_that(names(result)[[1]], matches("year_pair"))
    expect_that(class(result[,1]), matches("factor"))
    expect_that(names(result)[[2]], matches("MRS"))
    expect_that(class(result[,2]), matches("numeric"))
    expect_that(nrow(result), equals(23))
    expect_that(result[[1,2]], equals(332.9921, tolerance=.00001))
    expect_that(result[[23,2]], equals(670.6388, tolerance=.00001))
})