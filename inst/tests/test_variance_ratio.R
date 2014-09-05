context("calVR")

test_that("calVR function returns correct result", {
    # Ensure that tests 
    expect_that(length("a"), equals(1))
    
    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv("data/knz_001d.csv", sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))
})