context("synchrony")

test_that("synchrony loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))
    
    library(codyn)
    
    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv(system.file("data", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))
    
    # TODO: Test the synchrony function
    
    # TODO: Test the synch_onerep function
    
})