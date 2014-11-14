context("utilities")

test_that("utilities loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))

    library(codyn)

    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))

    # Does calComDat correctly transform a long dataframe into a matrix?
    dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")
    time_matrix <- codyn:::calComDat(data1 = dat1, "species", "year", "abundance")
    expect_true(is.data.frame(time_matrix))
#     knz_001d %>%
#       filter(subplot == "A_1") %>%
#       select(-subplot) %>%
#       codyn:::calComDat("species", "year", "abundance")
})
