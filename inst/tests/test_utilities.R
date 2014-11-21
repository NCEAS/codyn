context("utilities")

test_that("utilities loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))

    library(codyn)

    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))

    #take a subset
    dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")

    # Does calComDat correctly transform a long dataframe into a matrix?
    time_matrix <- calComDat(data1 = dat1, "species", "year", "abundance")
    expect_true(is.data.frame(time_matrix))
    # rownames should be timevar
    expect_true(all(as.numeric(rownames(time_matrix)) %in% dat1[["year"]]))
    # colnames should be species
    expect_true(all(colnames(time_matrix) %in% dat1[["species"]]))

})
