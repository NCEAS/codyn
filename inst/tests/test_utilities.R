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
    com_wide_df <- calComDat(data1 = dat1, "species", "year", "abundance")
    expect_true(is.data.frame(com_wide_df))
    # rownames should be timevar
    expect_true(all(as.numeric(rownames(com_wide_df)) %in% dat1[["year"]]))
    # colnames should be species
    expect_true(all(colnames(com_wide_df) %in% dat1[["species"]]))
})

test_that("Name checking works", {
  knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"),
                       sep=",", header=TRUE)
  expect_error(check_names(given = c("AAA", "year", "subplot", "abundance"),
                           data = knz_001d),
                 "data does not have name .*")
  expect_error(check_names(given = c("species", "BBB", "subplot", "abundance"),
                           data = knz_001d),
               "data does not have name .*")
  expect_error(check_names(given = c("species", "year", "CCC", "abundance"),
                           data = knz_001d),
               "data does not have name .*")
  expect_error(check_names(given = c("species", "year", "subplot", "DDD"),
                           data = knz_001d),
               "data does not have name .*")
})
