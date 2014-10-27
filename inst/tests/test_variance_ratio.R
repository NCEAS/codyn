context("variance ratio")


test_that("VR function returns correct result", {
    # Ensure that tests 
    expect_that(length("a"), equals(1))
    
    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv("../../data/knz_001d.csv", sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))
    myresults<-VR(knz_001d, "subplot", "species", "year", "abundance", 1)
    expect_that(class(myresults), equals("data.frame"))
    expect_that(nrow(myresults), equals(1))
    
    #test that it also works with alternate column names
    knz_001d2<-knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    myresults2<-VR(knz_001d2, "sub", "sp", "yr", "abund", 1)
    expect_that(sum(myresults2$VR), equals(sum(myresults$VR)))
    
    #test that it works even if there are additional unused columns
    knz_001d3<-knz_001d
    knz_001d3$site<-"KNZ"
    myresults3<-VR(knz_001d3, "subplot", "species", "year", "abundance", 1)
    expect_that(myresults3$VR, is_identical_to(myresults$VR))
  
    })