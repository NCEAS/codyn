context("variance ratio")

# Utility function to run calVR2 on multiple segments of a data frame
calVRs<-function(data1, rep, species, year, abundance){
    X <- split(data1, data1[rep])
    out<-lapply(X, FUN=calVR2, species, year, abundance)
    output<-cbind((names(out)), as.data.frame(unlist(out)))
    names(output)<-c(rep, "VR")
    row.names(output)<-NULL
    return(output)
}

test_that("calVR function returns correct result", {
    # Ensure that tests 
    expect_that(length("a"), equals(1))
    
    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv("../../data/knz_001d.csv", sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))
    myresults<-calVRs(knz_001d, "subplot", "species", "year", "abundance")
    expect_that(class(myresults), equals("data.frame"))
    expect_that(nrow(myresults), equals(20))
    
    #test that it also works with alternate column names
    knz_001d2<-knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    myresults2<-calVRs(knz_001d2, "sub", "sp", "yr", "abund")
    expect_that(sum(myresults2$VR), equals(sum(myresults$VR)))
    
    #test that it works even if there are additional unused columns
    knz_001d3<-knz_001d
    knz_001d3$site<-"KNZ"
    myresults3<-calVRs(knz_001d3, "subplot", "species", "year", "abundance")
    expect_that(myresults3, is_identical_to(myresults))
  
    })