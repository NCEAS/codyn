context("variance ratio")

test_that("varianceratio function returns correct result", {
    # Ensure that tests 
    expect_that(length("a"), equals(1))
    
    # Load our example data set
    data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    expect_that(names(knz_001d)[4], equals("abundance"))
    
    
    #give new column names
    knz_001d2 <- knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    
    #take a subset
    dat1 <- subset(knz_001d, knz_001d$subplot=="A_1")
    
    #rename the subset
    dat2<-dat1
    names(dat2)=c("sp", "yr", "sub", "abund")
    
    #make subplot a character
    dat3<-dat1
    dat3$subplot<-as.character(dat3$subplot)
    
    #take a two replicate subset
    dat4<-subset(knz_001d, subplot=="A_1"|subplot=="A_2")
    
    #make a species matrix
    datmat<-transpose_community(dat1, "year",  "species", "abundance")
    
    #test the class returned with default settings
    myresults<-varianceratio(knz_001d, time.var="year", species.var="species", 
                             abundance.var="abundance",  bootnumber=1, replicate="subplot")
    expect_equal(class(myresults), "data.frame")
    expect_equal(nrow(myresults), 1)
    expect_equal(myresults$VR, 1.01443, tolerance=0.00001)
    
    #test that it also works with alternate column names
    knz_001d2<-knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    myresults2<-varianceratio(knz_001d2,"yr", "sp", "abund", 1, "sub")
    expect_that(sum(myresults2$VR), equals(sum(myresults$VR)))
    myresults2.2<-varianceratio(knz_001d2, "yr", "sp",  "abund", 1, "sub", average.replicates=FALSE)
    expect_warning(varianceratio(knz_001d2, "yr", "sp",  "abund", 1, NA, average.replicates=FALSE))
    expect_warning(varianceratio(knz_001d2,"yr", "sp",  "abund", 1, NA, average.replicates=TRUE))
    
    #test that it works even if there are additional unused columns
    knz_001d3<-knz_001d
    knz_001d3$site<-"KNZ"
    myresults3<-varianceratio(knz_001d3, "year", "species",  "abundance",1, "subplot")
    expect_that(myresults3$VR, is_identical_to(myresults$VR))
  
    ##test that works with replicate=NA
    myresults4<-varianceratio(dat1, replicate=NA, bootnumber=10)
    
    #test that works with one, specified replicate
    myresults5<-varianceratio(dat1,replicate="subplot", bootnumber=10)
    
    #test that works with replicate as a character
    myresults6<-varianceratio(dat3,replicate="subplot", bootnumber=10)
    
    #test that gives same for a single rep if average reps is false
    myresults7<-varianceratio(dat3, bootnumber=1, replicate="subplot", average.replicates=FALSE)
    
    #test that all give the same VR value
    expect_that(myresults4$VR, is_identical_to(myresults5$VR))
    expect_that(myresults4$VR, is_identical_to(myresults6$VR))
    expect_that(myresults4$VR, is_identical_to(myresults7$VR))
    
    ##test varianceratio_matrixdata
    calVRresults<-varianceratio_matrixdata(datmat)
    expect_that(calVRresults, equals(myresults7$VR))
    expect_true(calVRresults>=0)
    expect_true(is.numeric(calVRresults))
    
    ##test varianceratio_longformdata
    lfresults<-varianceratio_longformdata(dat1,  time.var="year",species.var="species", abundance.var="abundance")
    expect_that(lfresults, equals(myresults7$VR))
    expect_true(lfresults>=0)
    expect_true(is.numeric(lfresults))
    
    ## test that bad data values are handled with graceful error messages
    bad_data <- knz_001d
    bad_data['abundance'][1,] <- 'dung'
    expect_error(varianceratio(bad_data, "year", "species", "abundance",  bootnumber=1, replicate="subplot"), "is.numeric")
})
