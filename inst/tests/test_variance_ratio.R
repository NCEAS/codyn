context("variance ratio")


test_that("varianceratio function returns correct result", {
    # Ensure that tests 
    expect_that(length("a"), equals(1))
    
    # Load our example data set
    
    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
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
    datmat<-transpose_community(dat1, "species", "year", "abundance")
    
    #test the class returned with default settings
    myresults<-varianceratio(knz_001d, "subplot", "species", "year", "abundance", 1)
    expect_that(class(myresults), equals("data.frame"))
    expect_that(nrow(myresults), equals(1))
    
    #test that it also works with alternate column names
    knz_001d2<-knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    myresults2<-varianceratio(knz_001d2, "sub", "sp", "yr", "abund", 1)
    expect_that(sum(myresults2$VR), equals(sum(myresults$VR)))
    myresults2.2<-varianceratio(knz_001d2, "sub", "sp", "yr", "abund", 1,averagereps=FALSE)
    myresults2.3<-varianceratio(knz_001d2, NA, "sp", "yr", "abund", 1,averagereps=FALSE)
    myresults2.4<-varianceratio(knz_001d2, NA, "sp", "yr", "abund", 1,averagereps=TRUE)
    
    
    #test that it works even if there are additional unused columns
    knz_001d3<-knz_001d
    knz_001d3$site<-"KNZ"
    myresults3<-varianceratio(knz_001d3, "subplot", "species", "year", "abundance", 1)
    expect_that(myresults3$VR, is_identical_to(myresults$VR))
  
    
    ##test that works with replicate=NA
    myresults4<-varianceratio(dat1, replicate=NA, bootnumber=10)
    
    #test that works with one, specified replicate
    myresults5<-varianceratio(dat1,replicate="subplot", bootnumber=10)
    
    #test that works with replicate as a character
    myresults6<-varianceratio(dat3,replicate="subplot", bootnumber=10)
    
    #test that gives same for a single rep if average reps is false
    myresults7<-varianceratio(dat3,replicate="subplot", bootnumber=10,averagereps=FALSE)
    
    #test that all give the same VR value
    expect_that(myresults4$VR, is_identical_to(myresults5$VR))
    expect_that(myresults4$VR, is_identical_to(myresults6$VR))
    expect_that(myresults4$VR, is_identical_to(myresults7$VR))
    
    
    ##test calVR
    calVRresults<-calVR(datmat)
    expect_that(calVRresults, equals(myresults7$VR))
    expect_true(calVRresults>=0)
    expect_true(is.numeric(calVRresults))
    
  
    ##test calVR_longformdata
    lfresults<-calVR_longformdata(dat1, species="species", year="year", abundance="abundance")
    expect_that(lfresults, equals(myresults7$VR))
    expect_true(lfresults>=0)
    expect_true(is.numeric(lfresults))
    
    
    
    
    
    
})
