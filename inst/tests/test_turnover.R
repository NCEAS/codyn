context("turnover")

test_that("turnover loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))
    
    library(codyn)
    
    # Load our example data set
    # data("knz_001d", package="codyn")  # This doesn't work for CSV files :(
    knz_001d <- read.csv(system.file("extdata", "knz_001d.csv", package="codyn"), sep=",", header=TRUE)
    expect_that(names(knz_001d)[4], equals("abundance"))
    

    #give new column names
    knz_001d2 <- knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    
    #add a random character and factor column
    knz_001d2$randcharacter<-"rchar"
    knz_001d2$randfactor<-as.factor(knz_001d2$randcharacter)
    
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
    
    #test that works regardless of whether parameter is specified or just ordered
    myresults<-turnover(knz_001d, replicate="subplot", species="species", year="year", abundance="abundance", metric="turnover")
    myresults2<-turnover(knz_001d, "subplot", "species", "year", "abundance", metric="turnover")    
    expect_that(myresults, is_identical_to(myresults2))  
    
    #test that turnover is the default metric
    myresults3<-turnover(knz_001d, replicate="subplot", species="species", year="year", abundance="abundance")
    expect_that(myresults, is_identical_to(myresults3))  
    
    #test that works with different column nmaes
    myresults4<-turnover(knz_001d2, replicate="sub", species="sp", year="yr", abundance="abund", metric="turnover")
    expect_that(sum(myresults$turnover), equals(sum(myresults4$turnover)))
    
    #test that works with different column orders if names specified
    myresults5<-turnover(knz_001d2, abundance="abund",replicate="sub", species="sp", year="yr",  metric="turnover")
    expect_that(sum(myresults$turnover), equals(sum(myresults5$turnover)))
    
    #test that gives an error if abundance is a character or factor
    myresults4<-turnover(knz_001d2, replicate="sub", species="sp", year="yr", abundance="randcharacter", metric="turnover")
    
    
    ##TODO: Decide how to handle factor and character abundance columns
    #Test that returns a warning when abundance is a character or factor column
    #expect_error(turnover(knz_001d2, replicate="sub", species="sp", year="yr", abundance="randcharacter", metric="turnover"))
    #expect_error(turnover(knz_001d2, replicate="sub", species="sp", year="yr", abundance="randfactor", metric="turnover"))
    
    


    
})

