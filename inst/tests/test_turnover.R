context("turnover")

test_that("turnover loads and returns correct result", {
    # Ensure that trivial tests work correctly
    expect_that(length("a"), equals(1))
    
    library(codyn)
    
    # Load our example data set
    data(knz_001d)
    expect_that(names(knz_001d)[4], equals("abundance"))
    
    #give new column names
    knz_001d2 <- knz_001d
    names(knz_001d2)=c("sp", "yr", "sub", "abund")
    
    #add a random character and factor column
    knz_001d2$randcharacter<-sample(letters, size = nrow(knz_001d2), replace = T)
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
    
    # test that calculation from turnover is correct and does not regress
    myresults<-turnover(df=knz_001d, time.var="year", species.var="species", abundance.var="abundance",  replicate.var="subplot", 
                        metric="total")
    expect_equal(class(myresults), "data.frame")
    expect_equal(nrow(myresults), 460)
    expect_equal(myresults[460,1], 0.47368421, tolerance=0.00001)
    expect_equal(sum(myresults[,1]), 116.2359, tolerance=0.00001)
    
    # 
    
    #test that works regardless of whether parameter is specified or just ordered
    myresults2<-turnover(df=knz_001d, replicate.var="subplot",
                         metric="total")    
    expect_identical(myresults, myresults2)
    
    #test that total is the default metric
    myresults3 <- turnover(df=knz_001d, time.var="year", species="species", abundance="abundance", replicate.var="subplot")
    expect_identical(myresults, myresults3)
    
    #test that works with different column nmaes
    myresults4 <- turnover(df=knz_001d2, time.var="yr", species="sp", abundance="abund", replicate="sub", metric="total")
    expect_equal(sum(myresults$total), sum(myresults4$total))
    
    #test that works with different column orders if names specified
    myresults5 <- turnover(df=knz_001d2, time.var="yr", abundance="abund", replicate="sub", metric="total", species="sp")
    expect_equal(sum(myresults$turnover), sum(myresults5$turnover))
    
    #test that gives an error if abundance is a character or factor
    expect_error(turnover(df=knz_001d2, time.var="yr", species.var="sp", abundance.var="randcharacter", replicate.var="sub"))
    expect_error(turnover(df=knz_001d2, time.var="yr", species.var="sp", abundance.var="randfactor", replicate.var="sub"))
    
    # test that stops if more than one record for a species within a year, in one replicate
    dat4 = rbind(dat1, dat1[nrow(dat1),])
    expect_error(turnover(dat4, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot"))
    
    dat5 = rbind(knz_001d, knz_001d[nrow(knz_001d),], knz_001d[1,])
    expect_error(turnover(dat5, time.var="year", species.var="species", abundance.var="abundance", replicate.var="subplot"))
    
    
})

