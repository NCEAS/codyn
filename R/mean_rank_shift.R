###Mean Rank Shift
library(reshape2)
library(plyr)
library(reshape)

#using Scott's version of SEV004 veg transect data
SEV004.sc <- read.csv("SEV004_sc.csv", header = TRUE, strip.white = TRUE)

veg <- SEV004.sc


###melting spp matrix
veg_long <- melt(veg, id.vars="year", value.name ="cover", variable.name = "species")
head(veg_long)


#subsetting data to only present data
mydat <- subset(veg_long, veg_long$cover>0)
#creating empty result dataframe (year and mrs) that will be populated via loop
pairwiseMRS <- data.frame(cbind(year=numeric(0), mrs=numeric(0)))
#creating year vector so not to refer to column in dataframe 
yr <- sort(unique(mydat$year))


#loop
#iterations = length of study -1 due to pairwise comparison of years
for (i in 2:(length(yr))){
	#subsetting via subsequent two year window
	subber <- subset(mydat, mydat$year == yr[i] | mydat$year == yr[i-1])
	
	#ID species present in both years
	subber$pres <- 1
	subbercount <- aggregate(pres~species, data=subber, sum)
	subboth <- subset(subbercount, subbercount$pres == 2) 
	
	#create dataframe with only the pairwise present species
	#merging original dataframe with list of pairwise present species
	subber2 <- merge(subber[1:3], subboth[1])
	#ordering dataframe by year, then by decresing cover 
	subber2 <- subber2[with(subber2, order(year, -cover)), ]
	#ranking species by cover per year (this is best done after removing non-present and non-pairwise species)
  subber2$rank <- with(subber2, ave(cover, year, FUN=function (x) rank(-x, ties.method ="min")))
	#taking the absolute difference in rank for each species between years 
	subber2$rankdiff <- with(subber2, ave(rank, species, FUN=function (x) abs(diff(x))))
	
	#initial subset of two year dataframe is now reduced to one year (ranks were replicated for each year, so it does not matter which year)
	subber3 <- subset(subber2, subber2$year == yr[i])
	#dropping dataframe format and taking the cross species rank average for the pairwise years 
	mrs <- mean(subber3$rankdiff)
	
	#populating empty result dataframe
	#taking appropriate year (second of the pair)
	year <- yr[i]
	#creating temporary tuple by binding year with mrs
	temp <- data.frame(cbind(year, mrs))
	#populating
	pairwiseMRS <- rbind(pairwiseMRS, temp)
}
	
	

	
	
