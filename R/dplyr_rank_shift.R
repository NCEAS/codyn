

# read in data ------------------------------------------------------------

veg <- read.csv("SEV004_sc.csv")


# libraries ---------------------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)

long_veg <- veg %>%
  gather(species, cover, -year) %>%
  arrange(year, cover) %>%
  filter(cover > 0)

## make list of years
years <- unique(long_veg$year) %>% sort

## make pairs of years
timechange <- data.frame(y1 = years[-length(years)],
                         y2 = years[-1])


timechange %>%
  mutate(yearpair = paste0(y1, "-", y2)) %>%
  left_join(select(long_veg, y1 = year, sp_y1 = species, cov_y1 = cover)) %>%
  left_join(select(long_veg, y2 = year, sp_y2 = species, cov_y2 = cover)) %>%
  group_by(yearpair) %>%
  do(inner_join(select(., species = sp_y1),
                select(., species = sp_y2)))


comboyears <- timechange %>%
  mutate(yearpair = paste0(y1, "-", y2)) %>%
  group_by(yearpair) %>%
  do(y1_list = filter(long_veg, year == .$y1) %>%
       select(species, cov1 = cover),
     y2_list = filter(long_veg, year == .$y2) %>%
       select(species, cov2 = cover)) %>%
  do(data.frame(yearpair = .$yearpair,
                inner_join(.$y1_list, .$y2_list)))

rankchange <- comboyears %>%
  ungroup %>%
  group_by(yearpair) %>%
  mutate(rank1 = percent_rank(cov1),
         rank2 = percent_rank(cov2),
         deltarank = abs(rank1 - rank2))

rankchange %>%
  ungroup %>%
  group_by(species) %>%
  summarize(mean_deltarank = mean(deltarank)) %>%
  arrange(mean_deltarank) %>%
  ggplot(aes(x = order(mean_deltarank, decreasing = TRUE), y = mean_deltarank)) + geom_point()

ggsave("nicefig.png")

# + stat_summary()


## bind years together

#
#
# yrs<-unique(mydat$year)
# for (i in 2:length(yrs)){
#   subber<-subset(mydat, mydat$year==yrs[i] | mydat$year==yrs[i-1])
#   ##ID species present in both years
#   subberpres<-subset(subber, subber$cover>0)
#   subberpres$pres<-1
#   subcount<-aggregate(pres~species, data=subberpres, sum)
#   subboth<-subset(subcount, subcount$pres==2)
#   #create dataframe with only the pairwise present species
#   subber2<-merge(subber, subboth[1])
#   subber2 <- subber2[with(subber2, order(year, -cover)),]
#   subber2$rank <- with(subber2, ave(cover, year, FUN=function (x) rank(-x, ties.method ="min")))
#   subber2$rankdiff <- with(subber2, ave(rank, species, FUN=function (x) abs(diff(x))))
#   subber3<-subset(subber2, subber2$year==yrs[i])
#   mrs<-mean(subber3[,"rankdiff"])
#   year<-yrs[i]
#   temp<-data.frame(cbind(year, mrs))
#   pairwiseMRS<-rbind(pairwiseMRS, temp)
# }
#
#
#
#
#   SEV004_long <- melt(SEV004.sc, id.vars=c("year"), value.name = "cover", variable.name = "species")
# head(SEV004_long)

werk <- 20

test_function <- function(A = 20){
  A + werk
}
