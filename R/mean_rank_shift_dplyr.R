#' Function for obtaining an example dataset
#'
#' @

getexample <- function(){
  dat <- read.csv("data/knz_001d.csv", header = TRUE, strip.white = TRUE)
  subset(dat$subplot == "A_1", drop = TRUE)
}




pairwiseMRS <- data.frame(cbind(year=numeric(0), mrs=numeric(0)))
yr <- sort(unique(veg$year))

for (i in 2:(length(yr))){
	veg %>%
	filter(cover>0) %>%
	filter(year == yr[i] | year == yr[i-1]) %>%

	mutate(pres = 1)
	subbercpount <- summarise(pres~species, sum)
	subboth <- subset(subbercount, subbercount$pres == 2)

}

clean_veg <- veg %>%
	gather(species, cover, -year) %>%
  arrange(year, cover) %>%
	filter(cover>0) %>%
	summarise(group_by(year), mean(cover))


percent_rank(veg$cover)
