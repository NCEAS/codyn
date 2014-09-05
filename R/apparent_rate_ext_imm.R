################
# Apparent rate of immigration and extinction

# format: pres/abs matrix

# count the number of times a species goes:
### "extinct": from 1 to 0
### immigrate: from 0 to 1

# If presence / absence is a vector over time, 
# the function "diff" (lag difference; the diff between a value and the preceding one) does a good job.

x<-c(0,1,1,0,0,0,1,0,1,0,1,0,0,1) # the species went extinct 4 times and re-colonised 5 times.

sum(diff(x)==1) #for immigration; it says how many time the lag difference = 1

sum(diff(x)==-1) # for extinction; how many time diff = -1
