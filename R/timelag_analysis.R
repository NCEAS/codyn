require(vegan)
require(car)
require(plotrix)

# creating distance matrix from species matrix
data.dm <- vegdist(spp_matrix, method="euclidean",
								diag=FALSE, upper=FALSE, na.rm=TRUE)
data.dm <- as.matrix(data.dm)	

# get the lags from distance matrix
get_lags = function(data.dm){

# label each row and each column.
rownums = row(data.dm)
colnums = col(data.dm)

# A mini-function to get all the elements that are lagged by i time steps,
# then put them into a long-form matrix with i in the first column and
# the matrix value in the second column.
get_lag_i = function(i){
    cbind(lag = i, value = data.dm[rownums == (colnums + i)])
  }

# apply get_lag_i to all lags from 1 to n-1
# replace n with number of columns 
  lag_list = lapply(
    1:(n-1),
    get_lag_i
  )

# squash all the lags from the list into one long-form matrix
  do.call(rbind, lag_list)
}

# final: matrix converted to long-form
dm.long <- data.frame(get_lags(data.dm))

# regression
lm.dm <- lm(value ~ lag, data=dm.long)
summary(lm.dm)



# plotting regression graph
# replace "xMin", "xMax", "yMin", "yMax" with your values 
par(mar=c(5, 5, 1, 1), xpd=TRUE)
plot(dm.long$value, dm.long$lag, type="n",
     xlab="time lag",
     ylab="distance",
		xlim=c(xMin,xMax), ylim=c(yMin,yMax))

# plotting points
points(dm.long$lag, dm.long$value, cex = 1.0, pch=16, col="#99CC66")

# plotting regression line
ablineclip(lm.dm, x1=0, x2=25, type="l", lty= 1, lwd=3.5, col="#000000")
# a=yStart, b=yEnd, x1=xStart, x2= xEnd
