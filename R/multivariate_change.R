#'@title Multivariate changes in composition and dispersion
#' @description 
#' @param df A data frame containing time, species, abundance and replicate columns and an optional column of treatment
#' @param time.var The name of the time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var the neame of the optional treatment column


multivariate_change <- function(df, time.var, species.var, abundance.var, replicate.var, treatment.var = NULL){
  
  df <- as.data.frame(df)
  
    if(is.null(treatment.var)){
  
      mult_com_change <- mult_change(df, time.var, species.var, abundance.var, replicate.var)
      
      }
  
  else {
    
  # calculate change for each treatment
  df[[treatment.var]] <- as.character(df[[treatment.var]])
  df <- df[order(df[[treatment.var]]),]
  X <- split(df, df[treatment.var])
  out <- lapply(X, FUN = mult_change, time.var, species.var, abundance.var, replicate.var)
  ID <- unique(names(out))
  out <- mapply(function(x, y) "[<-"(x, treatment.var, value = y) ,
                out, ID, SIMPLIFY = FALSE)
  mult_com_change <- do.call("rbind", out)
  
  
  }
  
  return(mult_com_change)
}

#######Private Functions

mult_change <- function(df, time.var, species.var, abundance.var, replicate.var){
  require(vegan)

  #get years
timestep <- sort(unique(df[[time.var]]))

#transpose data
df2<-subset(df, select = c(time.var, species.var, abundance.var, replicate.var))
df2$id <- paste(df2[[time.var]], df2[[replicate.var]], sep="##")
species<-codyn:::transpose_community(df2, 'id', species.var, abundance.var)
species$id <- row.names(species)
speciesid <- do.call(rbind.data.frame, strsplit(species$id, split="##"))
colnames(speciesid)[1] <- "time_forfixxyz"
colnames(speciesid)[2] <- replicate.var
speciesid[[time.var]]<-as.numeric(as.character(speciesid$time_forfixxyz))
speciesid.2 <- subset(speciesid, select = -time_forfixxyz)
species2 <- cbind(speciesid.2, species)
species3 <- subset(species2, select = -id)

#calculate bray-curtis dissimilarities
bc <- vegdist(species3[,3:ncol(species3)], method="bray")

#calculate distances of each plot to year centroid (i.e., dispersion)
disp <- betadisper(bc, species3[[time.var]], type="centroid")

#getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
cent_dist <- as.matrix(vegdist(disp$centroids, method="euclidean"))

##extracting only the comparisions we want year x to year x=1.
###year x+1
cent_dist_yrs <- data.frame(time1 = timestep[1:length(timestep)-1], time2 = timestep[2:length(timestep)],
                         composition_change = diag(cent_dist[2:nrow(cent_dist), 1:(ncol(cent_dist)-1)]))

#collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
disp2 <- data.frame(time=species2[[time.var]],
                 dist = disp$distances)

myformula <- as.formula(paste("dist", "~", "time"))
disp2.2<-aggregate(myformula, mean, data=disp2)

##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
disp_yrs <- data.frame(time2 = timestep[2:length(timestep)],
                    dispersion_change = diff(disp2.2$dist))

#merge together change in mean and dispersion data
distances <- merge(cent_dist_yrs, disp_yrs, by="time2")
distances$time_pair<-paste(distances$time1, distances$time2, sep="-")

distances<-subset(distances, select = c("time_pair", "composition_change", "dispersion_change"))
colnames(distances)[1]<-paste(time.var, "pair", sep="_")

return(distances)
}
