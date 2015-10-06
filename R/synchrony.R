#' A function to calculate species synchrony over time within multiple replicates
#'
#' @param df A dataframe containing replicate, time, species and abundance columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @param metric The synchrony metric to return. The default, "Loreau", returns synchrony as calculated by Loreau and de Mazancourt 2008.
#'        The alternative, "Gross", returns synchrony as calculated by Gross et al. 2014
#' @param replicate.var The name of the replicate column from df. Defaults to NA.
#' @return output The degree of species synchrony. If "Loreau", 1 is perfect synchrony and 0 is perfect asynchrony. 
#'        If "Gross", 1 is perfect synchrony and -1 is perfect asynchrony.
#' @examples 
#' data(knz_001d)
#' synchrony(knz_001d[knz_001d$subplot=="A_1",]) # for one subplot
#' synchrony(knz_001d, replicate.var = "subplot") # across all subplots
#' synchrony(knz_001d, replicate.var = "subplot", metric="Gross") # With Gross et al. 2014 metric.
#' @export
synchrony<-function(df, time.var="year", species.var="species", abundance.var="abundance", metric="Loreau", replicate.var=NA) {
  # check to see if there are actual replicates without specifying replicate.var
  if(is.na(replicate.var)==TRUE){
    check_single_onerep(df, time.var, species.var)  
    output <- synch_onerep(df, time.var, species.var, abundance.var, metric)
    } else {
      df[replicate.var] <- if(is.factor(df[[replicate.var]])) { 
        factor(df[[replicate.var]])
      } else {
        df[replicate.var]
      }
      check_single(df, time.var, species.var, replicate.var)
      X <- split(df, df[replicate.var])
      out <- lapply(X, FUN=synch_onerep, time.var, species.var, abundance.var, metric)
      reps <- unique(df[replicate.var])
      output <- cbind(reps, do.call("rbind", out))
      names(output) = c(replicate.var, "synchrony")
      row.names(output) <- NULL
    }
    return(output)
}

############################################################################
#
# Private functions: these are internal functions not intended for reuse.  
# Future package releases may change these without notice. External callers
# should not use them.
#
############################################################################

#' A function to calculate species synchrony over time within one replicate
#'
#' @param df A dataframe containing rep, time, species and abundance columns
#' @param time.var The name of the time column from df
#' @param species.var The name of the species column from df
#' @param abundance.var The name of the abundance column from df
#' @param metric The synchrony metric to return. The default, "Loreau", returns synchrony as calculated by Loreau and de Mazancourt 2008.
#'        The alternative, "Gross", returns synchrony as calculated by Gross et al. 2014
#' @return output The degree of species synchrony. If "Loreau", 1 is perfect synchrony and 0 is perfect asynchrony. 
#'        If "Gross", 1 is perfect synchrony and -1 is perfect asynchrony.

synch_onerep <- function(df, time.var, species.var, abundance.var, metric=c("Loreau", "Gross")) {
    metric = match.arg(metric) # for partial argument matching
    # check to make sure abundance is numeric data
    if(!is.numeric(df[,abundance.var])) { stop("Abundance variable is not numeric") }
    
    #remove any species that were never present. 
    df <- subset(df, abundance.var>0)
    #fill in 0s
    spplist <- unique(df[species.var])
    yearlist <- unique(df[time.var])
    fulllist <- expand.grid(species.var=spplist[[species.var]], time.var=yearlist[[time.var]])
    # recapture original names
    names(fulllist) = c(species.var, time.var)
    df2 <- merge(df[c(species.var, time.var, abundance.var)], fulllist, all.y = T)
    df2[is.na(df2)] <- 0  
    
    if(metric=="Loreau"){
      #calculate community variance
      XTformula <- as.formula(paste(abundance.var, "~", time.var, sep=""))
      XT <- aggregate(XTformula, data=df2, sum)
      #do this within rep
      varXT<-var(XT[abundance.var])
    
      #calculate species variance
      sdSppformula <- as.formula(paste(abundance.var, "~", species.var, sep=""))
      sdSpp <- aggregate(sdSppformula, data=df2, sd)
      varSpp <- sum((sdSpp[abundance.var])) * sum(sdSpp[abundance.var])
      
      #calculate synchrony
      synchrony <- as.numeric(varXT/varSpp)
      
    } else {
      if(metric=="Gross"){
        corout<-as.data.frame(cbind(species.var= as.character(), "sppcor"=as.numeric()))
        
        for (i in 1:nrow(spplist)){
          myspp<-as.character(spplist[[1]][i])
          focalspp<-df2[which(df2[species.var] == myspp),]
          com.focalspp<-transpose_community(focalspp, time.var, species.var, abundance.var)
          otherspp<-df2[which(df2[species.var] != myspp),]
          com.otherspp<-transpose_community(otherspp, time.var, species.var, abundance.var)
          agg.otherspp<-(rowSums(com.otherspp))
          sppcor<-cor(agg.otherspp, com.focalspp)
          ##need to add a default value for species that do not fluctuate at all
          subout<-as.data.frame(cbind(myspp, sppcor))
          names(subout)=c(species.var, "sppcor")
          subout$sppcor<-as.numeric(as.character(subout$sppcor))    
          corout<-rbind(corout, subout)
        }
        #average correlation for the community
        synchrony <- mean(corout$sppcor)
      }
    }
    
    return(synchrony)
}

