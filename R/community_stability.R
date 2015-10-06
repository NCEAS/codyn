#' @title Community Stability
#' @description   A function to calculate community stability over time
#' 
#' @param df the data frame to use in the calculation
#' @param time.var The name of the time column from df
#' @param abundance.var The name of the abundance column from df
#' @param replicate.var The name of the replicate column from df
#' @return The stability of community abundance, calculated as mean/standard deviation returns a single number or a dataframe with the following attributes:
#' \itemize{
#'  \item{replicate}{replicate name}
#'  \item{stability}{community stability for replicate}
#'  }
#' @details 
#' The input dataframe needs to contain columns for time and abundance; 
#' @export

community_stability<-function(df, time.var="year", abundance.var="abundance", replicate.var=NA){
    stopifnot(is.numeric(df[[time.var]]))
    stopifnot(is.numeric(df[[abundance.var]]))
    df<-df[which(df[[abundance.var]]>0),]

    if(is.na(replicate.var)==TRUE) {
      
        #sum abundance within a year
        aggform<-as.formula(paste(abundance.var, "~", time.var, sep=""))
        data2<-aggregate(aggform, data=df, sum)
        output<-stability_onerep(data2, abundance.var)
        
    } else {
    
        df<-df[order(df[[replicate.var]]),]  
        df[replicate.var]<-if(is.factor(df[[replicate.var]])==TRUE){factor(df[[replicate.var]])} else {df[replicate.var]}
  
        #sum abundance within a replicate and year
        aggform<-as.formula(paste(abundance.var, "~", replicate.var, "+", time.var, sep=""))
        data2<-aggregate(aggform, data=df, sum)
        data2<-data2[order(data2[[replicate.var]]),]  
        X<-split(data2, data2[replicate.var])
        out<-lapply(X, stability_onerep, abundance.var)
        reps<-unique(data2[replicate.var])
        output<-cbind(reps, do.call("rbind", out))
        names(output)=c(replicate.var, "stability")
        output<-subset(output, !is.na(output$stability))

    }
    
    row.names(output)<-NULL
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
#' @param df A dataframe containing x column
#' @param x The column to calculate stability on
#' @return Stability of x, calculated as the mean/sd
stability_onerep<-function(df,  x){
  return(mean(df[[x]], na.rm=T)/sd(df[[x]], na.rm=T))
}

