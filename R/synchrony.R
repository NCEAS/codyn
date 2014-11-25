#' A function to calculate species synchrony over time within multiple replicates
#'
#' @param data1 A dataframe containing rep, year, species and abundance columns
#' @param rep The name of the replicate column from data1
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return output The degree of species synchrony, where 1 is perfect synchrony and 0 is perfect asynchrony
#' @import reshape
#' @export
synchrony<-function(data1, rep, species, year, abundance) {
    X <- split(data1, data1[rep])
    out<-lapply(X, FUN=synch_onerep, species, year, abundance)
    reps<-unique(data1[rep])
    output<-cbind(reps, do.call("rbind", out))
    names(output)=c(rep, "synchrony")
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
#' @param data1 A dataframe containing rep, year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @return output The degree of species synchrony, where 1 is perfect synchrony and 0 is perfect asynchrony
synch_onerep<-function(data1, species, year, abundance) {
    #remove any species that were never present
    data1<-subset(data1, abundance>0)
    #fill in 0s
    spplist<-(unique(data1[species]))
    yearlist<-(unique(data1[year]))
    fulllist<-expand.grid(species=spplist[[species]], year=yearlist[[year]])
    data2<-merge(data1[c(species, year, abundance)], fulllist, all=T)
    data2[is.na(data2)]<-0  
  
    #calculate community variance
    XTformula<-as.formula(paste(abundance, "~", year, sep=""))
    XT<-aggregate(XTformula, data=data2, sum )
    #do this within rep
    varXT<-var(XT[abundance])
    
    #calculate species variance
    sdSppformula<-as.formula(paste(abundance, "~", species, sep=""))
    sdSpp<-aggregate(sdSppformula, data=data2, sd)
    varSpp<-sum((sdSpp[abundance]))*sum(sdSpp[abundance])
    
    #calculate synchrony
    synchrony<-as.numeric(varXT/varSpp)
    return(synchrony)
}

