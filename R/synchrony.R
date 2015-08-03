#' A function to calculate species synchrony over time within multiple replicates
#'
#' @param data1 A dataframe containing replicate, year, species and abundance columns
#' @param replicate The name of the replicate column from data1
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @param metric The synchrony metric to return. The default, "Loreau", returns synchrony as calculated by Loreau and de Mazancourt 2008.
#'        The alternative, "Gross", returns synchrony as calculated by Gross et al. 2014
#' @return output The degree of species synchrony. If "Loreau", 1 is perfect synchrony and 0 is perfect asynchrony. 
#'        If "Gross", 1 is perfect synchrony and -1 is perfect asynchrony.
#' @export
synchrony<-function(data1, replicate="replicate", species="species", year="year", abundance="abundance", metric="Loreau") {
  if(is.na(replicate)==TRUE){
    output<-synch_onerep(data1, species, year, abundance, metric)}else{
    data1[replicate]<-if(is.factor(data1[[replicate]])==TRUE){factor(data1[[replicate]])} else {data1[replicate]}
    X <- split(data1, data1[replicate])
    out<-lapply(X, FUN=synch_onerep, species, year, abundance, metric)
    reps<-unique(data1[replicate])
    output<-cbind(reps, do.call("rbind", out))
    names(output)=c(replicate, "synchrony")
    row.names(output)<-NULL
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
#' @param data1 A dataframe containing rep, year, species and abundance columns
#' @param year The name of the year column from data1
#' @param species The name of the species column from data1
#' @param abundance The name of the abundance column from data1
#' @param metric The synchrony metric to return. The default, "Loreau", returns synchrony as calculated by Loreau and de Mazancourt 2008.
#'        The alternative, "Gross", returns synchrony as calculated by Gross et al. 2014
#' @return output The degree of species synchrony. If "Loreau", 1 is perfect synchrony and 0 is perfect asynchrony. 
#'        If "Gross", 1 is perfect synchrony and -1 is perfect asynchrony.

synch_onerep<-function(data1, species, year, abundance, metric="Loreau") {
    #remove any species that were never present
    data1<-subset(data1, abundance>0)
    #fill in 0s
    spplist<-(unique(data1[species]))
    yearlist<-(unique(data1[year]))
    fulllist<-expand.grid(species=spplist[[species]], year=yearlist[[year]])
    data2<-merge(data1[c(species, year, abundance)], fulllist, all=T)
    data2[is.na(data2)]<-0  
  
    if(metric=="Loreau"){
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
    synchrony<-as.numeric(varXT/varSpp)} else{
      if(metric=="Gross"){
        corout<-as.data.frame(cbind(species=as.character(), "sppcor"=as.numeric()))
        for (i in 1:nrow(spplist)){
          myspp<-as.character(spplist[[1]][i])
          focalspp<-data2[which(data2[species] == myspp),]
          com.focalspp<-transpose_community(focalspp, species, year, abundance)
          otherspp<-data2[which(data2[species] != myspp),]
          com.otherspp<-transpose_community(otherspp, species, year, abundance)
          agg.otherspp<-(rowSums(com.otherspp))
          sppcor<-cor(agg.otherspp, com.focalspp)
          ##need to add a default value for species that do not fluctuate at all
          subout<-as.data.frame(cbind(myspp, sppcor))
          names(subout)=c(species, "sppcor")
          subout$sppcor<-as.numeric(as.character(subout$sppcor))    
          corout<-rbind(corout, subout)
        }
        
        #average correlation for the community
        synchrony<-mean(corout$sppcor)
      }
    }
    
    return(synchrony)
}

