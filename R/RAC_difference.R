#' @title Rank Abundance Curve Differences
#' @description Calculates differences between two samples for four comparable aspects of rank abundance curves (richness, evenness, rank, species composition). There are three ways differences can be calculated. 1) Between treatments within a block (note: block.var and treatment.var need to be specified). 2) Between treatments, pooling all replicates into a single species pool (note: pool = TRUE, treatment.var needs to be specified, and block.var will be NULL). 3) All pairwise combinations between all replicates (note: block.var = NULL, pool = FALSE and specifying treatment.var is optional. If treatment.var is specified, the treatment that each replicate belongs to will also be listed in the output).
#' @param df A data frame containing a species, abundance, and replicate columns and optional time, treatment, and block columns
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var The name of the optional treatment column
#' @param block.var The name of the optional block column
#' @param pool An argument to allow abundance values to be pooled within a treatment. The default value is "FALSE", a value of "TRUE" averages abundance of each species within a treatment at a given time point.
#' @return The RAC_difference function returns a data frame with the following attributes:
#' \itemize{
#'  \item{time.var: }{A column that has the same name and type as the time.var column, if time.var is specified.}
#'  \item{block.var: }{A column that has same name and type as the block.var column, if block.var is specified.}
#'  \item{replicate.var: }{A column that has same name and type as the replicate.var column, represents the first replicate being compared. Note, a replicate column will be returned only when pool is FALSE or block.var = NULL.}
#'  \item{replicate.var2: }{A column that has the same type as the replicate.var column, and is named replicate.var with a 2 appended to it, represents the second replicate being compared. Note, a replicate.var column will be returned only when pool is FALSE and block.var = NULL.}
#'  \item{treatment.var: }{A column that has same name and type as the treatment.var column, represents the first treatment being compared. A treatment.var column will be returned when pool is TRUE or block.var is present, or treatment.var is specified.}
#'  \item{treatment.var2: }{A column that has the same type as the treatment.var column, and is named treatment.var with a 2 appended to it, represents the second treatment being compared. A treatment.var column will be returned when pool is TRUE or block.var is present, or treatment.var is specified.}
#'  \item{richness_diff: }{A numeric column that is the difference between the compared samples (treatments or replicates) in species richness divided by the total number of species in both samples.}
#'  \item{evenness_diff: }{A numeric column of the difference between the compared samples (treatments or replicates) in evenness (measured using the EQ metric) divided by the total number of species in both samples.}
#'  \item{rank_diff: }{A numeric column of the average difference between the compared samples (treatments or replicates) in species' ranks divided by the total number of species in both samples. Species that are not present in both samples are given the S+1 rank in the sample it is absent in, where S is the number of species in that sample.}
#'  \item{species_diff: }{A numeric column of the number of species that are different between the compared samples (treatments or replicates) divided by the total number of species in both samples. This is equivelant to the Jaccard Index.}
#' }
#' @references Avolio et al. OUR PAPER
#' @examples 
#' data(pplots)
#' # With block and no time
#' df <- subset(pplots, year == 2002 & block < 3)
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                treatment.var = 'treatment',
#'                block.var = "block",
#'                replicate.var = "plot")
#' # With blocks and time
#' df <- subset(pplots, year < 2004 & block < 3)
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                treatment.var = 'treatment',
#'                block.var = "block",
#'                replicate.var = "plot",
#'                time.var = "year")
#' # Pooling by treatment no time
#' df <- subset(pplots, year == 2002)
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                treatment.var = 'treatment',
#'                pool = "YES",
#'                replicate.var = "plot")
#' # Pooling by treatment with time
#' df <- subset(pplots, year < 2004)
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                treatment.var = 'treatment',
#'                pool = "YES",
#'                replicate.var = "plot",
#'                time.var = "year")
#' # All pairwise replicates with treatment and no time
#' df <- subset(pplots, year == 2002 & plot %in% c(6, 25, 32))
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                replicate.var = "plot",
#'                treatment.var = "treatment")
#' # All pairwise replicates with treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                replicate.var = "plot",
#'                time.var = "year",
#'                treatment.var = "treatment")
#' # All pairwise replicates without treatment and no time
#' df <- subset(pplots, year == 2002 & plot %in% c(6, 25, 32))
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                replicate.var = "plot")
#' # All pairwise replicates without treatment
#' df <- subset(pplots, year < 2004 & plot %in% c(6, 25, 32))
#' RAC_difference(df = df,
#'                species.var = "species",
#'                abundance.var = "relative_cover",
#'                replicate.var = "plot",
#'                time.var = "year")
#' @export


RAC_difference <- function(df, time.var=NULL, species.var, abundance.var, replicate.var, treatment.var=NULL, pool="NO", block.var=NULL){

  # check no NAs in abundance column
  if(any(is.na(df[[abundance.var]]))) stop("Abundance column contains missing values")
  
  if(!is.null(block.var)){
    
    if(length(unique(df[[block.var]]))*length(unique(df[[treatment.var]])) != length(unique(df[[replicate.var]]))) stop("There is not one replicate per treatment in a block")
    
    myperms<-trt_perms(df, treatment.var)
    
     if(is.null(time.var)){
       
       #rank species in each replicate
       rep_trt <- unique(subset(df, select = c(replicate.var, treatment.var, block.var)))
       rankdf <- add_ranks_replicate(df, time.var=NULL, species.var, abundance.var, replicate.var)
       rankdf1 <- merge(rankdf, rep_trt, by=replicate.var)
      
       ## Create a second rankdf with a renamed treatment.var column
       rankdf2 <- rankdf1
       rankdf2[[paste(treatment.var, "2", sep = "")]] <- rankdf2[[treatment.var]]
       rankdf2[[treatment.var]] <- NULL
       
       ## Merge rankdf with all possible permutations of treatment combinations
       rankdfall <- merge(rankdf1, myperms, all.y = T)
       
       ## Merge the data together (for all possible permutations of treatments) within a block
       ranktog <- merge(rankdfall, rankdf2, by=c(species.var, block.var, paste(treatment.var, "2", sep="")))
      
       ## Create a variable to split on (block as well as unique treatment combos)
       ranktog$splitvariable = paste(ranktog[[block.var]], ranktog[[treatment.var]], ranktog[[paste(treatment.var, "2", sep = "")]], sep="##")
       
       ## Split the dataframe
       X <- split(ranktog, ranktog$splitvariable)
       
       ## Apply the  RAC function for differences
       out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
       ID <- unique(names(out))
       out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                     out, ID, SIMPLIFY = FALSE)
       output <- do.call("rbind", out)  
       
       ## Add in the identifying column names
       outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
       names(outnames) = c(block.var, treatment.var, paste(treatment.var, "2", sep=""))
       output$splitvariable <- NULL
       output <- cbind(outnames, output)
       
       }
    
    else{
      #rank species in each replicate
      rep_trt <- unique(subset(df, select = c(replicate.var, treatment.var, block.var)))
      rankdf <- add_ranks_replicate(df, time.var, species.var, abundance.var, replicate.var)
      rankdf1 <- merge(rankdf, rep_trt, by=replicate.var)
      
      ## Create a second rankdf with a renamed treatment.var column
      rankdf2 <- rankdf1
      rankdf2[[paste(treatment.var, "2", sep = "")]] <- rankdf2[[treatment.var]]
      rankdf2[[treatment.var]] <- NULL
      
      ## Merge rankdf with all possible permutations of treatment combinations
      rankdfall <- merge(rankdf1, myperms, all.y = T)
      
      ## Merge the data together (for all possible permutations of treatments) within a block
      ranktog <- merge(rankdfall, rankdf2, by=c(time.var, species.var, block.var, paste(treatment.var, "2", sep="")))
      
      ## Create a variable to split on (block, time, and unique treatment combos)
      ranktog$splitvariable = paste(ranktog[[time.var]], ranktog[[block.var]], ranktog[[treatment.var]], ranktog[[paste(treatment.var, "2", sep = "")]], sep="##")
      
      ## Split the dataframe
      X <- split(ranktog, ranktog$splitvariable)
      
      ## Apply the  RAC function for differences
      out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
      ID <- unique(names(out))
      out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                    out, ID, SIMPLIFY = FALSE)
      output <- do.call("rbind", out)  
      
      ## Add in the identifying column names
      outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
      names(outnames) = c(time.var, block.var, treatment.var, paste(treatment.var, "2", sep=""))
      output$splitvariable <- NULL
      output <- cbind(outnames, output)
      
    }
  }
  else{
    if(pool=="YES"){
      
      myperms<-trt_perms(df, treatment.var)
      
    if(is.null(time.var)){
      ##pool data into treatment and rank
      rankdf<-pool_replicates(df, time.var, species.var, abundance.var, replicate.var, treatment.var)
      
      ## Create a second rankdf with a renamed treatment.var column
      rankdf2 <- rankdf
      rankdf2[[paste(treatment.var, "2", sep = "")]] <- rankdf2[[treatment.var]]
      rankdf2[[treatment.var]] <- NULL
      
      ## Merge rankdf with all possible permutations of treatment combinations
      rankdfall <- merge(rankdf, myperms, all.y = T)
      
      ## Merge the data together (for all possible permutations of treatments) within a block
      ranktog <- merge(rankdfall, rankdf2, by=c(species.var, paste(treatment.var, "2", sep="")))
      
      ## Create a variable to split on (block as well as unique treatment combos)
      ranktog$splitvariable = paste(ranktog[[treatment.var]], ranktog[[paste(treatment.var, "2", sep = "")]], sep="##")
      
      ## Split the dataframe
      X <- split(ranktog, ranktog$splitvariable)
      
      ## Apply the  RAC function for differences
      out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
      ID <- unique(names(out))
      out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                    out, ID, SIMPLIFY = FALSE)
      output <- do.call("rbind", out)  
      
      ## Add in the identifying column names
      outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
      names(outnames) = c(treatment.var, paste(treatment.var, "2", sep=""))
      output$splitvariable <- NULL
      output <- cbind(outnames, output)
    }
      else{
        rankdf<-pool_replicates(df, time.var, species.var, abundance.var, replicate.var, treatment.var)
        
        ## Create a second rankdf with a renamed treatment.var column
        rankdf2 <- rankdf
        rankdf2[[paste(treatment.var, "2", sep = "")]] <- rankdf2[[treatment.var]]
        rankdf2[[treatment.var]] <- NULL
        
        ## Merge rankdf with all possible permutations of treatment combinations
        rankdfall <- merge(rankdf, myperms, all.y = T)
        
        ## Merge the data together (for all possible permutations of treatments) within a block
        ranktog <- merge(rankdfall, rankdf2, by=c(time.var, species.var, paste(treatment.var, "2", sep="")))
        
        ## Create a variable to split on (block, time, and unique treatment combos)
        ranktog$splitvariable = paste(ranktog[[time.var]], ranktog[[treatment.var]], ranktog[[paste(treatment.var, "2", sep = "")]], sep="##")
        
        ## Split the dataframe
        X <- split(ranktog, ranktog$splitvariable)
        
        ## Apply the  RAC function for differences
        out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
        ID <- unique(names(out))
        out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                      out, ID, SIMPLIFY = FALSE)
        output <- do.call("rbind", out)  
        
        ## Add in the identifying column names
        outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
        names(outnames) = c(time.var, treatment.var, paste(treatment.var, "2", sep=""))
        output$splitvariable <- NULL
        output <- cbind(outnames, output)
      }
    }
      else{
       if(is.null(treatment.var)){
        myperms<-rep_perms(df, replicate.var) 
         
         if(is.null(time.var)){
           rankdf <- add_ranks_replicate(df, time.var=NULL, species.var, abundance.var, replicate.var)

           ## Create a second rankdf with a renamed replicate.var column
           rankdf2 <- rankdf
           rankdf2[[paste(replicate.var, "2", sep = "")]] <- as.factor(rankdf2[[replicate.var]])
           rankdf2[[replicate.var]] <- NULL
           
           ## Merge rankdf with all possible permutations of treatment combinations
           rankdfall <- merge(rankdf, myperms, all.y = T)
           
           ## Merge the data together
           ranktog <- merge(rankdfall, rankdf2, by=c(species.var, paste(replicate.var, "2", sep="")))
           
           ## Create a variable to split on (each replicate combination)
           ranktog$splitvariable = paste(ranktog[[replicate.var]], ranktog[[paste(replicate.var, "2", sep = "")]], sep="##")
           
           ## Split the dataframe
           X <- split(ranktog, ranktog$splitvariable)
           
           ## Apply the  RAC function for differences
           out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
           ID <- unique(names(out))
           out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                         out, ID, SIMPLIFY = FALSE)
           output <- do.call("rbind", out)  
           
           ## Add in the identifying column names
           outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
           names(outnames) = c(replicate.var, paste(replicate.var, "2", sep=""))
           output$splitvariable <- NULL
           output <- cbind(outnames, output)
         }
      else{
        rankdf <- add_ranks_replicate(df, time.var, species.var, abundance.var, replicate.var)
        
        ## Create a second rankdf with a renamed replicate.var column
        rankdf2 <- rankdf
        rankdf2[[paste(replicate.var, "2", sep = "")]] <- as.factor(rankdf2[[replicate.var]])
        rankdf2[[replicate.var]] <- NULL
        
        ## Merge rankdf with all possible permutations of treatment combinations
        rankdfall <- merge(rankdf, myperms, all.y = T)
        
        ## Merge the data together
        ranktog <- merge(rankdfall, rankdf2, by=c(species.var, time.var, paste(replicate.var, "2", sep="")))
        
        ## Create a variable to split on (each replicate combination)
        ranktog$splitvariable = paste(ranktog[[time.var]], ranktog[[replicate.var]], ranktog[[paste(replicate.var, "2", sep = "")]], sep="##")
        
        ## Split the dataframe
        X <- split(ranktog, ranktog$splitvariable)
        
        ## Apply the  RAC function for differences
        out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
        ID <- unique(names(out))
        out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                      out, ID, SIMPLIFY = FALSE)
        output <- do.call("rbind", out)  
        
        ## Add in the identifying column names
        outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
        names(outnames) = c(time.var, replicate.var, paste(replicate.var, "2", sep=""))
        output$splitvariable <- NULL
        output <- cbind(outnames, output)
      }
          } 
      else{
        
     myperms <- rep_perms(df, replicate.var)
        
        #create replicate treatment for reference
        rep_trt <- unique(subset(df, select = c(replicate.var, treatment.var)))
        rep_trt2 <- rep_trt
        rep_trt2[[paste(replicate.var, "2", sep = "")]] <- as.factor(rep_trt2[[replicate.var]])
        rep_trt2[[replicate.var]] <- NULL
        rep_trt2[[paste(treatment.var, "2", sep = "")]] <- as.factor(rep_trt2[[treatment.var]])
        rep_trt2[[treatment.var]] <- NULL
        
        if(is.null(time.var)){
          rankdf <- add_ranks_replicate(df, time.var=NULL, species.var, abundance.var, replicate.var)
          
          ## Create a second rankdf with a renamed replicate.var column
          rankdf2 <- rankdf
          rankdf2[[paste(replicate.var, "2", sep = "")]] <- as.factor(rankdf2[[replicate.var]])
          rankdf2[[replicate.var]] <- NULL
          
          ## Merge rankdf with all possible permutations of treatment combinations
          rankdfall <- merge(rankdf, myperms, all.y = T)
          
          ## Merge the data together
          ranktog <- merge(rankdfall, rankdf2, by=c(species.var, paste(replicate.var, "2", sep="")))
          
          ## Create a variable to split on (each replicate combination)
          ranktog$splitvariable = paste(ranktog[[replicate.var]], ranktog[[paste(replicate.var, "2", sep = "")]], sep="##")
          
          ## Split the dataframe
          X <- split(ranktog, ranktog$splitvariable)
          
          ## Apply the  RAC function for differences
          out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
          ID <- unique(names(out))
          out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                        out, ID, SIMPLIFY = FALSE)
          output <- do.call("rbind", out)  
          
          ## Add in the identifying column names
          outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
          names(outnames) = c(replicate.var, paste(replicate.var, "2", sep=""))
          output$splitvariable <- NULL
          output <- cbind(outnames, output)
          
          output <- merge(output, rep_trt, by = replicate.var)
          output <- merge(output, rep_trt2, by = paste(replicate.var, "2", sep=""))
        }
        else{
          rankdf <- add_ranks_replicate(df, time.var, species.var, abundance.var, replicate.var)
          
          ## Create a second rankdf with a renamed replicate.var column
          rankdf2 <- rankdf
          rankdf2[[paste(replicate.var, "2", sep = "")]] <- as.factor(rankdf2[[replicate.var]])
          rankdf2[[replicate.var]] <- NULL
          
          ## Merge rankdf with all possible permutations of treatment combinations
          rankdfall <- merge(rankdf, myperms, all.y = T)
          
          ## Merge the data together
          ranktog <- merge(rankdfall, rankdf2, by=c(species.var, time.var, paste(replicate.var, "2", sep="")))
          
          ## Create a variable to split on (each replicate combination)
          ranktog$splitvariable = paste(ranktog[[time.var]], ranktog[[replicate.var]], ranktog[[paste(replicate.var, "2", sep = "")]], sep="##")
          
          ## Split the dataframe
          X <- split(ranktog, ranktog$splitvariable)
          
          ## Apply the  RAC function for differences
          out <- lapply(X, FUN=SERSp, "rank.x", "rank.y", paste(abundance.var, ".x", sep = ""),paste(abundance.var, ".y", sep = "")) 
          ID <- unique(names(out))
          out <- mapply(function(x, y) "[<-"(x, "splitvariable", value = y) ,
                        out, ID, SIMPLIFY = FALSE)
          output <- do.call("rbind", out)  
          
          ## Add in the identifying column names
          outnames <- data.frame(do.call('rbind', strsplit(as.character(output$splitvariable),'##',fixed=TRUE)))
          names(outnames) = c(time.var, replicate.var, paste(replicate.var, "2", sep=""))
          output$splitvariable <- NULL
          output <- cbind(outnames, output)
          
          output <- merge(output, rep_trt, by = replicate.var)
          output <- merge(output, rep_trt2, by = paste(replicate.var, "2", sep=""))
      }
      }
      }
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

# A function to calculate RAC difference between two samples 
# @param df a dataframe
# @param rank.var1 the name of the rank column at time 1
# @param rank.var1 the name of the rank column at time 2
# @param abundance.var1 the name of the abundance column at time 1
# @param abundance.var2 the name of the abundance column at time 2
SERSp <- function(df, rank.var1, rank.var2, abundance.var1, abundance.var2){

  df <- subset(df, df[[abundance.var1]]!=0 | df[[abundance.var2]]!=0)
    
  df <- subset(df, !is.na(df[[abundance.var1]]) & !is.na(df[[abundance.var2]]))
    
  #ricness and evenness differences
  e_t1 <- S(df[[abundance.var1]])
  s_t1 <- EQ(as.numeric(df[[abundance.var1]]))
  s_t2 <- S(df[[abundance.var2]])
  e_t2 <- EQ(as.numeric(df[[abundance.var2]]))
    
  sdiff <- abs(s_t1-s_t2)/nrow(df)
  ediff <- abs(e_t1-e_t2)/nrow(df)
  
  #Jaccard Index or Number of species not shared  
  spdiff <- df[df[[abundance.var1]] == 0|df[[abundance.var2]] == 0,]
  spdiffc <- nrow(spdiff)/nrow(df)
  
  #Mean Rank Difference
  mrsc_diff <- mean(abs(df[[rank.var1]]-df[[rank.var2]])/nrow(df))
    
  metrics <- data.frame(richness_diff=sdiff, evenness_diff=ediff, rank_diff=mrsc_diff, species_diff = spdiffc)
    
  return(metrics)
  }
