#' @title pool replicates into treatments and add ranks
#' @description pool replicates into treatments and add ranks
#' @param df A data frame containing an optional time column and requred species, abundance, replicate and optional treatment columns and optional block column
#' @param time.var The name of the optional time column 
#' @param species.var The name of the species column 
#' @param abundance.var The name of the abundance column 
#' @param replicate.var The name of the replicate column 
#' @param treatment.var The name of the optional treatment column
#' @param block.var The name of the optional block column


RAC_difference <- function(df, time.var=NULL, species.var, abundance.var, replicate.var, treatment.var=NULL, pool="NO", block.var=NULL){

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
  ######privte functions
  SERSp <- function(df, rank.var1, rank.var2, abundance.var1, abundance.var2){
    
    ## Remove 0s
    df <- subset(df, df[[abundance.var1]]!=0 | df[[abundance.var2]]!=0)
    
    ## Remove instances where abundance.var is NA, probably not necessary since we are doing over space and not time but in as a check #why is this necessary?
    df <- subset(df, !is.na(df[[abundance.var1]]) & !is.na(df[[abundance.var2]]))
    
    #ricness and evenness differences
    s_t1 <- S(df[[abundance.var1]])
    e_t1 <- EQ(as.numeric(df[[abundance.var1]]))
    s_t2 <- S(df[[abundance.var2]])
    e_t2 <- EQ(as.numeric(df[[abundance.var2]]))
    
    sdiff <- abs(s_t1-s_t2)/nrow(df)
    ediff <- abs(e_t1-e_t2)/nrow(df)
    
    spdiff <- df[df[[abundance.var1]] == 0|df[[abundance.var2]] == 0,]#i don't think this will work if all species are shared
    spdiffc <- nrow(spdiff)/nrow(df)
    
    mrsc_diff <- mean(abs(df[[rank.var1]]-df[[rank.var2]])/nrow(df))
    
    metrics <- data.frame(richness_diff=sdiff, evenness_diff=ediff, rank_diff=mrsc_diff, species_diff = spdiffc)
    
    return(metrics)
  }
