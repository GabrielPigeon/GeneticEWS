# © Gabriel Pigeon
# Version 1          date: 2017-07-14
# Script for functions used in the analysis

# source()
library(parallel)
#############################################################################.


# function to calculate CV
CV <- function(x, na.rm){
  ave<-mean(x, na.rm=na.rm)
  dev<-sd(x, na.rm=na.rm)
  return((dev/ave))
}



##function for rolling mean
rolling_mean <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- mean(x[1:i], na.rm=T);
  }
  return(result);
}



##function for rolling sd
rolling_sd <- function(x){
  k = length(x);
  result = rep(0, k);
  for(i in 1 : k){
    result[i] <- sd(x[1:i], na.rm=T);
  }
  return(result);
}

# function to generate composite signals
makeComposite <- function(time, indices,cum.max=NULL){
  #   time=thedata$norm.day;indices = thedata[,9:13] ;cum.max=NULL
  if(!length(time)==nrow(indices)){stop("time vector does no match indices matrix")}
  indicesM <- as.matrix(indices)
  if(!is.numeric(indicesM)){stop("problem with indices")}
  
  indicesM[is.nan(indicesM)] <- NA
  
  indices <- as.data.frame(indicesM)
  results <- cbind(time,indices)
  #
  #change the results data frame from wide to long format
  kk<-reshape2::melt(results, id=1)
  
  ##a vector of the indicators present in the results data frame
  indicators<-unique(kk$variable)
  
  ##calculate compositie EWS for single, pairwise, triplicate, etc.combinations
  ##split the data by tube
  
  ##save out results from loops here, using counter to specify the position in the list object finished.results
  finished.results<-NULL
  counter<-0
  
  ##loop for single predictors, double, etc
  ##calculate the number of indicators (ar1, cv, etc)
  no.pred<-min(cum.max,length(unique(kk$variable)))
  ##and their names
  indicators<-sort(unique(kk$variable))
  ##loop through from 1:the number of indicators to be included in the method
  for(h in 1:no.pred){
    #h=2
    ##made a data frame of the indicator names
    inds<-data.frame(indicators)
    ##if this method should include more than 1 indicator in it (i.e. if h>1)
    ## then loop through to the number of indicators to be included, adding a column to "inds" each time
    if(h!=1){
      for(o in 2:h){inds<-cbind(inds, indicators)}
      ##all possible combinations of inds
      combs.all<-as.matrix(expand.grid(inds))
      combs <- combs.all[apply(combs.all,1,function(x) length(unique(x)) )==h,]
      ##make a unique code for each composite indicator
      code <- apply(combs, 1, function(x) paste(sort(x), collapse = " + " ))
      combs <- as.data.frame(combs,stringsAsFactors=F)
      combs$code <- code
      combs <- combs[!duplicated(combs$code),]
      ##and remove any multiple occurances of the same combination of leading indicators
      rownames(combs)<-NULL
    }else{
      combs<-data.frame(indicators= as.character(inds$indicators),
                        code =as.character(inds$indicators),stringsAsFactors = F)
    }
    
    ##then loop through object combs
    result<-FALSE
    for(m in 1:length(combs[,1])){
      #m=20
      ##get results for each leading indicator from the experimental data
      sub.kk<-kk[(kk$variable%in%unlist(combs[m,])),]
      roll<-colSums(do.call("rbind", mclapply(split(sub.kk, list(sub.kk$variable)), function(d){
        if(length(d[,1])>0){
          return(d$value)
        }})))
      
      #make res dtf
      fin.kk <-  data.frame(time)
      ##create an id for the combined predictor
      fin.kk$id<-combs$code[m]
      ##calculate the rolling values
      ##add the rolling value
      fin.kk$value<-roll
      ##calculate the sd and mean of the rolling value
      fin.kk$roll.sd<-rolling_sd(fin.kk$value)
      fin.kk$roll.mean<-rolling_mean(fin.kk$value)
      fin.kk$roll.mean <- ifelse(is.nan(fin.kk$roll.mean),NA,fin.kk$roll.mean)
      fin.kk$roll.sd <- ifelse(is.nan(fin.kk$roll.sd),NA,fin.kk$roll.sd)
      
      ##has a result been generated?
      result<-TRUE
      counter<-counter+1
      
      ##if there is a result add it to the finished data frame, along with info on the number of predictors
      if(result==TRUE){
        fin.kk$no.of.pred<-h
        finished.results[[counter]]<-fin.kk}
    }
  }
  fr2 <- do.call("rbind",finished.results)
  
  fr2$upr2 <- fr2$roll.mean + 2 * fr2$roll.sd
  fr2$lwr2 <- fr2$roll.mean - 2 * fr2$roll.sd
  fr2$upr1 <- fr2$roll.mean + 1 * fr2$roll.sd
  fr2$lwr1 <- fr2$roll.mean - 1 * fr2$roll.sd
  return(fr2)
}

# find significant EWS signal and their date
find.ews <- function(Composite, sd=2,CrashDate=0){
  require(dplyr)
  
  earliestDate <- function(date,includeTF){
    min(date[includeTF==T],na.rm=T)
  }
  
  # stop if invalid sd value
  if(!is.numeric(sd) | !(length(sd) %in% c(1,nrow(Composite))) ) {
    stop("there is a problem with the sd")
  }
  # calculate limits given a certain limit
  Composite$myupr <- Composite$roll.mean + sd*Composite$roll.sd
  Composite$mylwr <- Composite$roll.mean - sd*Composite$roll.sd
  
  Composite$time<-Composite$time-CrashDate
  Res <- Composite %>% group_by(id) %>%
    mutate(TrueSignal=value>myupr,
           EWS_detect=sum(TrueSignal,na.rm=T)>0,
           EWS_Date=ifelse(EWS_detect==T,earliestDate(time,TrueSignal),NA)
    )
  
  return(Res)
  
}

# summarise Signals (is a signal detected and when)
summary.ews <- function(Signals){
  Detections <- Signals %>% group_by(id) %>% 
    summarise(no.of.pred=mean(no.of.pred),
              EWS_detect=sum(TrueSignal,na.rm=T)>0,
              EWS_Date=mean(EWS_Date)) %>% 
    arrange(EWS_Date)
  message(paste("detected" ,
                sum(Detections$EWS_detect),
                "signals out of",
                nrow(Detections),
                "tested\n", sum(Detections$EWS_Date<0,na.rm=T),"are predictive"
  ))
  
  return(Detections)
} 
