# load data and packages
load("cache/cleanData_20170127.Rdata")
library(tidyverse)
library(reshape2)
source("R/myfunc1.R")

# get single indices ---------------------------------------------------------------------
# remove unwanted things
dat<-df.pheno %>% select(-meanAFR,-sdAFR) %>% rename(Count=fem)

# get population growth smoother
ratio<-dat$Count[2:length(dat$Count)]/dat$Count[1:(length(dat$Count)-1)]
day<-dat$yr[1:length(dat$yr)-1]
lo.1<-loess(ratio~day)
x.seq<-seq(min(dat$yr),max(dat$yr))
pred.lo.1<-predict(lo.1, x.seq, se=TRUE)
lambs<-na.omit(data.frame(x.seq, pred.lo.1[1:2], row.names=NULL))
colnames(lambs) <- c("yr","smoother.pg","se.pg")
dat <- merge(dat,lambs)

#blank objects to save results
RES<-NULL
roll.cv<-NULL
roll.acf<-NULL
roll.ratio<-NULL
roll.ar<-NULL
roll.yrlingSize<-NULL
roll.return.rate<-NULL
roll.density.ratio<-NULL
roll.meanWt114AgeCorF <- NULL
roll.sdWt114AgeCorF <- NULL
roll.He <- NULL
roll.Ar <- NULL
roll.Fis <- NULL


##looped to calculate the rolling change
for(i in 2:length(dat$yr)){
  ##subset the population of interest up until dat i
  dat.t<-subset(dat, yr<=unique(sort(dat$yr))[i])
  
  ##calculate the CV at time t in the focal pop, relative to mean CV through time of that pop
  roll.cv[[i]]<-CV(dat.t$Count, TRUE)
  CV.t<-(CV(dat.t$Count, TRUE)-(mean(roll.cv, na.rm=T)))/sd(roll.cv, TRUE)
  
  ## for autocorrelation
  roll.acf[[i]]<-acf(dat.t$Count, lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]
  acf.t<-(acf(dat.t$Count, lag.max = 1, type = c("correlation"), plot = FALSE)$acf[2]-mean(roll.acf, na.rm=T))/sd(roll.acf, TRUE)
  
  ## for ar1
  if(length(which(diff(dat.t$Count)!=0))>0){ # minimum data  required to fitauto-regressive model
    roll.ar[[i]]<-ar.ols(dat.t$Count, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]
    ar.t<-(ar.ols(dat.t$Count, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar[1]-mean(roll.ar, na.rm=T))/sd(roll.ar, TRUE)
  }else{ar.t<-NA}
  
  ## for lamb:ewe ratio
  roll.ratio[[i]]<-mean(dat.t$lb.ewe.ratio, na.rm=T)
  ratio.t<- (mean(dat.t$lb.ewe.ratio, na.rm=T)-mean(roll.ratio, na.rm=T))/sd(roll.ratio, TRUE)
  
  ##for yrlingSize
  roll.yrlingSize[[i]]<-mean(dat.t$yrlingWt, na.rm=T)
  yrlingSize.t<-(mean(dat.t$yrlingWt, na.rm=T)-mean(roll.yrlingSize, na.rm=T))/sd(roll.yrlingSize, TRUE)
  
  ## for mean body size
  roll.meanWt114AgeCorF[[i]]<-mean(dat.t$meanWt114AgeCorF, na.rm=T)
  AgeCorsize.t<- (mean(dat.t$meanWt114AgeCorF, na.rm=T)-mean(roll.meanWt114AgeCorF, na.rm=T))/sd(roll.meanWt114AgeCorF, TRUE)
  
  ##for sd body size
  roll.sdWt114AgeCorF[[i]]<-mean(dat.t$sdWt114AgeCorF, na.rm=T)
  sdWt114AgeCorF.t<-(mean(dat.t$sdWt114AgeCorF, na.rm=T)-mean(roll.sdWt114AgeCorF, na.rm=T))/sd(roll.sdWt114AgeCorF, TRUE)
  
  ## for he
  roll.He[[i]]<-mean(dat.t$He, na.rm=T)
  He.t<- (mean(dat.t$He, na.rm=T)-mean(roll.He, na.rm=T))/sd(roll.He, TRUE)
  
  ## for Ar
  roll.Ar[[i]]<-mean(dat.t$Ar, na.rm=T)
  Ar.t<- (mean(dat.t$Ar, na.rm=T)-mean(roll.Ar, na.rm=T))/sd(roll.Ar, TRUE)
  
  ## for Fis
  roll.Fis[[i]]<-mean(dat.t$Fis, na.rm=T)
  Fis.t<- (mean(dat.t$Fis, na.rm=T)-mean(roll.Fis, na.rm=T))/sd(roll.Fis, TRUE)
  
  ##for density ratio
  if(length(which(diff(dat.t$Count)>0))>0){ # minimum data  required to fit density function
    spectfft <- spec.ar(dat.t$Count, n.freq = length(dat.t$Count), plot = FALSE, order = 1)
    roll.density.ratio[[	i]] <- spectfft$spec[1]/spectfft$spec[length(dat.t$Count)]
    den.t<-(roll.density.ratio[[i]]-(mean(roll.density.ratio, na.rm=T)))/sd(roll.density.ratio, TRUE)
  }else{den.t<-NA}
  
  ##for return rate
  roll.return.rate[[i]]<-(1/ar.t)
  ret.t<-(((1/ar.t)-(mean(roll.return.rate, na.rm=T)))/sd(roll.return.rate, TRUE))
  
  ##save results
  RES[[i]]<-data.frame(dat.t[i,], "cv"=CV.t, "acf"=acf.t, "ar1"=ar.t, "dr"=den.t,"rr"=ret.t, 
                       "ratio"=ratio.t, "yrlingSize"=yrlingSize.t,
                       "ageCorWt"=AgeCorsize.t,"ageCorWt.sd"=sdWt114AgeCorF.t,
                       "He.i"=He.t,"Ar.i"=Ar.t,"Fis.i"=Fis.t)
}
rm<-do.call("rbind", RES)
rm$norm.day<-(rm$yr-max(rm$yr))-1

# flip negative indicator so that bigger is bad
rm$rr<-rm$rr*-1
rm$ageCorWt<-rm$ageCorWt*-1
rm$yrlingSize <- rm$yrlingSize * -1
rm$ratio <- rm$ratio * -1
rm$He.i <- rm$He.i*-1
rm$Ar.i <- rm$Ar.i*-1

#remove nan
rm$He.i[is.nan(rm$He.i)] <- NA
rm$Ar.i[is.nan(rm$Ar.i)] <- NA
rm$Fis.i[is.nan(rm$Fis.i)] <- NA


# make graph of indicators
lab=data.frame(y=as.numeric(t(rm[rm$yr==1995,c(24:35)])))
lab$var=colnames(rm)[24:35]

ggplot(rm,aes(x=yr))+
  geom_path(aes(y=pg),color="black")+
  geom_path(aes(y=smoother.pg),color="black",linetype="dashed")+
  geom_path(aes(y=cv,color="cv"))+
  geom_path(aes(y=acf,color="acf"))+
  geom_path(aes(y=ar1,color="ar1"))+
  geom_path(aes(y=dr,color="density Ratio"))+
  geom_path(aes(y=rr,color="return-rate"))+
  geom_path(aes(y=ratio,color="ratio"))+
  geom_path(aes(y=ageCorWt,color="ageCorWt"))+
  geom_path(aes(y=ageCorWt.sd,color="ageCorWt.sd"))+
  geom_path(aes(y=yrlingSize,color="yrlingSize"))+
  geom_path(aes(y=He.i,color="He"))+
  geom_path(aes(y=Ar.i,color="Ar"))+
  geom_path(aes(y=Fis.i,color="Fis"))+
  geom_hline(yintercept = 1,color="darkgrey")+
  geom_vline(xintercept = 1992,color="darkgrey")+
  xlim(1976,1995)+  
  labs(x="Year", y="EWS value & limits",title="ram mountian")+theme_gray(16)+
  geom_label(data=lab,aes(x=1995,y=y,label=var),hjust=0)


# make composite and fins significant signals -----------------------------

rm2 <- makeComposite(rm$yr,rm[,c(24:35)],cum.max = 4)
rm.Signals <- find.ews(rm2,sd = 2,CrashDate = 1992)
rm.ewsSummary <- summary.ews(rm.Signals)

find.ews(rm2,sd = qnorm(1-(0.05/12)) ,CrashDate = 1992) %>% summary.ews()

#  plot single-double signals with limits
ggplot(rm2[rm2$no.of.pred<=2,],aes(x =time,y=value))+
  geom_path()+facet_wrap(~id)+
  geom_path(aes(y=lwr1),linetype="dotted")+geom_path(aes(y=upr1),linetype="dotted")+
  geom_path(aes(y=lwr2),linetype="dashed")+geom_path(aes(y=upr2),linetype="dashed")+
  labs(x="Year", y="EWS value & limits",title="Ram Mountain")+theme_gray(16)+
  geom_vline(xintercept = 1992,color="darkgrey")

# plot significant signals
ggplot(rm2[rm2$id %in% rm.ewsSummary[rm.ewsSummary$EWS_detect==T,]$id,],
       aes(x =time,y=value))+
  geom_path()+facet_wrap(~id)+
  geom_path(aes(y=lwr1),linetype="dotted")+geom_path(aes(y=upr1),linetype="dotted")+
  geom_path(aes(y=lwr2),linetype="dashed")+geom_path(aes(y=upr2),linetype="dashed")+
  labs(x="Year", y="EWS value & limits",title="Ram Mountain")+theme_gray(16)+
  geom_vline(xintercept = 1992,color="darkgrey")

# plot genetic signals
ggplot(rm2[rm2$id %in% c("He.i","Ar.i","Fis.i"),], aes(x =time,y=value))+
  geom_path()+facet_wrap(~id)+
  geom_path(aes(y=lwr1),linetype="dotted")+geom_path(aes(y=upr1),linetype="dotted")+
  geom_path(aes(y=lwr2),linetype="dashed")+geom_path(aes(y=upr2),linetype="dashed")+
  labs(x="Year", y="EWS value & limits",title="Genetic EWS at Ram Mountain")+
  theme_gray(16)+
  geom_vline(xintercept = 1992,color="darkgrey")

