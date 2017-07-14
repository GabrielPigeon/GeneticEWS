# © Gabriel Pigeon
# Version 1          date: 2017-07-14
# Script for Data management and cleaning
# goes from Raw data to cleaned data
# saved cleaned data in "cache"

source("R/myfunc1.R")
library(gabtool)
#############################################################################.

# load data ---------------------------------------------------------------
df.RamMtnpop <- readxl::read_excel("data/ram/RamMtnpop.xls")
df.wt <- read_csv("data/ram/adjwt1972-2016_noMaxDeltaJ.csv")
df.surv <-  read_csv("data/ram/surv.csv")
RamGen <- read_csv("data/ram/genetics.csv")

# add afr -----------------------------------------------------------------
ageOfRepro <- function(reproduced,age){
    tmp=na.omit(reproduced*age)
    tmp=tmp[!tmp==0]
    return(tmp)
  }

df.surv <- df.surv %>% mutate(reproduced=code.sr %in% c(1:5)) %>% group_by(ID) %>% mutate(afr=min(ageOfRepro(reproduced,age)))
df.surv[df.surv$cohort<=1972,"afr"] <- NA
df.surv[df.surv$sex=="male","afr"] <- NA
df.surv[is.infinite(df.surv$afr),"afr"] <- NA
df.surv[df.surv$cohort<=1972,"afr"] <- NA

df.surv$ageReprod <- df.surv$age
df.surv[df.surv$reproduced==F,"ageReprod"] <- NA
df.surv[df.surv$sex=="male","ageReprod"] <- NA
# get yrly pheno ---------------------------------------------------
# *** mean age ,mean age of reproducer , mean afr ,sd afr   ---------

surv.yrly <- df.surv %>% group_by(yr) %>% 
  summarise(meanAge=mean(age),sdAge=sd(age),
            meanAgeReprod=mean(ageReprod,na.rm=T), sdAgeReprod=sd(ageReprod,na.rm=T),
            meanAFR=mean(afr,na.rm=T),sdAFR=sd(afr,na.rm=T)  )
surv.yrly[surv.yrly$yr<max(df.surv[df.surv$cohort==1972,"yr"]),c("meanAFR","sdAFR")] <- NA
surv.yrly[surv.yrly$yr>=2006,c("meanAFR","sdAFR")] <- NA


# *** juv to adult ratio   ---------------------------------------
ratio <- df.surv %>% group_by(yr) %>% 
  summarise(nb.lb=sum(age==0),nb.ewe=sum((age>=2 & sex=="female")),lb.ewe.ratio=nb.lb/nb.ewe) %>% 
  select(yr,lb.ewe.ratio)

# *** mean wt ,sd wt  ---------------------------------------------
maxDelaJJ <- 31
df.wt <- df.wt %>% mutate(yr=year) %>% select(-year)
df.wt <- df.surv %>% select(ID,sex,yr,age) %>% right_join(df.wt)

df.wt12 <- df.wt %>%  filter(JJ<25) %>% mutate(wt12=adjwt.mixed)
df.wt12[is.na(df.wt12$wt12),"wt12"] <- df.wt12[is.na(df.wt12$wt12),"adjwt.pop"]
df.wt12[is.na(df.wt12$wt12),"wt12"] <- df.wt12[is.na(df.wt12$wt12),"adjwt.actual"]
df.wt12[df.wt12$DeltaJJ>maxDelaJJ,"wt12"] <- NA



df.wt114 <- df.wt %>%  filter(JJ>25) %>% mutate(wt114=adjwt.mixed)
df.wt114[is.na(df.wt114$wt114),"wt114"] <- df.wt114[is.na(df.wt114$wt114),"adjwt.pop"]
df.wt114[is.na(df.wt114$wt114),"wt114"] <- df.wt114[is.na(df.wt114$wt114),"adjwt.actual"]
df.wt114[df.wt114$DeltaJJ>maxDelaJJ,"wt114"] <- NA

df.wtB <- full_join(df.wt12[,-(5:14)],df.wt114[,-(5:14)])
df.wtB <- df.wtB %>% group_by(age,sex) %>% 
  mutate(wt12AgeCor=wt12 - mean(wt12,na.rm=T), wt114AgeCor=wt114 - mean(wt114,na.rm=T))

# wt.yrly <- df.wtB %>% group_by(yr) %>% 
#   summarise(meanWt12=mean(wt12,na.rm=T),sdWt12=sd(wt12,na.rm=T),
#             meanWt114=mean(wt114,na.rm=T),sdWt114=sd(wt114,na.rm=T),
#             meanWt12AgeCor=mean(wt12AgeCor,na.rm=T),sdWt12AgeCor=sd(wt12AgeCor,na.rm=T),
#             meanWt114AgeCor=mean(wt114AgeCor,na.rm=T),sdWt114AgeCor=sd(wt114AgeCor,na.rm=T) )
wt.yrlyF <- df.wtB %>% filter(sex=="female") %>% group_by(yr) %>% 
  summarise(meanWt12F=mean(wt12,na.rm=T),sdWt12F=sd(wt12,na.rm=T),
            meanWt114F=mean(wt114,na.rm=T),sdWt114F=sd(wt114,na.rm=T),
            meanWt12AgeCorF=mean(wt12AgeCor,na.rm=T),sdWt12AgeCorF=sd(wt12AgeCor,na.rm=T),
            meanWt114AgeCorF=mean(wt114AgeCor,na.rm=T),sdWt114AgeCorF=sd(wt114AgeCor,na.rm=T) )

# *** mean + sd yearling wt   ---------------------------------------

yrlingwt <- df.wt114 %>% filter(age==1) %>% group_by(yr) %>% 
  summarise(yrlingWt=mean(wt114,na.rm=T))


# output yrly file --------------------------------------------------------
df.pheno <- df.RamMtnpop %>% select(yr,fem,total) 
df.pheno$pg <- c(exp(diff(log(df.pheno$total))),NA)
df.pheno <- df.pheno %>% full_join(surv.yrly) %>% #full_join(wt.yrly) %>% 
  full_join(wt.yrlyF) %>% left_join(ratio) %>% left_join(yrlingwt)

df.pheno <- df.pheno %>% mutate_if(is.numeric,function(x){ifelse(is.nan(x),yes = NA,no = x)})
df.pheno <- df.pheno %>% filter(yr>=1975 & yr<=2015)

#merge genetic
df.pheno <- left_join(df.pheno,RamGen)

# get pop growth and collapse date ----------------------------------------





# output clean data -------------------------------------------------------
# remove useless object
rm(df.wtB,df.wt12,df.wt114,df.wt,ageOfRepro,wt.yrlyF,RamGen)


# save cleaned data, ready for analysis
save(df.pheno,df.surv,df.RamMtnpop,maxDelaJJ,file = "cache/cleanData_20170127.Rdata")
