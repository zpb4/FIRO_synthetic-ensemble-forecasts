#Converts syn-GEFS HEFS output files to R data arrays for processing
#Require
library(stringr)
library(abind)

#reads raw data and also outputs processed data to 'data' repository
setwd('FIRO_synthetic-ensemble-forecasts/data')

ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day') #same length
ixx3<-as.POSIXlt(ix3)
m<-length(ix3)
num_ens<-10

locs<-c('LAMC1','UKAC1','HOPC1','HOPC1L')

#for(l in 1:length(locs)){  Note: current dataset only has LAMC location for plot comparisons, but other sites are available
for(l in 1:1){
  for(s in 1:num_ens){
    loc<-locs[l]
    folder<-paste('syn-gefs_ensembles/ens',s,'/',sep='')
    #folder<-paste('syn-gefs_ensembles_',loc,'/ens',s,'/',sep='')
  
    #1) Process raw .csv data
    i<-1
    fname<-paste(folder,ixx3[i]$year+1900,str_pad(ixx3[i]$mon+1,2,'left','0'),str_pad(ixx3[i]$mday,2,'left','0'),'12_',loc,'_SQIN_hourly.csv',sep='')
    dat1 <- read.csv(fname)
    dat1[2:337,2:62]<-NA
  
    hefs_raw<-array(NA,c(m,336,61))
  
    for(i in 1:m){
      fname<-paste(folder,ixx3[i]$year+1900,str_pad(ixx3[i]$mon+1,2,'left','0'),str_pad(ixx3[i]$mday,2,'left','0'),'12_',loc,'_SQIN_hourly.csv',sep='')
      dat <- tryCatch(read.csv(fname),error = function (e){print(paste(fname,'no data'));return(dat1)},warning = function (w){print(paste(fname,'no data'));return(dat1)})
      lst_dat<-dat[2:337,2:62]
      for(j in 1:61){
        out<-lst_dat[[j]]
        #out<-levels(lst_dat[,j])[lst_dat[,j]]
        out2<-as.numeric(out)
        hefs_raw[i,,j]<-as.numeric(out2)
      }
    }
  
  saveRDS(hefs_raw,paste(folder,'hefs_raw_syn-gefs.rds',sep=''))
  
  #2)Calculate a daily mean from hourly values synched at 12Z
  idx<-cbind(seq(1,336,24),seq(24,336,24))
  
  hefs_dly_mean<-array(NA,c(m,14,61))
  
  for(i in 1:14){
    mn<-apply(hefs_raw[,idx[i,1]:idx[i,2],],c(1,3),function(x){mean(x,na.rm=T)})
    hefs_dly_mean[,i,]<-mn
  }
  
  saveRDS(hefs_dly_mean,paste(folder,'hefs_dly_mean_syn-gefs_12z.rds',sep=''))
  
  #3) Rearrange arrays to synch forecasts with observations days
  hefs_ens_forc<-array(0,c(61,m,14))
  
  for(i in 1:14){
    int<-array(0,c(61,m))
    int[,(i+1):m]<-t(hefs_dly_mean[,i,])[,-c((m-i+1):m)]
    hefs_ens_forc[,,i]<-int
  }
  
  na_fun<-function(x){for(i in 1:length(x)){
    if(is.na(x[i])==T){x[i]<-x[i-1]+((x[i+1]-x[i-1])/2)}}
    return(x)}
  
  for(i in 1:61){
    hefs_ens_forc[i,,]<-apply(hefs_ens_forc[i,,],2,na_fun)
  }
  
  saveRDS(hefs_ens_forc,paste(folder,'hefs_ens_forc_syn-gefs_12z.rds',sep=''))
  
  #4) Calculate ensemble mean values and save
  hefs_ens_mean_forc<-apply(hefs_ens_forc,c(2,3),mean)
  
  saveRDS(hefs_ens_mean_forc,paste(folder,'hefs_ens_mean_forc_syn-meteo_12z.rds',sep=''))
  
  }
}

rm(list=ls());gc()

##################################################END###################################################
