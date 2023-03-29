
plt_cumul_ensembles<-function(nsamps,nsets,lds,site,forc,show_leads,show_xlab,show_ylab,show_yaxslab,show_legend,plt_par,plot_samp,
                        plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                        sgefs_samp_seed,shefs_samp_seed,event,outfile){

#nsets<-10
#nsamps<-10
#lds<-c(6)
#site<-'lamc'
#forc<-'hefs'
#show_leads<-T
#show_xlab<-T
#show_ylab<-T
#show_yaxslab
#show_legend<-T
#plot_samp<-T
#plot_median<-F
#plot_mean<-T
#bounds<-c(0.05,0.95)
#sgefs_seed<-1
#shefs_set_seed<-1
#shefs_seed<-1
#hefs_samp_seed<-1
#sgefs_samp_seed<-1
#shefs_samp_seed<-1
#event<-'1986-02-18'
#outfile<-'v9_wtsamp30'

setwd('h:/firo_lamc/')
library(scales)
library(abind)

leads<-14
ens_num<-61

cfs_taf<-2.29568411*10**-5 * 86400 / 1000
roll_sum<-function(x){out<-c();out[1]<-x[1];
for(i in 1:(length(x)-1)){out[i+1]<-out[i]+x[i+1]};return(out)}
ubnds<-function(x){y<-sort(x);return(y[round(bounds[2]*ens_num)])}
lbnds<-function(x){y<-sort(x);return(y[round(bounds[1]*ens_num)])}

source('plot/forecast_verification_functions.R')
inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]

ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

hefs_forc<-readRDS(paste('data/hefs_',site,'_act-meteo/hefs_ens_forc_act-meteo_12z.rds',sep=''))

set.seed(shefs_set_seed);ens_set<-sample(1:nsets,1)

syn_hefs<-readRDS(paste('z:/syn_forecast_test/output/',outfile,'/syn_hefs_flow_8510_',outfile,'_ens_',ens_set,'_12z.rds',sep=''))

set.seed(shefs_seed);ens_shefs<-sample(1:nsamps,1)
syn_hefs<-syn_hefs[,,,ens_shefs]

set.seed(sgefs_seed);ens_sgefs<-sample(1:nsamps,1)
syn_gefs<-readRDS(paste('data/hefs_',site,'_syn-meteo_v3/ens',ens_sgefs,'/hefs_ens_forc_syn-meteo_12z.rds',sep=''))

locs<-c('hopc','lamc','ukac')
site_idx<-which(locs==site)
idx_locs<-list(1:14,15:28,29:42)
shefs_idx<-idx_locs[[which(locs==site)]]

hefs<-hefs_forc[,hefs_idx,]*1000
sgefs<-syn_gefs[,hefs_idx,]*1000
shefs<-syn_hefs[shefs_idx,,]

hefs_out<-array(NA,c(dim(shefs)[2],(dim(shefs)[1]+1),ens_num))
shefs_out<-array(NA,c(dim(shefs)[2],(dim(shefs)[1]+1),ens_num))
sgefs_out<-array(NA,c(dim(shefs)[2],(dim(shefs)[1]+1),ens_num))

hefs_out[,1,]<-matrix(rep(obs[,site_idx],ens_num),ncol=ens_num,byrow=F)
shefs_out[,1,]<-matrix(rep(obs[,site_idx],ens_num),ncol=ens_num,byrow=F)
sgefs_out[,1,]<-matrix(rep(obs[,site_idx],ens_num),ncol=ens_num,byrow=F)

for(i in 1:14){
  hefs_out[,(i+1),]<-t(hefs[1:ens_num,c((i+1):dim(shefs)[2],rep(dim(shefs)[2],(i))),i])
  shefs_out[,(i+1),]<-shefs[i,c((i+1):dim(shefs)[2],rep(dim(shefs)[2],(i))),1:ens_num]
  sgefs_out[,(i+1),]<-t(sgefs[1:ens_num,c((i+1):dim(shefs)[2],rep(dim(shefs)[2],(i))),i])
}

#Ensemble Plots--------------------------------------------------------------
#combine plot
#4 panel plot
#plot HEFS
#par(mfrow=c(1,1),mar=c(2,3,0,1),mgp=c(2,0.75,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(2.5,2.5,2,0.5))
if(forc=='hefs'){
  synf<-hefs_out
  clr<-'blue'
  labs<-c('HEFS','HEFS mean','HEFS median')
  flab<-'HEFS'
  set.seed(hefs_samp_seed);synf_samp<-sample(1:ens_num,1)
}

if(forc=='shefs'){
  synf<-shefs_out
  clr<-'dark orange'
  labs<-c('sHEFS','sHEFS mean','sHEFS median')
  flab<-'syn-HEFS'
  set.seed(shefs_samp_seed);synf_samp<-sample(1:ens_num,1)
}

if(forc=='sgefs'){
  synf<-sgefs_out
  clr<-'purple'
  labs<-c('sGEFS','sGEFS mean','sGEFS median')
  flab<-'syn-GEFS'
  set.seed(sgefs_samp_seed);synf_samp<-sample(1:ens_num,1)
}

idx<-which(ix==event)
  
ylm<-1.5*max(roll_sum(synf[(idx-i):(idx+14-i),1,1]) * cfs_taf)

for(i in lds){
  plot(0:14,roll_sum(synf[(idx-i):(idx+14-i),1,1]) * cfs_taf,type='l',lwd=3,ylim=c(0,ylm),
       xlab='',ylab='Cumulative Flow (TAF)',main='',axes=F)
  axis(2,at=seq(0,ylm,20),labels=F)
  axis(1,at=0:14,labels = F)
  box(which='plot')
  for(j in 1:ens_num){
    lines(0:14,roll_sum(synf[idx-i,,j]) * cfs_taf,col='light gray')
  }
  cumul_flow<-apply(synf[idx-i,,],2,roll_sum)
  lines(0:14,apply(cumul_flow,1,ubnds) * cfs_taf,lwd=3,lty=3,col='dark gray')
  lines(0:14,apply(cumul_flow,1,lbnds) * cfs_taf,lwd=3,lty=3,col='dark gray')
  if(plot_mean==T){
    lines(0:14,apply(cumul_flow,1,mean) * cfs_taf,col=clr,lwd=3,lty=2)
    if(show_legend==T){
      legend(leg_pos,c('Obs',labs[1],labs[2]),
           col=c('black','light gray',clr,'dark gray'),lwd=c(3,1,3,3),lty=c(1,1,2,3),cex=2)}
  }
  if(plot_median==T){
    lines(0:14,apply(cumul_flow,1,median) * cfs_taf,col=clr,lwd=3,lty=3)
    if(show_legend==T){
      legend(leg_pos,c('Obs',labs[1],labs[3]),
          col=c('black','light gray',clr,'dark gray'),lwd=c(3,1,3,3),lty=c(1,1,2,3),cex=2)}
  }
  if(plot_samp==T){
    lines(0:14,roll_sum(synf[(idx-i),,synf_samp]) * cfs_taf,col='magenta')
  }
  lines(0:14,roll_sum(synf[(idx-i):(idx+14-i),1,1]) * cfs_taf,lwd=3,col='black')
  abline(v=i,col='red',lty=3)
  abline(v=0,col='black',lty=3)
  text(i,0.9*ylm,ix[idx],adj=0,col='red',cex=2)
  text(0,0.9*ylm,ix[idx-i],adj=0,col='black',cex=2)
  if(show_ylab==T){mtext(flab,side=2,line=5,font=2,cex=2.5)}
  if(show_leads==T){mtext(paste('Lead',i),side=3,line=1,font=2,cex=2.5)}
  if(show_xlab==T){
    axis(1,at=0:14,labels = 0:14)
    mtext('Forecast Day',side=1,line=3,font=2,cex=1.5)}
  if(show_yaxslab==T){
    axis(2,at=seq(0,ylm,20),seq(0,ylm,20))}
}

}


########################################END################################################