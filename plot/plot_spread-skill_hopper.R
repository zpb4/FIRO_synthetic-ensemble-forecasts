
plt_sprd_skill<-function(nsamps,nsets,lds,lpcnt,upcnt,show_leads,show_xlab,plt_par,outfile,plt_lim){

#nsamps<-10
#nsets<-10
#lds<-c(1,3,5,10)
#lpcnt<-0
#upcnt<-100
#show_leads<-F
#show_xlab<-F
#plt_par<-T
#outfile<-'v10_wtsamp30'
#plt_lim<-1500

  
setwd('h:/firo_lamc/')
library(scales)
library(abind)

leads<-14
ens_num<-61
pct<-c(lpcnt,upcnt)*0.01
pct[pct==0]<-0.001

source('plot/forecast_verification_functions.R')
inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]

ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')
ix2<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ixx2<-as.POSIXlt(ix2)

hefs_forc_lamc<-readRDS('data/hefs_lamc_act-meteo/hefs_ens_forc_act-meteo_12z.rds')

ens_set<-sample(1:nsets,1)
syn_hefs_lamc<-readRDS(paste('z:/syn_forecast_test/output/',outfile,'/syn_hefs_flow_8510_',outfile,'_ens_',ens_set,'_12z.rds',sep=''))

syn_gefs_lamc<-array(NA,c(ens_num,length(ix3),leads,nsamps))

for(i in 1:nsamps){
  syn_gefs<-readRDS(paste('data/hefs_lamc_syn-meteo_v3/ens',i,'/hefs_ens_forc_syn-meteo_12z.rds',sep=''))
  syn_gefs_lamc[,,,i]<-syn_gefs
}


hefs_lam<-hefs_forc_lamc[,hefs_idx,]*1000
sgefs_lam<-syn_gefs_lamc[,hefs_idx,,]*1000
shefs_lam<-syn_hefs_lamc[15:28,,,]

#cold season index with quantile split
dec_mar<-sort(c(which(ixx2$mon>10),which(ixx2$mon<3)))
ob_ss<-sort(obs$LAMC1[dec_mar],index.return=T)
cut_ix<-(round(length(ob_ss$x)*pct[1])):(round(length(ob_ss$x)*pct[2]))
ss_ixx<-sort(ob_ss$ix[cut_ix])
ss_ix<-dec_mar[ss_ixx]
ss_idx<-ss_ix


hefs_lamc<-hefs_lam[,ss_idx,]
shefs_lamc<-shefs_lam[,ss_idx,,]
sgefs_lamc<-sgefs_lam[,ss_idx,,]
obs_lamc<-obs$LAMC1[ss_idx]

####################calculate##################################

var_hefs<-apply(hefs_lamc,c(2,3),st2_ens)
mn_hefs<-apply(hefs_lamc,c(2,3),mn_ens)

var_shefs<-array(NA,c(nsamps,leads,length(ss_idx)))
var_sgefs<-array(NA,c(nsamps,leads,length(ss_idx)))

mn_shefs<-array(NA,c(nsamps,leads,length(ss_idx)))
mn_sgefs<-array(NA,c(nsamps,leads,length(ss_idx)))

for(n in 1:nsamps){
  for(j in 1:leads){
    for(k in 1:length(ss_idx)){
      var_shefs[n,j,k]<-st2_ens(shefs_lamc[j,k,,n])
      mn_shefs[n,j,k]<-mn_ens(shefs_lamc[j,k,,n])

      var_sgefs[n,j,k]<-st2_ens(sgefs_lamc[,k,j,n])
      mn_sgefs[n,j,k]<-mn_ens(sgefs_lamc[,k,j,n])
    }
  }
}

bin_spd_mse_hefs<-array(NA,c(leads,2,20))
bin_spd_mse_shefs<-array(NA,c(nsamps,leads,2,20))
bin_spd_mse_sgefs<-array(NA,c(nsamps,leads,2,20))

for(k in 1:leads){
  ct_hefs<-cut(rank(var_hefs[,k],ties.method = 'first'),20,labels=1:20)
  for(j in 1:20){
    hefs_var_sset<-var_hefs[which(ct_hefs==j),k]
    hefs_mn_sset<-mn_hefs[which(ct_hefs==j),k]
    bin_spd_mse_hefs[k,1,j]<-bin_spread(hefs_var_sset,ens_num)
    bin_spd_mse_hefs[k,2,j]<-bin_mse(hefs_mn_sset,obs_lamc[which(ct_hefs==j)],ens_num)
  }
}

for(n in 1:nsamps){
  for(k in 1:leads){
    ct_shefs<-cut(rank(var_shefs[n,k,],ties.method = 'first'),20,labels=1:20)
    ct_sgefs<-cut(rank(var_sgefs[n,k,],ties.method = 'first'),20,labels=1:20)
    for(j in 1:20){
      shefs_var_sset<-var_shefs[n,k,which(ct_shefs==j)]
      shefs_mn_sset<-mn_shefs[n,k,which(ct_shefs==j)]
      sgefs_var_sset<-var_sgefs[n,k,which(ct_sgefs==j)]
      sgefs_mn_sset<-mn_sgefs[n,k,which(ct_sgefs==j)]
      bin_spd_mse_shefs[n,k,1,j]<-bin_spread(shefs_var_sset,ens_num)
      bin_spd_mse_shefs[n,k,2,j]<-bin_mse(shefs_mn_sset,obs_lamc[which(ct_shefs==j)],ens_num)
      bin_spd_mse_sgefs[n,k,1,j]<-bin_spread(sgefs_var_sset,ens_num)
      bin_spd_mse_sgefs[n,k,2,j]<-bin_mse(sgefs_mn_sset,obs_lamc[which(ct_sgefs==j)],ens_num)
    }
  }
}

####################plot##################################

if(plt_par==T){
  par(mfrow=c(1,length(lds)),mar=c(3,4,0,1),mgp=c(2,0.5,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3.5,3.5,0.5))
}
ylb<-c('Ensemble Mean Error',rep("",length(lds)-1))

alp_orange<-alpha('orange',alpha = 0.15)
alp_purp<-alpha('purple',alpha = 0.1)

for(i in 1:length(lds)){
  plot(bin_spd_mse_hefs[lds[i],1,],bin_spd_mse_hefs[lds[i],2,],type='b',main='',
       xlab='Ensemble Spread',ylab=ylb[i],xlim=c(0,plt_lim),ylim=c(0,plt_lim),lwd=2,pch=15,axes=F)
  axis(1,at=seq(0,1500,500),labels=seq(0,1500,500))
  axis(2,at=seq(0,1500,500),labels=F)
  if(i==1){
    axis(2,at=seq(0,1500,500),labels = seq(0,1500,500))
  }
  for(n in 1:nsamps){
    lines(bin_spd_mse_sgefs[n,i,1,],bin_spd_mse_shefs[n,i,2,],col=alp_orange,lwd=2)
    lines(bin_spd_mse_sgefs[n,i,1,],bin_spd_mse_sgefs[n,i,2,],col=alp_purp,lwd=2)
  }
  #lines(bin_spd_mse_sgefs[1,i,1,],bin_spd_mse_shefs[1,i,2,],col='orange',lwd=1,type='b',pch=19)
  #lines(bin_spd_mse_sgefs[1,i,1,],bin_spd_mse_sgefs[1,i,2,],col='purple',lwd=1,type='b',pch=18)
  lines(apply(bin_spd_mse_sgefs[,i,1,],2,median),apply(bin_spd_mse_shefs[,i,2,],2,median),col='orange',lwd=1,type='b',pch=19)
  lines(apply(bin_spd_mse_sgefs[,i,1,],2,median),apply(bin_spd_mse_sgefs[,i,2,],2,median),col='purple',lwd=1,type='b',pch=18)
  lines(bin_spd_mse_hefs[lds[i],1,],bin_spd_mse_hefs[lds[i],2,],type='b',lwd=2,pch=15)
  if(i==1){
    legend('topleft',c('HEFS','syn-HEFS','syn-GEFS'),col=c('black','orange','purple'),lwd=c(2,1,1),pch=c(15,19,18),cex=1.75)}
  abline(0,1,col='gray',lwd=2,lty=2)
  if(show_leads==T){
    mtext(paste('Lead',lds[i]),side=3,font=2,line=1,cex=2.5)}
  box(which='plot')
}

}

########################################END################################################