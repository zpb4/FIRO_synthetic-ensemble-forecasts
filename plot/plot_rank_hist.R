
plt_rank_hist<-function(m,lds,lpcnt,upcnt,show_leads,syn_plot,sgefs_typ,plt_par,show_xlab){
setwd('h:/firo_lamc/')
library(scales)

inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]
#m<-1
leads<-14
ens_num<-61
#lds<-c(1,3,5,10)
#lpcnt<-90
#upcnt<-100
pct<-c(lpcnt,upcnt)*0.01
pct[pct==0]<-0.001
#show_leads<-T
#show_xlab<-T
#syn_plot<-'sGEFS'
#plt_par<-'off'


#shefs_typ<-'loessv3_md_prinorm-lm-glob-mn-corr'

source('plot/forecast_verification_functions.R')

ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')
ix2<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ixx2<-as.POSIXlt(ix2)

hefs_forc_lamc<-readRDS('data/hefs_lamc_act-meteo/hefs_ens_forc_act-meteo_12z.rds')
syn_hefs_lamc<-readRDS(paste('data/lam_flow/ens_full_out/syn_hefs_flow_8510_',shefs_typ,'_ens_',m,'_12z.rds',sep=''))
syn_gefs_lamc<-readRDS('data/hefs_lamc_syn-meteo_v3/ens1/hefs_ens_forc_syn-meteo_12z.rds')

hefs_lam<-hefs_forc_lamc[,hefs_idx,]*1000
sgefs_lam<-syn_gefs_lamc[,hefs_idx,]*1000
shefs_lam<-syn_hefs_lamc[15:28,,]

#cold season index with quantile split
dec_mar<-sort(c(which(ixx2$mon>10),which(ixx2$mon<3)))
ob_ss<-sort(obs$LAMC1[dec_mar],index.return=T)
cut_ix<-(round(length(ob_ss$x)*pct[1])):(round(length(ob_ss$x)*pct[2]))
ss_ixx<-sort(ob_ss$ix[cut_ix])
ss_ix<-dec_mar[ss_ixx]
ss_idx<-ss_ix


hefs_lamc<-hefs_lam[,ss_idx,]
shefs_lamc<-shefs_lam[,ss_idx,]
sgefs_lamc<-sgefs_lam[,ss_idx,]
obs_lamc<-obs$LAMC1[ss_idx]

####################HEFS##################################

rnk_hefs<-array(NA,c(leads,length(ss_idx)))

for(j in 1:leads){
  for(i in 1:length(ss_idx)){
    rnk_hefs[j,i]<-ens_rank(hefs_lamc[,i,j],obs_lamc[i])
  }
}

rnk_hist_hefs<-array(NA,c(leads,ens_num+1))

for(j in 1:leads){
  for(i in 1:(ens_num+1)){
    rnk_hist_hefs[j,i]<-length(which(rnk_hefs[j,]==i))
  }
}

cumul_hist_hefs<-t(apply(rnk_hist_hefs,1,roll_sum))

####################syn-GEFS##################################

rnk_sgefs<-array(NA,c(leads,length(ss_idx)))

for(j in 1:leads){
  for(i in 1:length(ss_idx)){
    rnk_sgefs[j,i]<-ens_rank(sgefs_lamc[,i,j],obs_lamc[i])
  }
}

rnk_hist_sgefs<-array(NA,c(leads,ens_num+1))

for(j in 1:leads){
  for(i in 1:(ens_num+1)){
    rnk_hist_sgefs[j,i]<-length(which(rnk_sgefs[j,]==i))
  }
}

cumul_hist_sgefs<-t(apply(rnk_hist_sgefs,1,roll_sum))

####################syn-HEFS##################################

rnk_shefs<-array(NA,c(leads,length(ss_idx)))

for(j in 1:leads){
  for(i in 1:length(ss_idx)){
    rnk_shefs[j,i]<-ens_rank(shefs_lamc[j,i,],obs_lamc[i])
  }
}

rnk_hist_shefs<-array(NA,c(leads,ens_num+1))

for(j in 1:leads){
  for(i in 1:(ens_num+1)){
    rnk_hist_shefs[j,i]<-length(which(rnk_shefs[j,]==i))
  }
}

cumul_hist_shefs<-t(apply(rnk_hist_shefs,1,roll_sum))

#Cumulative Rank Histograms--------------------------------------------------------------
#combine plot

if(plt_par==T){
par(mfrow=c(1,length(lds)),mar=c(3,4,0,1),mgp=c(2,0.5,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3.5,3.5,0.5))
}

if(show_xlab==T){xlb<-'Ensemble Rank'} else{xlb<-''}
ylb<-c('Density',rep("",length(lds)-1))

alp_orange<-alpha('orange',alpha = 0.4)
alp_purp<-alpha('purple',alpha = 0.2)
my.gray<-alpha('gray',alpha=.7)

if(syn_plot=='sHEFS'){
  clr<-alp_orange
}

if(syn_plot=='sGEFS'){
  clr<-alp_purp
}

#par(mfrow=c(1,4))
for(i in 1:length(lds)){
  rank_hefs<-c()
  rank_syn<-c()
  for(j in 1:(ens_num+1)){
    if(syn_plot=='sHEFS'){
      rnks<-rep(j,rnk_hist_shefs[lds[i],j])}
    if(syn_plot=='sGEFS'){
      rnks<-rep(j,rnk_hist_sgefs[lds[i],j])}
    rank_syn<-c(rank_syn,rnks)
    rnks<-rep(j,rnk_hist_hefs[lds[i],j])
    rank_hefs<-c(rank_hefs,rnks)
  }
  hist(rank_hefs,breaks=seq(0.5,62.5,1),main='',xlab=xlb,
       ylab=ylb[i],freq=T,ylim=c(0,7.5*length(rank_hefs)/(ens_num+1)),col='gray')
  hist(rank_syn,breaks=seq(0.5,62.5,1),add=T,col=clr,border=my.gray)
  abline(h=round(length(rank_hefs)/(ens_num+1)),col='light blue',lty=2,lwd=2)
  legend('topleft',c('HEFS','SYN','UNIF'),col=c('gray',clr,'light blue'),lwd=c(6,6,2),lty=c(1,1,2),cex=1)
  #ri<-Ent(rnk_hist_hefs[lds[i],])
  #ch<-chi2_comp(rnk_hist_hefs[lds[i],],rnk_hist_syn_loess[lds[i],])
  #ch2<-format(round(ch[[1]],2),nsmall=2)
  #rix<-format(round(ri,5),nsmall=5)
  #pval<-format(round(ch[[2]],5),nsmall=5)
  #text(40,7*length(rank_hefs)/(ens_num+1),bquote(~chi^2==.(ch2)),cex=1.2,adj=0)
  #text(40,6*length(rank_hefs)/(ens_num+1),bquote(~p==.(pval)),cex=1.2,adj=0)
  #text(40,5*length(rank_hefs)/(ens_num+1),bquote(~phi==.(rix)),cex=1.2,adj=0)
  if(i==1){mtext(paste(lpcnt,'-',upcnt,'%',sep=''),side=2,line=4,font=2,cex=2.5)}
  if(show_leads==T){mtext(paste('Lead',lds[i]),side=3,font=2,cex=2.5)}
}
}

########################################END################################################