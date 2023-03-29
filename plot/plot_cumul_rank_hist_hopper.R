
plt_cumul_rank_hist<-function(nsamps,nsets,lds,lpcnt,upcnt,show_leads,show_xlab,show_xaxs,show_legend,plt_par,outfile){

#nsamps<-10
#nsets<-10
#lds<-c(1,3,5,10)
#lpcnt<-90
#upcnt<-100
#show_leads<-F
#show_xlab<-F
#plt_par<-T
#outfile<-'v10_wtsamp30'
  
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

rnk_sgefs<-array(NA,c(nsamps,leads,length(ss_idx)))

for(n in 1:nsamps){
  for(j in 1:leads){
    for(i in 1:length(ss_idx)){
      rnk_sgefs[n,j,i]<-ens_rank(sgefs_lamc[,i,j,n],obs_lamc[i])
    }
  }
}

rnk_hist_sgefs<-array(NA,c(nsamps,leads,ens_num+1))

for(n in 1:nsamps){
  for(j in 1:leads){
    for(i in 1:(ens_num+1)){
      rnk_hist_sgefs[n,j,i]<-length(which(rnk_sgefs[n,j,]==i))
    }
  }
}

cumul_hist_sgefs<-array(NA,c(nsamps,leads,ens_num+1))

for(n in 1:nsamps){
  cumul_hist_sgefs[n,,]<-t(apply(rnk_hist_sgefs[n,,],1,roll_sum))
}

####################syn-HEFS##################################

rnk_shefs<-array(NA,c(nsamps,leads,length(ss_idx)))

for(n in 1:nsamps){
  for(j in 1:leads){
    for(i in 1:length(ss_idx)){
      rnk_shefs[n,j,i]<-ens_rank(shefs_lamc[j,i,,n],obs_lamc[i])
    }
  }
}

rnk_hist_shefs<-array(NA,c(nsamps,leads,ens_num+1))

for(n in 1:nsamps){
  for(j in 1:leads){
    for(i in 1:(ens_num+1)){
      rnk_hist_shefs[n,j,i]<-length(which(rnk_shefs[n,j,]==i))
    }
  }
}

cumul_hist_shefs<-array(NA,c(nsamps,leads,ens_num+1))

for(n in 1:nsamps){
  cumul_hist_shefs[n,,]<-t(apply(rnk_hist_shefs[n,,],1,roll_sum))
}

#Cumulative Rank Histograms--------------------------------------------------------------
#combine plot

if(plt_par==T){
  par(mfrow=c(1,length(lds)),mar=c(3,4,0,1),mgp=c(2,0.5,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3.5,3.5,0.5))
}

ylb<-c('Cumulative Density',rep("",length(lds)-1))
alp_orange<-alpha('orange',alpha = 0.15)
alp_purp<-alpha('purple',alpha = 0.1)

for(i in 1:length(lds)){
  plot(1:(ens_num+1),cumul_hist_hefs[lds[i],]/length(ss_idx),main='',
       xlab='',ylab=ylb[i],xlim=c(0,62),ylim=c(0,1),type='l',axes=F)
  if(nsamps>1){
    for(n in 1:nsamps){
      lines(1:(ens_num+1),cumul_hist_shefs[n,lds[i],]/length(ss_idx),col=alp_orange,lwd=2)
      lines(1:(ens_num+1),cumul_hist_sgefs[n,lds[i],]/length(ss_idx),col=alp_purp,lwd=2)
    }
  }
  #lines(1:(ens_num+1),cumul_hist_shefs[1,lds[i],]/length(ss_idx),col='dark orange',lwd=2)
  #lines(1:(ens_num+1),cumul_hist_sgefs[1,lds[i],]/length(ss_idx),col='purple',lwd=2)
  lines(1:(ens_num+1),apply(cumul_hist_shefs[,lds[i],],2,median)/length(ss_idx),col='dark orange',lwd=2)
  lines(1:(ens_num+1),apply(cumul_hist_sgefs[,lds[i],],2,median)/length(ss_idx),col='purple',lwd=2)
  lines(1:(ens_num+1),cumul_hist_hefs[lds[i],]/length(ss_idx),lwd=3)
  abline(0,(1/ens_num),col='light gray',lty=2,lwd=2)
  axis(2,at=seq(0,1,0.2),labels=F)
  axis(1,at=seq(0,62,10),labels=F)
  if(i==1){
    text(1,0.9,bquote(~n==.(length(ss_idx))),cex=2,adj=0,font=3)
    axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2))}
  if(i == length(lds) & show_legend==T){legend('topleft',c('HEFS','syn-HEFS','syn-GEFS'),col=c('black','orange','purple'),lwd=c(2,2,2),cex=2)}
  if(show_xlab==T){
    mtext('Ensemble Rank',side=1,line=2.25,font=2,cex=1.4)}
  if(show_xaxs==T){
    axis(1,at=seq(0,62,10),labels = seq(0,62,10))
  }
  #ri<-Ent(rnk_hist_hefs[lds[i],])
  #ri2<-Ent(rnk_hist_shefs[lds[i],])
  #ri3<-Ent(rnk_hist_sgefs[lds[i],])
  #ch<-chi2_comp(rnk_hist_hefs[lds[i],],rnk_hist_sgefs[lds[i],])
  #ch2<-format(round(ch[[1]],2),nsmall=1)
  #rix<-format(round(ri,3),nsmall=3)
  #rix2<-format(round(ri2,3),nsmall=3)
  #rix3<-format(round(ri3,3),nsmall=3)
  #pval<-format(round(ch[[2]],3),nsmall=4)
  #text(40,0.25,bquote(~chi^2==.(ch2)),cex=1.5,adj=0)
  #text(40,0.15,bquote(~p==.(pval)),cex=1.5,adj=0)
  #text(40,0.05,bquote(~phi==.(rix)),cex=1.5,adj=0)
  #text(40,0.25,bquote(~phi==.(rix)),cex=1.5,adj=0)
  #text(40,0.15,bquote(~phi==.(rix2)),cex=1.5,adj=0,col='orange')
  #text(40,0.05,bquote(~phi==.(rix3)),cex=1.5,adj=0,col='purple')
  if(i==1){mtext(paste(lpcnt,'-',upcnt,'%',sep=''),side=2,line=5,font=2,cex=2.5)}
  if(show_leads==T){mtext(paste('Lead',lds[i]),side=3,font=2,line=1,cex=2.5)}
  box(which='plot')
}

}

########################################END################################################