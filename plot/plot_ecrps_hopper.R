
plt_ecrps<-function(nsamps,nsets,lds,lpcnt,upcnt,show_leads,show_xlab,plt_par,outfile,ylm,yinc){

#nsamps<-10
#nsets<-10
#lds<-c(1,3,5,10)
#lpcnt<-90
#upcnt<-100
#show_leads<-F
#show_xlab<-F
#plt_par<-T
#outfile<-'v10_wtsamp30'
#ylm<-2000
#yinc<-500
  
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

ecrps_hefs<-array(NA,c(leads,length(ss_idx)))
ecrps_shefs<-array(NA,c(nsamps,leads,length(ss_idx)))
ecrps_sgefs<-array(NA,c(nsamps,leads,length(ss_idx)))

for(i in 1:leads){
  for(k in 1:length(ss_idx)){
    ecrps_hefs[i,k]<-eCRPS(hefs_lamc[,k,i],obs_lamc[k])
  }
}

for(n in 1:nsamps){
  for(i in 1:leads){
    for(k in 1:length(ss_idx)){
      ecrps_shefs[n,i,k]<-eCRPS(shefs_lamc[i,k,,n],obs_lamc[k])
      ecrps_sgefs[n,i,k]<-eCRPS(sgefs_lamc[,k,i,n],obs_lamc[k])
    }
  }
}

#Cumulative Rank Histograms--------------------------------------------------------------
#combine plot

if(plt_par==T){
  par(mar=c(5.5,3,3,3),cex.main = 1,cex.axis = 1,cex.lab = 1,mfrow=c(1,4))
}

ylb<-c('Ensemble CRPS',rep("",length(lds)-1))
if(show_xlab==T){xlb<-'Ensemble Rank'} else{xlb<-''}

for(i in 1:length(lds)){
  boxplot(ecrps_hefs[lds[i],],as.vector(ecrps_shefs[,lds[i],]),as.vector(ecrps_sgefs[,lds[i],]),
          main = '',
          at = c(0.5,2,3.5),
          names = c('','',''),
          las = 1,
          boxwex = 0.5,
          ylab = ylb[i],
          notch = F,
          col = c('gray','orange','purple'),
          ylim = c(0,ylm),
          xlim = c(0,4),
          outline = F,
          axes = F
  )
  if(show_xlab==T){axis(1,at=c(0.5,2,3.5),labels = c('HEFS','sHEFS','sGEFS'))} else{axis(1,at=c(0.5,2,3.5),labels=F)}
  if(i==1){
    text(0.1,(0.9*ylm),bquote(~n==.(length(ss_idx))),cex=2,adj=0,font=3)
    axis(2,at=seq(0,ylm,yinc),labels=seq(0,ylm,yinc))} else{axis(2,at=seq(0,ylm,yinc),labels=F)}
  if(i==1){mtext(paste(lpcnt,'-',upcnt,'%',sep=''),side=2,line=5,font=2,cex=2.5)}
  if(show_leads==T){mtext(paste('Lead',lds[i]),side=3,line=1,font=2,cex=2.5)}
  box(which='plot')
}

}

########################################END################################################