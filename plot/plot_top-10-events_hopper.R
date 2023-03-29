
plt_top10_events<-function(nsamps,nsets,lds,show_xlab,show_extremes,show_mn_xlab,show_var_xlab,calc_shefs,outfile_gefs,outfile_hefs,syn_forc_typ,mn_lim,vr_lim,max_idx){

#nsamps<-10
#nsets<-10
#lds<-c(1,3,5,10)
#show_xlab<-F
#show_extremes<-T
#calc_shefs<-F
#syn_forc_typ<-'shefs'
#show_mn_xlab<-F
#show_var_xlab<-T
#outfile_hefs<-'v9_wtsamp30'
#outfile_gefs<-'v3'
#mn_lim<-3
#vr_lim<-4
#max_idx<-c("2005-12-31", "1995-01-09", "1986-02-18", "1993-01-21", "1997-01-02",
           #"2004-02-17", "1998-02-20","1995-03-10", "2002-12-21","2006-03-06")

  
setwd('h:/firo_lamc/')
library(scales)
library(abind)
library(RColorBrewer)

leads<-14
ens_num<-61
cfs_taf<-2.29568411*10**-5 * 86400 / 1000

source('plot/forecast_verification_functions.R')
inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]

ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

lamc_hefs<-readRDS('data/hefs_lamc_act-meteo/hefs_dly_mean_act-meteo_12z.rds')
lamc_hefs<-lamc_hefs[hefs_idx,,]*1000

roll_sum<-function(x){out<-c();out[1]<-x[1];
for(i in 1:(length(x)-1)){out[i+1]<-out[i]+x[i+1]};return(out)}

####################calculate##################################

if(calc_shefs==T & syn_forc_typ=='shefs'){
m<-nsamps * nsets

mns<-array(NA,c(length(max_idx),m,leads))
meds<-array(NA,c(length(max_idx),m,leads))
vr<-array(NA,c(length(max_idx),m,leads))

cumul_mns<-array(NA,c(length(max_idx),m,leads))
cumul_meds<-array(NA,c(length(max_idx),m,leads))
cumul_var<-array(NA,c(length(max_idx),m,leads))

samp_vec<-matrix(1:(nsets*nsamps),ncol=nsamps,byrow=T)

for(s in 1:nsets){
    syn_forc_in<-readRDS(paste('z:/syn_forecast_test/output/',outfile_hefs,'/syn_hefs_flow_8510_',outfile_hefs,'_ens_',s,'_12z.rds',sep=''))

  
  for(n in 1:nsamps){
    syn_forc<-syn_forc_in[,,,n]

    syn_lam<-syn_forc[15:28,,]
  
    syn_lam_out<-array(NA,c(dim(syn_hop)[2],(dim(syn_lam)[1]+1),ens_num))
    syn_lam_out[,1,]<-matrix(rep(obs[,2],ens_num),ncol=ens_num,byrow=F)

  
    for(l in 1:14){
      syn_lam_out[,(l+1),]<-syn_lam[l,c((l+1):dim(syn_lam)[2],rep(dim(syn_lam)[2],(l))),1:ens_num]
    }
  
    for(k in 1:length(max_idx)){
      idx<-which(ix==max_idx[k])
    
      mns[k,samp_vec[s,n],]<-apply(syn_lam[,idx,],1,mean)*cfs_taf
      meds[k,samp_vec[s,n],]<-apply(syn_lam[,idx,],1,median)*cfs_taf
      vr[k,samp_vec[s,n],]<-apply(syn_lam[,idx,]*cfs_taf,1,var)
    
      for(j in 1:leads){
        mn<-apply(syn_lam_out[idx-j,,],2,roll_sum)
        cumul_mns[k,samp_vec[s,n],j]<-mean(mn[j+1,])*cfs_taf
        cumul_meds[k,samp_vec[s,n],j]<-median(mn[j+1,])*cfs_taf
        cumul_var[k,samp_vec[s,n],j]<-var(mn[j+1,]*cfs_taf)
      }
    }
  }
}
}

if(syn_forc_typ=='sgefs'){
  m<-nsamps * nsets
  
  mns<-array(NA,c(length(max_idx),m,leads))
  meds<-array(NA,c(length(max_idx),m,leads))
  vr<-array(NA,c(length(max_idx),m,leads))
  
  cumul_mns<-array(NA,c(length(max_idx),m,leads))
  cumul_meds<-array(NA,c(length(max_idx),m,leads))
  cumul_var<-array(NA,c(length(max_idx),m,leads))
  
  samp_vec<-matrix(1:(nsets*nsamps),ncol=nsamps,byrow=T)
  
  for(s in 1:nsets){
    for(n in 1:nsamps){
      syn_forc<-readRDS(paste('data/hefs_lamc_syn-meteo_',outfile_gefs,'/ens',n,'/hefs_ens_forc_syn-meteo_12z.rds',sep=''))

      syn_lam_in<-syn_forc[,hefs_idx,]*1000
      syn_lam<-array(NA,c(leads,length(hefs_idx),ens_num))
      
      for(i in 1:leads){
        syn_lam[i,,]<-t(syn_lam_in[,,i])
      }
      
      syn_lam_out<-array(NA,c(dim(syn_lam)[2],(dim(syn_lam)[1]+1),ens_num))
      
      syn_lam_out[,1,]<-matrix(rep(obs[,2],ens_num),ncol=ens_num,byrow=F)
      
      for(l in 1:14){
        syn_lam_out[,(l+1),]<-syn_lam[l,c((l+1):dim(syn_lam)[2],rep(dim(syn_lam)[2],(l))),1:ens_num]
      }
      
      
      for(k in 1:length(max_idx)){
        idx<-which(ix==max_idx[k])
        
        mns[k,samp_vec[s,n],]<-apply(syn_lam[,idx,],1,mean)*cfs_taf
        meds[k,samp_vec[s,n],]<-apply(syn_lam[,idx,],1,median)*cfs_taf
        vr[k,samp_vec[s,n],]<-apply(syn_lam[,idx,]*cfs_taf,1,var)
        
        for(j in 1:leads){
          mn<-apply(syn_lam_out[idx-j,,],2,roll_sum)
          cumul_mns[k,samp_vec[s,n],j]<-mean(mn[j+1,])*cfs_taf
          cumul_meds[k,samp_vec[s,n],j]<-median(mn[j+1,])*cfs_taf
          cumul_var[k,samp_vec[s,n],j]<-var(mn[j+1,]*cfs_taf)
        }
      }
    }
  }
}

if(calc_shefs==F & syn_forc_typ=='shefs'){
  mns<-readRDS(paste('z:/syn_forecast_test/output/',outfile_hefs,'/mns_out.rds',sep=''))
  meds<-readRDS(paste('z:/syn_forecast_test/output/',outfile_hefs,'/meds_out.rds',sep=''))
  vr<-readRDS(paste('z:/syn_forecast_test/output/',outfile_hefs,'/vr_out.rds',sep=''))
  
  cumul_mns<-readRDS(paste('z:/syn_forecast_test/output/',outfile_hefs,'/cumul_mns_out.rds',sep=''))
  cumul_meds<-readRDS(paste('z:/syn_forecast_test/output/',outfile_hefs,'/cumul_meds_out.rds',sep=''))
  cumul_var<-readRDS(paste('z:/syn_forecast_test/output/',outfile_hefs,'/cumul_var_out.rds',sep=''))
  
}


#---------------------------Calculate HEFS metrics---------------------------------------------

mns_hefs<-array(NA,c(length(max_idx),leads))
meds_hefs<-array(NA,c(length(max_idx),leads))
var_hefs<-array(NA,c(length(max_idx),leads))

cumul_mns_hefs<-array(NA,c(length(max_idx),leads))
cumul_meds_hefs<-array(NA,c(length(max_idx),leads))
cumul_var_hefs<-array(NA,c(length(max_idx),leads))

hefs_out<-abind(matrix(rep(obs$LAMC1,ens_num),ncol=ens_num,byrow = F),lamc_hefs,along=2)

for(k in 1:length(max_idx)){
  idx<-which(ix==max_idx[k])
  
  for(s in 1:leads){
    mns_hefs[k,s]<-mean(lamc_hefs[idx-s,s,])*cfs_taf
    meds_hefs[k,s]<-median(lamc_hefs[idx-s,s,])*cfs_taf
    var_hefs[k,s]<-var(lamc_hefs[idx-s,s,]*cfs_taf)
    
    mn<-apply(hefs_out[idx-s,,],2,roll_sum)
    cumul_mns_hefs[k,s]<-mean(mn[s+1,])*cfs_taf
    cumul_meds_hefs[k,s]<-median(mn[s+1,])*cfs_taf
    cumul_var_hefs[k,s]<-var(mn[s+1,]*cfs_taf)
  }
}

####################plot##################################

#cumulative combined plot - top 10 events
rng<-1.5

if(show_extremes==T){
  rng<-0
}

cumul_obs<-array(NA,c(leads,length(max_idx)))

for(i in 1:length(max_idx)){
  idx<-which(ix==max_idx[i])
  
  for(j in 1:leads){
    cumul_obs[j,i]<-sum(obs$LAMC1[(idx-j):idx])*cfs_taf
  }
}

m<-nsets * nsamps

bplot_arr<-cbind(t(cumul_mns[,,lds[1]]/matrix(rep(cumul_mns_hefs[,lds[1]],m),ncol=m,byrow=F)),
                 t(cumul_mns[,,lds[2]]/matrix(rep(cumul_mns_hefs[,lds[2]],m),ncol=m,byrow=F)),
                 t(cumul_mns[,,lds[3]]/matrix(rep(cumul_mns_hefs[,lds[3]],m),ncol=m,byrow=F)),
                 t(cumul_mns[,,lds[4]]/matrix(rep(cumul_mns_hefs[,lds[4]],m),ncol=m,byrow=F)))  


g<-brewer.pal(10,'Set3')

boxplot(bplot_arr,
        main = '',#paste('Cumulative Ensemble Mean Comparison: 10 Extreme Events'),
        at = c(seq(0.5,5,0.5),seq(7,11.5,0.5),seq(13.5,18,0.5),seq(20,24.5,0.5)),
        names = c(rep("",6),paste(lds[1],'day Lead',sep='-'),rep("",8),paste(lds[2],'day Lead',sep='-'),rep("",8),paste(lds[3],'day Lead',sep='-'),rep("",8),paste(lds[4],'day Lead',sep='-'),rep("",6)),
        las = 1,
        boxwex = .3,
        notch = F,
        col=rep(g,4),
        xlab='',
        ylab='',
        ylim = c(0,mn_lim),
        xlim = c(0,30),
        outline = F,
        axes=F,
        range=rng
)
points(c(seq(0.5,5,0.5),seq(7,11.5,0.5),seq(13.5,18,0.5),seq(20,24.5,0.5)),
       c(cumul_obs[lds[1],]/cumul_mns_hefs[,lds[1]],cumul_obs[lds[2],]/cumul_mns_hefs[,lds[2]],
         cumul_obs[lds[3],]/cumul_mns_hefs[,lds[3]],cumul_obs[lds[4],]/cumul_mns_hefs[,lds[4]]),
       col='black',pch=17,cex=1.5)

abline(h=1,lty=2,col='dark green')
mtext('Cumulative Inflow Ratio',side=2,line=3,font=2,cex=2)

if(show_mn_xlab==T){
  axis(1,at=c(3,9.5,16,22.5),labels = c(paste(lds[1],'day Lead',sep='-'),paste(lds[2],'day Lead',sep='-'),paste(lds[3],'day Lead',sep='-'),paste(lds[4],'day Lead',sep='-')))
}else{axis(1,at=c(3,9.5,16,22.5),labels =F)}

axis(2,las=2,at=seq(0,2.5,0.5),labels = seq(0,2.5,0.5))
box(which='plot')
legend('topright',c(max_idx,'HEFS'),lwd=c(rep(6,10),1),lty=c(rep(1,10),2),text.font=c(rep(1,10),2),col=c(g,'dark green'),cex=1.5)
legend('topleft','Observed',pch=17,col='black',cex=1.5)


#---------------------------------------------------------------------------------------
#cumulative variance
bplot_arr_var<-cbind(t(cumul_var[,,lds[1]]/matrix(rep(cumul_var_hefs[,lds[1]],m),ncol=m,byrow=F)),
                     t(cumul_var[,,lds[2]]/matrix(rep(cumul_var_hefs[,lds[2]],m),ncol=m,byrow=F)),
                     t(cumul_var[,,lds[3]]/matrix(rep(cumul_var_hefs[,lds[3]],m),ncol=m,byrow=F)),
                     t(cumul_var[,,lds[4]]/matrix(rep(cumul_var_hefs[,lds[4]],m),ncol=m,byrow=F)))  


boxplot(bplot_arr_var,
        main = '',#paste('Cumulative Ensemble Variance Comparison: 10 Extreme Events'),
        at = c(seq(0.5,5,0.5),seq(7,11.5,0.5),seq(13.5,18,0.5),seq(20,24.5,0.5)),
        names = c(rep("",6),paste(lds[1],'day Lead',sep='-'),rep("",8),paste(lds[2],'day Lead',sep='-'),rep("",8),paste(lds[3],'day Lead',sep='-'),rep("",8),paste(lds[4],'day Lead',sep='-'),rep("",6)),
        las = 1,
        boxwex = .3,
        notch = F,
        col=rep(g,4),
        xlab='',
        ylab='',
        ylim = c(0,vr_lim),
        xlim = c(0,30),
        outline = F,
        axes=F,
        range=rng
)
abline(h=1,lty=2,col='dark green')

if(show_var_xlab==T){
  axis(1,at=c(3,9.5,16,22.5),labels = c(paste(lds[1],'day Lead',sep='-'),paste(lds[2],'day Lead',sep='-'),paste(lds[3],'day Lead',sep='-'),paste(lds[4],'day Lead',sep='-')))
}else{axis(1,at=c(3,9.5,16,22.5),labels =F)}

axis(2,las=2,at=seq(0,4,0.5),labels = seq(0,4,0.5))
mtext('Cumulative Variance Ratio',side=2,line=3,font=2,cex=2)
box(which='plot')
legend('topright',c(max_idx,'HEFS'),lwd=c(rep(6,10),1),lty=c(rep(1,10),2),text.font=c(rep(1,10),2),col=c(g,'dark green'),cex=1.5)

}

########################################END################################################