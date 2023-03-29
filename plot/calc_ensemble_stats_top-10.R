library(doParallel)
library(abind)

parallel::detectCores()

n.cores <- parallel::detectCores()-2
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

print(paste('start',Sys.time()))

nsamps<-10
nsets<-1
lds<-c(1,3,5,10)
outfile_hefs<-'corr_v9_sumrsamp30'
plt_lim<-3
max_idx<-c("2005-12-31", "1995-01-09", "1986-02-18", "1993-01-21", "1997-01-02",
             "2004-02-17", "1998-02-20","1995-03-10", "2002-12-21","2006-03-06")
  
leads<-14
ens_num<-61
cfs_taf<-2.29568411*10**-5 * 86400 / 1000
  
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

inf<-read.csv('~/syn_forecast_test/data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]
  
roll_sum<-function(x){out<-c();out[1]<-x[1];
for(i in 1:(length(x)-1)){out[i+1]<-out[i]+x[i+1]};return(out)}

####################calculate##################################
  
mns<-array(NA,c(length(max_idx),nsamps,leads))
meds<-array(NA,c(length(max_idx),nsamps,leads))
vr<-array(NA,c(length(max_idx),nsamps,leads))
  
cumul_mns<-array(NA,c(length(max_idx),nsamps,leads))
cumul_meds<-array(NA,c(length(max_idx),nsamps,leads))
cumul_var<-array(NA,c(length(max_idx),nsamps,leads))

out<-foreach(s = 1:nsets,.inorder=F) %dopar% {
  syn_forc_in<-readRDS(paste('~/syn_forecast_test/output/',outfile_hefs,'/syn_hefs_flow_8510_',outfile_hefs,'_ens_',s,'_12z.rds',sep=''))

  for(n in 1:nsamps){
    syn_forc<-syn_forc_in[,,,n]
      
    syn_hop<-syn_forc[1:14,,]
    syn_lam<-syn_forc[15:28,,]
    syn_uka<-syn_forc[29:42,,]
      
    syn_hop_out<-array(NA,c(dim(syn_hop)[2],(dim(syn_hop)[1]+1),ens_num))
    syn_lam_out<-array(NA,c(dim(syn_hop)[2],(dim(syn_hop)[1]+1),ens_num))
    syn_uka_out<-array(NA,c(dim(syn_hop)[2],(dim(syn_hop)[1]+1),ens_num))
      
    syn_hop_out[,1,]<-matrix(rep(obs[,1],ens_num),ncol=ens_num,byrow=F)
    syn_lam_out[,1,]<-matrix(rep(obs[,2],ens_num),ncol=ens_num,byrow=F)
    syn_uka_out[,1,]<-matrix(rep(obs[,3],ens_num),ncol=ens_num,byrow=F)
      
    for(l in 1:14){
      syn_hop_out[,(l+1),]<-syn_hop[l,c((l+1):dim(syn_hop)[2],rep(dim(syn_hop)[2],(l))),1:ens_num]
      syn_lam_out[,(l+1),]<-syn_lam[l,c((l+1):dim(syn_hop)[2],rep(dim(syn_hop)[2],(l))),1:ens_num]
      syn_uka_out[,(l+1),]<-syn_uka[l,c((l+1):dim(syn_hop)[2],rep(dim(syn_hop)[2],(l))),1:ens_num]
    }
      
    for(k in 1:length(max_idx)){
      idx<-which(ix==max_idx[k])
        
      mns[k,n,]<-apply(syn_lam[,idx,],1,mean)*cfs_taf
      meds[k,n,]<-apply(syn_lam[,idx,],1,median)*cfs_taf
      vr[k,n,]<-apply(syn_lam[,idx,]*cfs_taf,1,var)
        
      for(j in 1:leads){
        mn<-apply(syn_lam_out[idx-j,,],2,roll_sum)
        cumul_mns[k,n,j]<-mean(mn[j+1,])*cfs_taf
        cumul_meds[k,n,j]<-median(mn[j+1,])*cfs_taf
        cumul_var[k,n,j]<-var(mn[j+1,]*cfs_taf)
      }
    }
  }
  return(list(mns,meds,vr,cumul_mns,cumul_meds,cumul_var))
}


m<-nsets*nsamps

mns_out<-array(NA,c(length(max_idx),m,leads))
meds_out<-array(NA,c(length(max_idx),m,leads))
vr_out<-array(NA,c(length(max_idx),m,leads))

cumul_mns_out<-array(NA,c(length(max_idx),m,leads))
cumul_meds_out<-array(NA,c(length(max_idx),m,leads))
cumul_var_out<-array(NA,c(length(max_idx),m,leads))

samp_vec<-matrix(1:(nsets*nsamps),ncol=nsamps,byrow=T)

for(i in 1:nsets){
  mns_out[,samp_vec[i,],]<-out[[i]][[1]]
  meds_out[,samp_vec[i,],]<-out[[i]][[2]]
  vr_out[,samp_vec[i,],]<-out[[i]][[3]]
  
  cumul_mns_out[,samp_vec[i,],]<-out[[i]][[4]]
  cumul_meds_out[,samp_vec[i,],]<-out[[i]][[5]]
  cumul_var_out[,samp_vec[i,],]<-out[[i]][[6]]
}

saveRDS(mns_out,paste('~/syn_forecast_test/output/',outfile_hefs,'/mns_out.rds',sep=''))
saveRDS(meds_out,paste('~/syn_forecast_test/output/',outfile_hefs,'/meds_out.rds',sep=''))
saveRDS(vr_out,paste('~/syn_forecast_test/output/',outfile_hefs,'/vr_out.rds',sep=''))

saveRDS(cumul_mns_out,paste('~/syn_forecast_test/output/',outfile_hefs,'/cumul_mns_out.rds',sep=''))
saveRDS(cumul_meds_out,paste('~/syn_forecast_test/output/',outfile_hefs,'/cumul_meds_out.rds',sep=''))
saveRDS(cumul_var_out,paste('~/syn_forecast_test/output/',outfile_hefs,'/cumul_var_out.rds',sep=''))

print(paste('end',Sys.time()))

rm(list=ls());gc()

###############################END###################################
  