###################################################################################
# Code to process meteorologica TMAX, TMIN, and PRECIP data
# 
# Supports manuscript: 'A generalized, multivariate approach to generate synthetic 
# short-to-medium range hydro-meteorological forecasts across locations, variables, 
# and lead times', Authors - Brodeur, Z, and Steinschneider, S.
# 
# Created by Zach Brodeur, 15 December 2020
# -Updated 18 April 2021: Standardized and clarified for replicability
#####################################################################################

#Required packages and drive
setwd('h:/firo_lamc/temp_precip/')

cond_mean<-function(Q,f){
  cmean <- mean(f) + (cov(Q,f)/(sd(Q))^2) * (Q - mean(Q))
  return(cmean)
}

#--------------------------------------------------------------------------------
#1) Process observed data
obs<-read.csv('data/Russian River Old Historical Forcings - CNRFC_2019_Forcings.csv')
ix<-seq(as.Date('1948-10-02'),as.Date('2017-09-30'),'day')

map_hop<-obs$HOPC1LOF
map_hop<-levels(map_hop)[map_hop]
map_hop<-as.numeric(map_hop[4:(length(map_hop)-2)])
map_hop<-map_hop*25.4 #in to mm conversion
map_hop_mat<-matrix(map_hop,ncol=4,byrow=T)
row.names(map_hop_mat)<-as.character(ix)
colnames(map_hop_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(map_hop_mat,'data/basic/map_hop_mat.rds')

mat_hop<-obs$HOPC1LOF.1
mat_hop<-levels(mat_hop)[mat_hop]
mat_hop<-as.numeric(mat_hop[4:(length(mat_hop)-2)]) 
mat_hop<-((mat_hop - 32) * (5/9))
mat_hop_mat<-matrix(mat_hop,ncol=4,byrow=T)
row.names(mat_hop_mat)<-as.character(ix)
colnames(mat_hop_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(mat_hop_mat,'data/basic/mat_hop_mat.rds')

map_lam<-obs$LAMC1HOF
map_lam<-levels(map_lam)[map_lam]
map_lam<-as.numeric(map_lam[4:(length(map_lam)-2)])
amp_lam<-map_lam*25.4
map_lam_mat<-matrix(map_lam,ncol=4,byrow=T)
row.names(map_lam_mat)<-as.character(ix)
colnames(map_lam_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(map_lam_mat,'data/basic/map_lam_mat.rds')

mat_lam<-obs$LAMC1HOF.1
mat_lam<-levels(mat_lam)[mat_lam]
mat_lam<-as.numeric(mat_lam[4:(length(mat_lam)-2)]) 
mat_lam<-((mat_lam - 32) * (5/9))
mat_lam_mat<-matrix(mat_lam,ncol=4,byrow=T)
row.names(mat_lam_mat)<-as.character(ix)
colnames(mat_lam_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(mat_lam_mat,'data/basic/mat_lam_mat.rds')

map_uka<-obs$UKAC1HOF
map_uka<-levels(map_uka)[map_uka]
map_uka<-as.numeric(map_uka[4:(length(map_uka)-2)])
map_uka<-map_uka*25.4
map_uka_mat<-matrix(map_uka,ncol=4,byrow=T)
row.names(map_uka_mat)<-as.character(ix)
colnames(map_uka_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(map_uka_mat,'data/basic/map_uka_mat.rds')

mat_uka<-obs$UKAC1HOF.1
mat_uka<-levels(mat_uka)[mat_uka]
mat_uka<-as.numeric(mat_uka[4:(length(mat_uka)-2)]) 
mat_uka<-((mat_uka - 32) * (5/9))
mat_uka_mat<-matrix(mat_uka,ncol=4,byrow=T)
row.names(mat_uka_mat)<-as.character(ix)
colnames(mat_uka_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(mat_uka_mat,'data/basic/mat_uka_mat.rds')

map_tot<-obs$Area.Weighted
map_tot<-levels(map_tot)[map_tot]
map_tot<-as.numeric(map_tot[4:(length(map_tot)-2)])
map_tot<-map_tot*25.4
map_tot_mat<-matrix(map_tot,ncol=4,byrow=T)
row.names(map_tot_mat)<-as.character(ix)
colnames(map_tot_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(map_tot_mat,'data/basic/map_tot_mat.rds')

mat_tot<-obs$Area.Weighted.1
mat_tot<-levels(mat_tot)[mat_tot]
mat_tot<-as.numeric(mat_tot[4:(length(mat_tot)-2)]) 
mat_tot<-((mat_tot - 32) * (5/9))
mat_tot_mat<-matrix(mat_tot,ncol=4,byrow=T)
row.names(mat_tot_mat)<-as.character(ix)
colnames(mat_tot_mat)<-c('00:00Z','06:00Z','12:00Z','18:00Z')
saveRDS(mat_tot_mat,'data/basic/mat_tot_mat.rds')

ob_map<-matrix(rep(map_tot_mat,17),ncol=68,byrow=F)[,-c(1,66:68)]
ob_mat<-matrix(rep(mat_tot_mat,17),ncol=68,byrow=F)[,-c(1,66:68)]

#--------------------------------------------------------------------------------
#2) Process GEF v10 Forecast Data and assign to variables
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-12-31'),'day') #same length
ix3<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx3<-as.POSIXlt(ix3)
cold<-sort(c(which(ixx3$mon>8),which(ixx3$mon<3)))
warm<-1:length(ixx3)
warm<-warm[-c(cold)]

gefs_tmax<-read.csv('data/gefs_tmax_6h_8517.csv')
gefs_tmin<-read.csv('data/gefs_tmin_6h_8517.csv')
gefs_prcp<-read.csv('data/gefs_prcp_6h_8517.csv')

gefs_tmax_fit<-gefs_tmax[which(ix2=='1985-01-01'):which(ix2=='2017-09-30'),]
dim(gefs_tmax_fit)

gefs_tmin_fit<-gefs_tmin[which(ix2=='1985-01-01'):which(ix2=='2017-09-30'),]
dim(gefs_tmin_fit)

gefs_prcp_fit<-gefs_prcp[which(ix2=='1985-01-01'):which(ix2=='2017-09-30'),]
dim(gefs_prcp_fit)

tmax_err<-gefs_tmax_fit
tmin_err<-gefs_tmin_fit
prcp_err<-gefs_prcp_fit

tmax_forc<-gefs_tmax_fit
tmin_forc<-gefs_tmin_fit
prcp_forc<-gefs_prcp_fit

mat_comb<-matrix(rep(mat_tot_mat,16),ncol=64)
mat_comb<-cbind(mat_comb,mat_comb[,1])
mat_comb<-mat_comb[,-c(1)]
colnames(mat_comb)<-seq(6,384,6)

map_comb<-matrix(rep(map_tot_mat,16),ncol=64)
map_comb<-cbind(map_comb,map_comb[,1])
map_comb<-map_comb[,-c(1)]
colnames(map_comb)<-seq(6,384,6)

#--------------------------------------------------------------------------------
#3) Calculate matrices of forecast error

#errors are observed - forecast
obs_mat_fit<-mat_comb[which(ix=='1985-01-01'):which(ix=='2017-09-30'),]
dim(obs_mat_fit)

obs_map_fit<-map_comb[which(ix=='1985-01-01'):which(ix=='2017-09-30'),]
dim(obs_map_fit)

tmax_err[,2:4]<-obs_mat_fit[,1:3]-gefs_tmax_fit[,2:4]
tmin_err[,2:4]<-obs_mat_fit[,1:3]-gefs_tmin_fit[,2:4]
prcp_err[,2:4]<-obs_map_fit[,1:3]-gefs_prcp_fit[,2:4]

tmax_forc[,2:4]<-gefs_tmax_fit[,2:4]
tmin_forc[,2:4]<-gefs_tmin_fit[,2:4]
prcp_forc[,2:4]<-gefs_prcp_fit[,2:4]

idx1<-seq(4,(15*4),4)
idx2<-seq(5,(15*4+1),4)

for(i in 1:15){
  tmx<-gefs_tmax_fit[1:(length(ix3)-i),idx2[i]:(idx2[i]+3)]
  tmx_f<-rbind(matrix(0,ncol=4,nrow=i),as.matrix(tmx))
  tmax_err[,idx2[i]:(idx2[i]+3)]<-obs_mat_fit[,idx1[i]:(idx1[i]+3)] - tmx_f
  tmax_forc[,idx2[i]:(idx2[i]+3)]<-tmx_f
  tmax_err[1:i,idx2[i]:(idx2[i]+3)]<-0
  
  tmn<-gefs_tmin_fit[1:(length(ix3)-i),idx2[i]:(idx2[i]+3)]
  tmn_f<-rbind(matrix(0,ncol=4,nrow=i),as.matrix(tmn))
  tmin_err[,idx2[i]:(idx2[i]+3)]<-obs_mat_fit[,idx1[i]:(idx1[i]+3)] - tmn_f
  tmin_forc[,idx2[i]:(idx2[i]+3)]<-tmn_f
  tmin_err[1:i,idx2[i]:(idx2[i]+3)]<-0
  
  prp<-gefs_prcp_fit[1:(length(ix3)-i),idx2[i]:(idx2[i]+3)]
  prp_f<-rbind(matrix(0,ncol=4,nrow=i),as.matrix(prp))
  prcp_err[,idx2[i]:(idx2[i]+3)]<-obs_map_fit[,idx1[i]:(idx1[i]+3)] - prp_f
  prcp_forc[,idx2[i]:(idx2[i]+3)]<-prp_f
  prcp_err[1:i,idx2[i]:(idx2[i]+3)]<-0
}

tmax_err[,65]<-obs_mat_fit[,64]-c(rep(0,16),gefs_tmax_fit[1:(length(ix3)-16),65])
tmax_err[1:16,65]<-0
tmin_err[,65]<-obs_mat_fit[,64]-c(rep(0,16),gefs_tmin_fit[1:(length(ix3)-16),65])
tmin_err[1:16,65]<-0
prcp_err[,65]<-obs_map_fit[,64]-c(rep(0,16),gefs_prcp_fit[1:(length(ix3)-16),65])
prcp_err[1:16,65]<-0

prcp_cmean<-array(NA,dim(obs_map_fit))
prcp_err_cmean<-array(NA,dim(prcp_err))

prcp_ls<-array(NA,dim(obs_map_fit))
prcp_err_ls<-array(NA,dim(prcp_err))
                             
for(i in 1:64){
  cmean_cold<-cond_mean(obs_map_fit[cold,i],prcp_forc[cold,(i+1)])
  cmean_warm<-cond_mean(obs_map_fit[warm,i],prcp_forc[warm,(i+1)])
  cmean_cold[obs_map_fit[cold,i]==0]<-0
  cmean_warm[obs_map_fit[warm,i]==0]<-0
  cmean_cold[cmean_cold<0]<-0
  cmean_warm[cmean_warm<0]<-0
  prcp_cmean[cold,i]<-cmean_cold
  prcp_cmean[warm,i]<-cmean_warm
  prcp_err_cmean[,(i+1)]<-prcp_cmean[,i]-prcp_forc[,(i+1)]
}

#par(mfrow=c(4,4))

for(i in 1:64){
  loess_fit_cold<-loess(prcp_forc[cold,(i+1)]~obs_map_fit[cold,i],span=1,degree = 1,family = 'symmetric',
        control=loess.control(surface='interpolate'))
  loess_fit_warm<-loess(prcp_forc[warm,(i+1)]~obs_map_fit[warm,i],span=1,degree = 1,family = 'gaussian',
                        control=loess.control(surface='interpolate'))
  ls_cold<-predict(loess_fit_cold,obs_map_fit[cold,i])
  ls_warm<-predict(loess_fit_warm,obs_map_fit[warm,i])
  ls_cold[ls_cold<0]<-0
  ls_warm[ls_warm<0]<-0
  prcp_ls[cold,i]<-ls_cold
  prcp_ls[warm,i]<-ls_warm
  prcp_err_ls[,(i+1)]<-prcp_ls[,i]-prcp_forc[,(i+1)]
  #plot(obs_map_fit[cold,i],prcp_forc[cold,(i+1)])
  #lines(obs_map_fit[cold,i],cmean_cold,col='red')
}

#i<-64
#plot(obs_map_fit[cold,i],prcp_forc[cold,(i+1)],xlim=c(0,10),ylim=c(0,10))
#points(obs_map_fit[cold,i],prcp_cmean[cold,i],col='red')


#debias
tmax_mean<-array(NA,c(2,64))
tmin_mean<-array(NA,c(2,64))

tmax_mean[1,]<-apply(tmax_err[cold,2:65],2,mean)
tmax_err_db<-tmax_err
tmax_err_db[cold,2:65]<-tmax_err[cold,2:65] - matrix(rep(tmax_mean[1,],length(cold)),ncol=64,byrow=T)

tmax_mean[2,]<-apply(tmax_err[warm,2:65],2,mean)
tmax_err_db[warm,2:65]<-tmax_err[warm,2:65] - matrix(rep(tmax_mean[2,],length(warm)),ncol=64,byrow=T)

apply(tmax_err_db[,2:65],2,mean)

tmin_mean[1,]<-apply(tmin_err[cold,2:65],2,mean)
tmin_err_db<-tmin_err
tmin_err_db[cold,2:65]<-tmin_err[cold,2:65] - matrix(rep(tmin_mean[1,],length(cold)),ncol=64,byrow=T)

tmin_mean[2,]<-apply(tmin_err[warm,2:65],2,mean)
tmin_err_db[warm,2:65]<-tmin_err[warm,2:65] - matrix(rep(tmin_mean[2,],length(warm)),ncol=64,byrow=T)

apply(tmin_err_db[,2:65],2,mean)


#3d. Create precip occurrence array (day of, day prior across 3 locations)
map_hop_fit<-map_hop_mat[which(ix=='1985-01-01'):which(ix=='2017-09-30'),]
map_hop_prior<-map_hop_mat[which(ix=='1984-12-31'):which(ix=='2017-09-29'),]
map_lam_fit<-map_lam_mat[which(ix=='1985-01-01'):which(ix=='2017-09-30'),]
map_lam_prior<-map_lam_mat[which(ix=='1984-12-31'):which(ix=='2017-09-29'),]
map_uka_fit<-map_uka_mat[which(ix=='1985-01-01'):which(ix=='2017-09-30'),]
map_uka_prior<-map_uka_mat[which(ix=='1984-12-31'):which(ix=='2017-09-29'),]

prcp_occ_arr<-cbind(map_hop_fit,map_hop_prior,map_lam_fit,map_lam_prior,map_uka_fit,map_uka_prior)
prcp_occ_arr[prcp_occ_arr>0]<-1
prcp_samp_arr<-cbind(map_hop_fit,map_hop_prior,map_lam_fit,map_lam_prior,map_uka_fit,map_uka_prior)

#-------------------------------------------------------------------------------------
#4) Save data
saveRDS(ob_mat,'data/basic/ob_mat.rds')
saveRDS(ob_map,'data/basic/ob_map.rds')
saveRDS(obs_mat_fit,'data/basic/obs_mat_fit.rds')
saveRDS(obs_map_fit,'data/basic/obs_map_fit.rds')

saveRDS(gefs_prcp,'data/basic/gefs_prcp.rds')
saveRDS(gefs_prcp_fit,'data/basic/gefs_prcp_fit.rds')
saveRDS(gefs_tmax,'data/basic/gefs_tmax.rds')
saveRDS(gefs_tmax_fit,'data/basic/gefs_tmax_fit.rds')
saveRDS(gefs_tmin,'data/basic/gefs_tmin.rds')
saveRDS(gefs_tmin_fit,'data/basic/gefs_tmin_fit.rds')

saveRDS(tmax_forc,'data/basic/tmax_forc.rds')
saveRDS(tmin_forc,'data/basic/tmin_forc.rds')
saveRDS(prcp_forc,'data/basic/prcp_forc.rds')

saveRDS(prcp_err,'data/basic/prcp_resids.rds')
saveRDS(tmax_err,'data/basic/tmax_resids.rds')
saveRDS(tmin_err,'data/basic/tmin_resids.rds')
saveRDS(tmax_err_db,'data/basic/tmax_resids_ub.rds')
saveRDS(tmin_err_db,'data/basic/tmin_resids_ub.rds')

saveRDS(prcp_occ_arr,'data/basic/prcp_occ_arr.rds')
saveRDS(prcp_samp_arr,'data/basic/prcp_samp_arr.rds')
saveRDS(tmax_mean,'data/basic/mn_tmax.rds')
saveRDS(tmin_mean,'data/basic/mn_tmin.rds')

saveRDS(prcp_cmean,'data/basic/prcp_cmean.rds')
saveRDS(prcp_err_cmean,'data/basic/prcp_cmean_resids.rds')

saveRDS(prcp_ls,'data/basic/prcp_loess.rds')
saveRDS(prcp_err_ls,'data/basic/prcp_loess_resids.rds')
#-----------------------------------------------------------------
#5) Historical and Climatological Data Processing
#Processes climatological averages used for skill score calculations
#Same basic procedure as above


#--------------------------------------------------------------------------

#6) Calculate Climatology
ix<-seq(as.Date('1948-10-02'),as.Date('2017-09-30'),'day')
ixx<-as.POSIXlt(ix)
cold<-sort(c(which(ixx$mon>8),which(ixx$mon<3)))
warm<-1:length(ixx)
warm<-warm[-c(cold)]


#Calculate climatology for each day of year
ob_mat<-readRDS('data/basic/ob_mat.rds')
ob_map<-readRDS('data/basic/ob_map.rds')
tmax_mean<-readRDS('data/basic/mn_tmax.rds')
tmin_mean<-readRDS('data/basic/mn_tmin.rds')

tmax_mat_cold<-matrix(rep(tmax_mean[1,1:4],length(cold)),ncol=4,byrow=T)
tmax_mat_warm<-matrix(rep(tmax_mean[2,1:4],length(warm)),ncol=4,byrow=T)
tmin_mat_cold<-matrix(rep(tmin_mean[1,1:4],length(cold)),ncol=4,byrow=T)
tmin_mat_warm<-matrix(rep(tmin_mean[2,1:4],length(warm)),ncol=4,byrow=T)

ob_tmin<-ob_mat
ob_tmax<-ob_mat

ob_tmin[cold,]<-ob_mat[cold,]-matrix(rep(tmin_mat_cold,16),ncol=64,byrow=F)
ob_tmin[warm,]<-ob_mat[warm,]-matrix(rep(tmin_mat_warm,16),ncol=64,byrow=F)
ob_tmax[cold,]<-ob_mat[cold,]-matrix(rep(tmax_mat_cold,16),ncol=64,byrow=F)
ob_tmax[warm,]<-ob_mat[warm,]-matrix(rep(tmax_mat_warm,16),ncol=64,byrow=F)

saveRDS(ob_tmin,'ob_tmin.rds')
saveRDS(ob_tmin,'ob_tmax.rds')

hist_climo<-array(NA,c(3,length(ix),64))

#for(i in 1:length(ixx)){
  #hist_climo[1,i,]<-apply(ob_map[which(ixx$mon==ixx[i]$mon & ixx$mday==ixx[i]$mday),],2,mean)
  #hist_climo[2,i,]<-apply(ob_tmax[which(ixx$mon==ixx[i]$mon & ixx$mday==ixx[i]$mday),],2,mean)
  #hist_climo[3,i,]<-apply(ob_tmin[which(ixx$mon==ixx[i]$mon & ixx$mday==ixx[i]$mday),],2,mean)
#}

for(i in 1:length(ixx)){
  hist_climo[1,i,]<-apply(ob_map[which(ixx$mon==ixx[i]$mon & ixx$mday==ixx[i]$mday),],2,mean)
  hist_climo[2,i,]<-apply(ob_mat[which(ixx$mon==ixx[i]$mon & ixx$mday==ixx[i]$mday),],2,mean)
  hist_climo[3,i,]<-apply(ob_mat[which(ixx$mon==ixx[i]$mon & ixx$mday==ixx[i]$mday),],2,mean)
}

gefs_idx<-which(ix=='1985-01-01'):which(ix=='2017-09-30')
climo<-hist_climo[,gefs_idx,]

ob_tmin_fit<-ob_tmin[gefs_idx,]
ob_tmax_fit<-ob_tmax[gefs_idx,]

saveRDS(ob_tmin_fit,'data/basic/ob_tmin_fit.rds')
saveRDS(ob_tmax_fit,'data/basic/ob_tmax_fit.rds')

saveRDS(hist_climo,'data/basic/hist_climo.rds')
saveRDS(climo,'data/basic/climo.rds')

rm(list=ls());gc()

#################################END#############################