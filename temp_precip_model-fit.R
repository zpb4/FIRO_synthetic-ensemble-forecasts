###################################################################################
# Code to fit Multivariate Generalized Likelihood/Skew Generalized Error Distribution
# (GL/SGED) synthetic forecast model to TMAX, TMIN, and PRECIP data
# 
# Supports manuscript: 'A generalized, multivariate approach to generate synthetic 
# short-to-medium range hydro-meteorological forecasts across locations, variables, 
# and lead times', Authors - Brodeur, Z, and Steinschneider, S.
# 
# Created by Zach Brodeur, 15 December 2020
# -Updated 18 April 2021: Standardized and clarified for replicability
##################################################################################

#required packages and drive
setwd('h:/firo_lamc/temp_precip/')
library(fGarch)
library(logisticPCA)
library(BigVAR)


#define cold (ONDJFM) and warm (AMJJAS) seasons
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]

src<-'data/mag-samp_loess'
#------------------------------------------------------------------------------
#1) kNN sample based on day of/day prior precipitation occurrence array defined in 
#temp_precip_data-process. This will be used later in sampling the a_t indices, but
#can be done at any point in the process since the precip occurrence arrays have
#already been defined

#1a. Fit logistic PCA model to precipitation occurrence array (day of, day prior)
prcp_samp_arr<-readRDS('data/basic/prcp_samp_arr.rds')
pc<-prcomp(prcp_samp_arr,center=F,scale.=F)
eig_val<-pc$sdev^2
var_expl<-c()
for(i in 1:length(eig_val)){
  var_expl[i]<-sum(eig_val[1:i])/sum(eig_val)*100
}
plot(1:length(eig_val),var_expl,ylim=c(0,100),type='l')
abline(v=7,col='red')

pc_fit<-pc$x[,1:8]
pc_samp<-prcp_samp_arr%*%pc$rotation
pc_samp<-pc_samp[,1:8]
pc_size=2

#obs_agg<-apply(obprcp_arr[,4:7],1,sum)

#plot(pc$x[,4],obs_agg)
#KNN sampling params
rng<-15 #+/- range from specified day in calendar year to sample from, so if Jan-15
#for example, will only look for matches in Jan-1 to Jan-30 across all years.
knn<-round(sqrt(length(cold))) #square root of n

#for memory considerations, only 100 samples at a time possible; saved 10 times over gets to 1000
n<-1
m<-10

#weighted kernel sampling as for streamflow
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn)

#1b. kNN sampling routine to create (n x m) # of time series of kNN sampled indices based on
#'likeness' between precip occurrence observation at time t against one in the 'fitted' period
knn_idx<-array(NA,c(length(ix2),m))

for(p in 1:n){
  for(k in 1:m){
    #define +/- range subset to sample from
    for(i in 1:length(ix2)){
      #watch for leap years!
      if(ixx2[i]$mon==1 & ixx2[i]$mday==29) ctr<-which(ixx2$mon==ixx2[i]$mon & ixx2$mday==ixx2[i-1]$mday) else
        ctr<-which(ixx2$mon==ixx2[i]$mon & ixx2$mday==ixx2[i]$mday)
      lb<-ctr-rng; lb[lb<1] <- 1
      ub<-ctr+rng; ub[ub>length(ix2)] <- length(ix2)
      sset<-c()
      for(j in 1:length(lb)){
        sset<-c(sset,lb[j]:ub[j])
      }
      #calculate PC vector for current time index
      if(pc_size==1){
      pc_vec<-pc_samp[i]
      pc_mat<-rep(pc_vec,length(sset))
      sub<-pc_fit[sset]
      dist<-(sub-pc_mat)
      } 
      
      else{pc_vec<-pc_samp[i,]
      #repeat to length of time subset for simple matrix operations
      pc_mat<-matrix(rep(pc_vec,length(sset)),ncol=length(pc_vec),byrow=T)
      #slice out all the PC vectors for the subset
      sub<-pc_fit[sset,]
      #calculate Euclidean distance between sample and all examples in subset
      dist<-sqrt(apply((sub-pc_mat)^2,1,sum))}
      #random sample from 0 distance examples if more than one
      if(length(dist[dist==0])>1) {
        id<-sample(which(dist==0),1);
        fidx<-sset[id];
        knn_idx[i,k]<-fidx} else {
        #otherwise weighted sample for all distances > 0
        x<-sort(dist)
        x<-x[1:knn] #top 30 values are KNN
        s<-sample(x,1,prob=wts)
        id<-which(dist==s)
        if(length(id)>1) {id <- sample(id,1)} #resample for duplicate values
        fidx<-sset[id]
        knn_idx[i,k]<-fidx}
    }
    print(paste(k,Sys.time()))
  }
 saveRDS(knn_idx,paste(src,'/knn_idx_',p,'.rds',sep=''))
}

print(paste('S1 complete',Sys.time()))
#-------------------------------------------------------------------
#2) BigVar model
#VAR model to create decorrelated residuals
#follows implementation in fol_fit-model

ar<-3

uc_resids<-array(NA,c(3,length(ix2),64))
var_coefs_prcp<-array(NA,c(2,64,64*3))

prcp_resids<-readRDS('data/basic/prcp_loess_resids.rds')
tmax_resids<-readRDS('data/basic/tmax_resids_ub.rds')
tmin_resids<-readRDS('data/basic/tmin_resids_ub.rds')

#Calculate uncorrelated matrices
sn<-list(cold,warm,cold,warm,cold,warm)
rresids<-list(prcp_resids,prcp_resids,tmax_resids,tmax_resids,tmin_resids,tmin_resids)
idx<-c(1,1,2,2,3,3)

#6 indices for precip-cold, precip-warm, tmax-cold,...)
for(i in 1:2){
  rresid_mat<-rresids[[i]][sn[[i]],2:65]
  rresid_mat<-as.matrix(rresid_mat)
  #bvar_fit<-BigVAR.fit(rresid_mat,p=ar,struct = "Basic",lambda = 10,intercept = F)
  #bvar<-bvar_fit[,2:193,1]
  #lagmatrix<-cbind(rresid_mat[(ar):(length(rresid_mat[,1])-1),],rresid_mat[(ar-1):(length(rresid_mat[,1])-2),],rresid_mat[(ar-2):(length(rresid_mat[,1])-3),])
  #var_coefs_prcp[i,,]<-bvar
  #resids<-lagmatrix %*% t(bvar) #calculate decorrelation residuals from bvar_fit
  #mat<-matrix(0,dim(rresid_mat)[1],dim(rresid_mat)[2])
  #mat[(ar+1):dim(mat)[1],]<-resids #resids set to 0 for first columns = max lag
  #uc_mat<-rresid_mat - mat
  #uc_resids[idx[i],sn[[i]],]<-uc_mat
  
  m1 = constructModel(rresid_mat, p = ar, struct = "Basic", gran = c(25, 10),IC = F,
                      verbose = T, VARX = list(),MN = F, separate_lambdas = F,intercept = F)
  m1_res = cv.BigVAR(m1)
  var_coefs_prcp[i,,]<-m1_res@betaPred[,2:length(m1_res@betaPred[1,])]
  resids<-m1_res@fitted #calculate decorr residuals from 
  mat<-matrix(0,dim(rresid_mat)[1],dim(rresid_mat)[2])
  mat[(ar+1):dim(mat)[1],]<-resids #resids set to 0 for first columns = max lag
  uc_mat<-rresid_mat - mat
  uc_resids[idx[i],sn[[i]],]<-uc_mat
}

saveRDS(var_coefs_prcp,paste(src,'/var_coefs_prcp.rds',sep=''))

var_coefs_temp<-array(NA,c(2,128,128*3))
rr_idx<-rbind(c(3,5),c(4,6))
idx<-rbind(c(2,3),c(2,3))

for(i in 1:2){
  rresid_mat<-cbind(rresids[[rr_idx[i,1]]][sn[[i]],2:65],rresids[[rr_idx[i,2]]][sn[[i]],2:65])
  rresid_mat<-as.matrix(rresid_mat)
  #bvar_fit<-BigVAR.fit(rresid_mat,p=ar,struct = "Basic",lambda = 1,intercept = F)
  #bvar<-bvar_fit[,2:385,1]
  #lagmatrix<-cbind(rresid_mat[(ar):(length(rresid_mat[,1])-1),],rresid_mat[(ar-1):(length(rresid_mat[,1])-2),],rresid_mat[(ar-2):(length(rresid_mat[,1])-3),])
  #var_coefs_temp[i,,]<-bvar
  #resids<-lagmatrix %*% t(bvar) #calculate decorrelation residuals from bvar_fit
  #mat<-matrix(0,dim(rresid_mat)[1],dim(rresid_mat)[2])
  #mat[(ar+1):dim(mat)[1],]<-resids #resids set to 0 for first columns = max lag
  #uc_mat<-rresid_mat - mat
  
  m1 = constructModel(rresid_mat, p = ar, struct = "Basic", gran = c(25, 10),IC = F,
                      verbose = T, VARX = list(),MN = F, separate_lambdas = F,intercept = F)
  m1_res = cv.BigVAR(m1)
  var_coefs_temp[i,,]<-m1_res@betaPred[,2:length(m1_res@betaPred[1,])]
  resids<-m1_res@fitted #calculate decorr residuals from 
  mat<-matrix(0,dim(rresid_mat)[1],dim(rresid_mat)[2])
  mat[(ar+1):dim(mat)[1],]<-resids #resids set to 0 for first columns = max lag
  uc_mat<-rresid_mat - mat
  
  uc_resids[idx[i,1],sn[[i]],]<-uc_mat[,1:64]
  uc_resids[idx[i,2],sn[[i]],]<-uc_mat[,65:128]
}

saveRDS(var_coefs_temp,paste(src,'/var_coefs_temp.rds',sep=''))
saveRDS(uc_resids,paste(src,'/uc_resids.rds',sep=''))
rm(prcp_resids,tmax_resids,tmin_resids,uc_resids)

print(paste('S2 complete',Sys.time()))
#-----------------------------------------------------------------------------------------------------
#3) GL SGED Model
#monthly model fits for each lead time and variable
#follows implementation in fol_fit-model

obprcp_arr<-readRDS('data/basic/obs_map_fit.rds')
prcp_loess<-readRDS('data/basic/prcp_loess.rds')
fprcp_arr<-readRDS('data/basic/prcp_forc.rds')
obtmax_arr<-readRDS('data/basic/obs_mat_fit.rds')
obtmax_arr<-obtmax_arr+273
obtmin_arr<-readRDS('data/basic/obs_mat_fit.rds')
obtmin_arr<-obtmin_arr+273
uc_resids<-readRDS(paste(src,'/uc_resids.rds',sep=''))

gl_par_arr<-array(NA,c(6,64,4))
at_arr<-array(0,c(3,length(ix2),64))
at_nep_arr<-array(0,c(3,length(ix2),64))
sn<-list(cold,warm)

source('GL_maineqs.R')

#1 is cold season, 2 is warm season fits
for(k in 1:2){
  for(i in 1:64){
    uc_sim<-uc_resids[1,sn[[k]],i]
    fp<-obprcp_arr[sn[[k]],i]
    fp_cm<-prcp_loess[sn[[k]],i]
    fp2<-fprcp_arr[sn[[k]],i+1]
    uc_sim<-uc_sim[fp2>0] # only fit GL_SGED to precip forecasts > 0
    sim<-fp_cm[fp2>0]
    gl_mle<-optim(par=c(.5,.5,0,1),fn=GL_fun_noscale_obsxxx,inflow=sim,sim_inflow=sim,et=uc_sim,
                  method = 'L-BFGS-B',lower = c(0.01,0,-0.99,0.1),upper = c(3,1,1,10),
                  control = list(maxit=100000))
    
    at<-a_txx(uc_sim,sigma_t(gl_mle$par[1],gl_mle$par[2],sim))
    gl_par_arr[k,i,]<-c(gl_mle$par[1],gl_mle$par[2],gl_mle$par[3],gl_mle$par[4])
    at_ins<-at_arr[1,sn[[k]],i]
    at_ins[fp2>0]<-at #a_t calculated only for non-zero forecasts
    at_arr[1,sn[[k]],i]<-at_ins
    #calculate NEPs from fitted GL_SGED parameters
    at_ins[fp2>0]<-psged(at,mean=0,sd=1,nu=(2/(1 + gl_mle$par[3])),xi=gl_mle$par[4])
    at_nep_arr[1,sn[[k]],i]<-at_ins
    
    uc_sim<-uc_resids[2,sn[[k]],i]
    sim<-obtmax_arr[sn[[k]],i] / 100 #decrease scaling of temperature in K
    gl_mle<-optim(par=c(.5,.5,0,1),fn=GL_fun_noscale_obsxxx,inflow=sim,sim_inflow=sim,et=uc_sim,
                  method = 'L-BFGS-B',lower = c(0.01,0,-0.99,0.1),upper = c(3,1,1,10),
                  control = list(maxit=100000))
    
    at<-a_txx(uc_sim,sigma_t(gl_mle$par[1],gl_mle$par[2],sim))
    gl_par_arr[k+2,i,]<-c(gl_mle$par[1],gl_mle$par[2],gl_mle$par[3],gl_mle$par[4])
    at_arr[2,sn[[k]],i]<-at
    #calculate NEPs from fitted GL_SGED parameters
    at_nep_arr[2,sn[[k]],i]<-psged(at,mean=0,sd=1,nu=(2/(1 + gl_mle$par[3])),xi=gl_mle$par[4])
    
    uc_sim<-uc_resids[3,sn[[k]],i]
    sim<-obtmin_arr[sn[[k]],i] / 100#decrease scaling of temperature in K
    gl_mle<-optim(par=c(.5,.5,0,1),fn=GL_fun_noscale_obsxxx,inflow=sim,sim_inflow=sim,et=uc_sim,
                  method = 'L-BFGS-B',lower = c(0.01,0,-0.99,0.1),upper = c(3,1,1,10),
                  control = list(maxit=100000))
    
    at<-a_txx(uc_sim,sigma_t(gl_mle$par[1],gl_mle$par[2],sim))
    gl_par_arr[k+4,i,]<-c(gl_mle$par[1],gl_mle$par[2],gl_mle$par[3],gl_mle$par[4])
    at_arr[3,sn[[k]],i]<-at
    #calculate NEPs from fitted GL_SGED parameters
    at_nep_arr[3,sn[[k]],i]<-psged(at,mean=0,sd=1,nu=(2/(1 + gl_mle$par[3])),xi=gl_mle$par[4])
  }
}

saveRDS(gl_par_arr,paste(src,'/gl_par_arr.rds',sep=''))
saveRDS(at_arr,paste(src,'/at_arr.rds',sep=''))
saveRDS(at_nep_arr,paste(src,'/at_nep_arr.rds',sep=''))

#rm(list=ls());gc()

print(paste('S3 complete',Sys.time()))
#------------------------------------------------------------------
#4) Precip Occurrence for kNN sampled arrays
#assigns occurrence/non-occurrence in forecasts based on sampled index in kNN
#In other words, imposes forecast positive/negative structure from sampled days of the
#'fitted' record

#load forecast data for PRECIP
fprcp_arr<-readRDS('prcp_forc.rds')

#fitted period date vector
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)

#arrays of 1/0 for occurrence/non-occurrence of precip in actual forecasts
fprcp_occ_arr<-fprcp_arr[,2:65]
fprcp_occ_arr[fprcp_occ_arr>0]<-1

#define occurrence/non-occurrence in forecasts for each kNN sampled index
#n<-1
#m<-100
knn_prcp_occ<-array(0,c(m,length(ix2),64))

for(s in 1:n){
  knn_idx<-readRDS(paste(src,'/knn_idx_',s,'.rds',sep=''))
  for(i in 1:m){
    #each index in kNN sampled array associated with occurrence/non-occurrences in forecasts
    knn_prcp_occ[i,,]<-as.matrix(fprcp_occ_arr[knn_idx[,i],])
  }
  print(paste(s,Sys.time()))
  saveRDS(knn_prcp_occ,paste(src,'/knn_prcp_occ_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

print(paste('S4 complete',Sys.time()))
#------------------------------------------------------------------
#5) Create a_t 'sampling' arrays from matrices of a_t values for each variable,
#location, lead time. a_t values are resample according to kNN indices from part 1

#5a. PRECIP

#load fitted a_t for precip
at_arr<-readRDS(paste(src,'/at_arr.rds',sep=''))

#date vector for fitted period
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)

#create 100 samples (100 samples 10 times over for data processing constraints)
#n<-1
#m<-10

#assign empirical a_t entries for each kNN index
at_samp_arr_prcp<-array(0,c(m,length(ix2),64))

for(s in 1:n){
  #kNN sampled indices
  knn_idx<-readRDS(paste(src,'/knn_idx_',s,'.rds',sep=''))
  for(i in 1:m){
    #assign a_t value to each sample index
    at_samp_arr_prcp[i,,]<-at_arr[1,knn_idx[,i],]
  }
  print(paste(s,Sys.time()))
  saveRDS(at_samp_arr_prcp,paste(src,'/at_samp_arr_prcp_',s,'.rds',sep=''))
}

rm(at_samp_arr_prcp)

#5b. TMAX
#assign empirical a_t entries for each kNN index
at_samp_arr_tmax<-array(0,c(m,length(ix2),64))

for(s in 1:n){
  knn_idx<-readRDS(paste(src,'/knn_idx_',s,'.rds',sep=''))
  for(i in 1:m){
    at_samp_arr_tmax[i,,]<-at_arr[2,knn_idx[,i],]
  }
  print(paste(s,Sys.time()))
  saveRDS(at_samp_arr_tmax,paste(src,'/at_samp_arr_tmax_',s,'.rds',sep=''))
}

rm(at_samp_arr_tmax)

#5c. TMIN
#assign empirical a_t entries for each kNN index
at_samp_arr_tmin<-array(0,c(m,length(ix2),64))

for(s in 1:n){
  knn_idx<-readRDS(paste(src,'/knn_idx_',s,'.rds',sep=''))
  for(i in 1:m){
    at_samp_arr_tmin[i,,]<-at_arr[3,knn_idx[,i],]
  }
  saveRDS(at_samp_arr_tmin,paste(src,'/at_samp_arr_tmin_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

print(paste('S5 complete',Sys.time()))
#----------------------------------------------------------------------------------------------------------------------
#6) Create new arrays of a_t by randomly sampling from fitted GL_SGED for each location, variable, lead time,
#and rearranging them according to rank structure in kNN sampled a_t arrays from step 5

#6a. PRECIP
#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
#n<-1
#m<-10

#load GL_SGED parameters
gl_par_arr<-readRDS(paste(src,'/gl_par_arr.rds',sep=''))

#For each timeseries of kNN indices, generate a_t from fitted GL_SGED distribution and
#then Schaake shuffle to match empirical a_t rank structure
for(s in 1:n){
  at_samp_arr_prcp<-readRDS(paste(src,'/at_samp_arr_prcp_',s,'.rds',sep=''))
  knn_at_samp_arr_prcp<-array(0,c(m,length(ix2),64))
  for(i in 1:m){
    for(k in 1:64){
      #cold season section
      at_prcp1<-at_samp_arr_prcp[i,cold,k]
      at_prcp<-at_prcp1[at_prcp1!=0] #rank only non-zero values
      at_prcp_ranks<-rank(at_prcp,ties.method = 'random')
      #generate new randomly sampled a_t
      new_at<-rsged(length(at_prcp),mean=0,sd=1,nu=(2/(1+gl_par_arr[1,k,3])),xi=gl_par_arr[1,k,4])
      #define ranks of new a_t
      new_at_ranks<-rank(new_at,ties.method = 'random')
      at_repl<-c()
      #rearrange new a_t to match rank structure in empirical data
      for(j in 1:length(new_at)){
        at_repl[j]<-new_at[which(new_at_ranks==at_prcp_ranks[j])]
      }
      knn_ins<-knn_at_samp_arr_prcp[i,cold,k]
      #replace non-zero data in array with resampled/shuffled a_t
      knn_ins[at_prcp1!=0]<-at_repl
      knn_at_samp_arr_prcp[i,cold,k]<-knn_ins
    
      #warm season section, same process as above
      at_prcp1<-at_samp_arr_prcp[i,warm,k]
      at_prcp<-at_prcp1[at_prcp1!=0]
      at_prcp_ranks<-rank(at_prcp,ties.method = 'random')
      new_at<-rsged(length(at_prcp),mean=0,sd=1,nu=(2/(1+gl_par_arr[2,k,3])),xi=gl_par_arr[2,k,4])
      new_at_ranks<-rank(new_at,ties.method = 'random')
      at_repl<-c()
      for(j in 1:length(new_at)){
        at_repl[j]<-new_at[which(new_at_ranks==at_prcp_ranks[j])]
      }
      knn_ins<-knn_at_samp_arr_prcp[i,warm,k]
      knn_ins[at_prcp1!=0]<-at_repl
      knn_at_samp_arr_prcp[i,warm,k]<-knn_ins
    }
  }
  print(paste(s,Sys.time()))
  saveRDS(knn_at_samp_arr_prcp,paste(src,'/knn_at_samp_arr_prcp_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

#6b. TMAX
#TMAX follows precip, except no requirement to delineate non-zero a_t

#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
#n<-1
#m<-10

#load GL_SGED parameters
gl_par_arr<-readRDS(paste(src,'/gl_par_arr.rds',sep=''))

for(s in 1:n){
  at_samp_arr_tmax<-readRDS(paste(src,'/at_samp_arr_tmax_',s,'.rds',sep=''))
  knn_at_samp_arr_tmax<-array(0,c(m,length(ix2),64))
  for(i in 1:m){
    for(k in 1:64){
      #cold season section
      at_tmax<-at_samp_arr_tmax[i,cold,k]
      at_tmax_ranks<-rank(at_tmax,ties.method = 'random')
      new_at<-rsged(length(at_tmax),mean=0,sd=1,nu=(2/(1+gl_par_arr[3,k,3])),xi=gl_par_arr[3,k,4])
      new_at_ranks<-rank(new_at,ties.method = 'random')
      at_repl<-c()
      for(j in 1:length(new_at)){
        at_repl[j]<-new_at[which(new_at_ranks==at_tmax_ranks[j])]
      }
      knn_at_samp_arr_tmax[i,cold,k]<-at_repl
    
      #warm season section
      at_tmax<-at_samp_arr_tmax[i,warm,k]
      at_tmax_ranks<-rank(at_tmax,ties.method = 'random')
      new_at<-rsged(length(at_tmax),mean=0,sd=1,nu=(2/(1+gl_par_arr[4,k,3])),xi=gl_par_arr[4,k,4])
      new_at_ranks<-rank(new_at,ties.method = 'random')
      at_repl<-c()
      for(j in 1:length(new_at)){
        at_repl[j]<-new_at[which(new_at_ranks==at_tmax_ranks[j])]
      }
      knn_at_samp_arr_tmax[i,warm,k]<-at_repl
    }
  }
  print(paste(s,Sys.time()))
  saveRDS(knn_at_samp_arr_tmax,paste(src,'/knn_at_samp_arr_tmax_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

#6c. TMIN
#TMIN follows TMAX

#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
#n<-1
#m<-10

#load GL_SGED parameters
gl_par_arr<-readRDS(paste(src,'/gl_par_arr.rds',sep=''))

for(s in 1:n){
  at_samp_arr_tmin<-readRDS(paste(src,'/at_samp_arr_tmin_',s,'.rds',sep=''))
  knn_at_samp_arr_tmin<-array(0,c(m,length(ix2),64))
  for(i in 1:m){
    for(k in 1:64){
      at_tmin<-at_samp_arr_tmin[i,cold,k]
      at_tmin_ranks<-rank(at_tmin,ties.method = 'random')
      new_at<-rsged(length(at_tmin),mean=0,sd=1,nu=(2/(1+gl_par_arr[5,k,3])),xi=gl_par_arr[5,k,4])
      new_at_ranks<-rank(new_at,ties.method = 'random')
      at_repl<-c()
      for(j in 1:length(new_at)){
        at_repl[j]<-new_at[which(new_at_ranks==at_tmin_ranks[j])]
      }
      knn_at_samp_arr_tmin[i,cold,k]<-at_repl
    
      at_tmin<-at_samp_arr_tmin[i,warm,k]
      at_tmin_ranks<-rank(at_tmin,ties.method = 'random')
      new_at<-rsged(length(at_tmin),mean=0,sd=1,nu=(2/(1+gl_par_arr[6,k,3])),xi=gl_par_arr[6,k,4])
      new_at_ranks<-rank(new_at,ties.method = 'random')
      at_repl<-c()
      for(j in 1:length(new_at)){
        at_repl[j]<-new_at[which(new_at_ranks==at_tmin_ranks[j])]
      }
      knn_at_samp_arr_tmin[i,warm,k]<-at_repl
    }
  }
  print(paste(s,Sys.time()))
  saveRDS(knn_at_samp_arr_tmin,paste(src,'/knn_at_samp_arr_tmin_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

print(paste('S6 complete',Sys.time()))
#---------------------------------------------------------------------------------------------------------------------
#7) Calculate NEP for each kNN a_t sample array from step 6 using fitted GL_SGED

gl_par_arr<-readRDS(paste(src,'/gl_par_arr.rds',sep=''))
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
#n<-1
#m<-10

knn_at_nep_arr<-array(0,c(m,3,length(ix2),64))
sn<-list(cold,warm,cold,warm,cold,warm)
idx<-c(1,1,2,2,3,3)

for(s in 1:n){
  knns<-c(rep(paste(src,'/knn_at_samp_arr_prcp_',s,'.rds',sep=''),2),
          rep(paste(src,'/knn_at_samp_arr_tmax_',s,'.rds',sep=''),2),
          rep(paste(src,'/knn_at_samp_arr_tmin_',s,'.rds',sep=''),2))
  for(i in 1:m){
    for(j in 1:6){
      knn_at_samp_arr<-readRDS(knns[j])
      for(k in 1:64){
        at<-knn_at_samp_arr[i,sn[[j]],k]
        at1<-at[at!=0]
        nep_at<-psged(at1,mean=0,sd=1,nu=(2/(1+gl_par_arr[j,k,3])),xi=gl_par_arr[j,k,4])
        ins<-rep(0,length(at))
        ins[at!=0]<-nep_at
        knn_at_nep_arr[i,idx[j],sn[[j]],k]<-ins
      }
    }
  }
  print(paste(s,Sys.time()))
  saveRDS(knn_at_nep_arr,paste(src,'/knn_at_nep_arr_',s,'.rds',sep=''))
}


rm(list=ls());gc()

#################################################END####################################################################