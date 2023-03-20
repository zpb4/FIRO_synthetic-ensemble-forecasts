###################################################################################
# Code to synthetically generate multiple forecasts of TMAX, TMIN, and PRECIP
# --Synthetic forecasts in fitted period for reliability/skill analyses
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

#-----------------------------------------------------------------------------------------
#1) Create synthetic forecast errors for each variable

#1a. PRECIP
#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
sn<-list(cold,warm,cold,warm,cold,warm)
idx<-c(1,1,2,2,3,3)
ar<-3

n<-1
m<-10

src<-'data/mag-samp_loess'
out<-'out/mag-samp_loess'
#load data
syn_res_prcp<-array(0,c(m,length(ix2),64))
syn_resids_prcp<-array(0,c(m,length(ix2),64))
source('GL_maineqs.R')
obprcp_arr<-readRDS('data/basic/obs_map_fit.rds')
prcp_loess<-readRDS('data/basic/prcp_loess.rds')
var_coefs_prcp<-readRDS(paste(src,'/var_coefs_prcp.rds',sep=''))
gl_par_arr<-readRDS(paste(src,'/gl_par_arr.rds',sep=''))

#create (n x m) no. of synthetic error matrices
for(s in 1:n){
  #load Schaake shuffled and kNN sampled a_t
  knn_at_samp_arr_prcp<-readRDS(paste(src,'/knn_at_samp_arr_prcp_',s,'.rds',sep=''))
  for(d in 1:m){
    for(l in 1:2){
      mat<-knn_at_samp_arr_prcp[d,sn[[l]],]
      coeff<-var_coefs_prcp[l,,]
      ob<-obprcp_arr[sn[[l]],]
      ob_cm<-prcp_loess[sn[[l]],]
      syn_resid_mat<-matrix(0,ncol=64,nrow=dim(mat)[1])
      syn_resid_mat2<-matrix(0,ncol=64,nrow=dim(mat)[1])
      #generate forecast errors across all dimensions
      for(j in 1:64){
        for(k in (ar+1):dim(mat)[1]){
          syn_resid_mat[k,j]<-t(matrix(c(syn_resid_mat[(k-1),],syn_resid_mat[(k-2),],syn_resid_mat[(k-3),]))) %*% matrix(coeff[j,]) + 
            sigma_t(gl_par_arr[l,j,1],gl_par_arr[l,j,2],ob_cm[k,j])*mat[k,j]
          syn_res2<-syn_resid_mat[k,j]
          #resample from a_t distribution if error would produce negative forecast
          #'corrected' errors
          while(syn_res2>ob_cm[k,j]){
            ats<-rsged(1,mean=0,sd=1,nu=(2/(1+gl_par_arr[l,j,3])),xi=gl_par_arr[l,j,4])
            syn_res2<-sigma_t(gl_par_arr[l,j,1],gl_par_arr[l,j,2],ob_cm[k,j]) * ats
          }
          syn_resid_mat2[k,j]<-syn_res2
        }
      }
      syn_res_prcp[d,sn[[l]],]<-syn_resid_mat #uncorrected errors
      syn_resids_prcp[d,sn[[l]],]<-syn_resid_mat2 #save corrected errors
    }
    print(paste(d,Sys.time()))
  }
  saveRDS(syn_resids_prcp,paste(out,'/syn_resids_prcp_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

print(paste('S1a complete',Sys.time()))

#1b. TEMP
#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
sn<-list(cold,warm,cold,warm,cold,warm)
sn2<-c(1,2,1,2,1,2)
idx<-c(1,1,2,2,3,3)
ar<-3
#n<-1
#m<-10

#load data
syn_resids_tmax<-array(0,c(m,length(ix2),64))
syn_resids_tmin<-array(0,c(m,length(ix2),64))
source('GL_maineqs.R')
obtmax_arr<-readRDS('data/basic/obs_mat_fit.rds')
obtmax_arr<-obtmax_arr+273
var_coefs_temp<-readRDS(paste(src,'/var_coefs_temp.rds',sep=''))
gl_par_arr2<-readRDS(paste(src,'/gl_par_arr.rds',sep=''))
mn_tmax<-readRDS('data/basic/mn_tmax.rds') #previously calculated means to rebias errors
mn_tmin<-readRDS('data/basic/mn_tmin.rds') #previously calculated means to rebias errors

#create (n x m) no. of synthetic error matrices
for(s in 1:n){
  #load Schaake shuffled and kNN sampled a_t
  knn_at_samp_arr_tmax<-readRDS(paste(src,'/knn_at_samp_arr_tmax_',s,'.rds',sep=''))
  knn_at_samp_arr_tmin<-readRDS(paste(src,'/knn_at_samp_arr_tmin_',s,'.rds',sep=''))
  for(d in 1:m){
    for(l in 1:2){
      mat<-cbind(knn_at_samp_arr_tmax[d,sn[[l]],],knn_at_samp_arr_tmin[d,sn[[l]],])
      coeff<-var_coefs_temp[l,,]
      gl_par_arr<-rbind(gl_par_arr2[l+2,,],gl_par_arr2[l+4,,])
      ob<-cbind(obtmax_arr[sn[[l]],],obtmax_arr[sn[[l]],])
      syn_resid_mat<-matrix(0,ncol=128,nrow=dim(mat)[1])
      #create synthetic errors across all dimensions
      for(j in 1:128){
        for(k in (ar+1):dim(mat)[1]){
          syn_resid_mat[k,j]<-t(matrix(c(syn_resid_mat[(k-1),],syn_resid_mat[(k-2),],syn_resid_mat[(k-3),]))) %*% matrix(coeff[j,]) + 
            sigma_t(gl_par_arr[j,1],gl_par_arr[j,2],(0.01*ob[k,j]))*mat[k,j]
          #if(j>64){
            #while((ob[k,j] - syn_resid_mat[k,j])>(ob[k,j-64] - syn_resid_mat[k,j-64])){
              #ats<-rsged(1,mean=0,sd=1,nu=(2/(1+gl_par_arr[j,3])),xi=gl_par_arr[j,4])
              #syn_resid_mat[k,j]<-sigma_t(gl_par_arr[j,1],gl_par_arr[j,2],ob[k,j]) * ats
            #}
          #}
        }
      }
      syn_resids_tmax[d,sn[[l]],]<-syn_resid_mat[,1:64] + matrix(rep(mn_tmax[sn2[[l]],],length(sn[[l]])),ncol=64,byrow=T) #add mean to 'rebias' errors
      syn_resids_tmin[d,sn[[l]],]<-syn_resid_mat[,65:128] + matrix(rep(mn_tmin[sn2[[l]],],length(sn[[l]])),ncol=64,byrow=T) 
    }
    print(paste(d,Sys.time()))
  }
  saveRDS(syn_resids_tmax,paste(out,'/syn_resids_tmax_',s,'.rds',sep=''))
  saveRDS(syn_resids_tmin,paste(out,'/syn_resids_tmin_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

print(paste('S1b complete',Sys.time()))

#--------------------------------------------------------------------------------------------------------------------
#2) Create synthetic forecasts with synthetic errors
#2a.PRECIP
#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
sn<-list(cold,warm)
#n<-1
#m<-10
syn_forc_prcp<-array(0,c(m,length(ix2),64))

#load data
obprcp_arr<-readRDS('data/basic/obs_map_fit.rds')
prcp_loess<-readRDS('data/basic/prcp_loess.rds')

#Note: PRECIP has an extra step so that occurrence/non-occurrence structure is captured
#appropriately in synthetic forecasts

#2a. part 1: Make raw synthetic forecasts (n x m) times over
for(s in 1:n){
  syn_resids_prcp<-readRDS(paste(out,'/syn_resids_prcp_',s,'.rds',sep=''))
  for(d in 1:m){
    for(k in 1:2){
      #subtract forecast errors from observations to make forecasts
      syn_forc_prcp[d,sn[[k]],]<-prcp_loess[sn[[k]],] - syn_resids_prcp[d,sn[[k]],]
    }
  }
  print(paste(s,Sys.time()))
  saveRDS(syn_forc_prcp,paste(out,'/syn_forc_prcp_',s,'.rds',sep=''))
}

#2a. part 2: Impose empirical occurrence/non-occurrence structure on result
for(s in 1:n){
  #array of occurrence/non-occurrence entries for each sample series
  knn_prcp_occ<-readRDS(paste(src,'/knn_prcp_occ_',s,'.rds',sep=''))
  syn_forc_prcp<-readRDS(paste(out,'/syn_forc_prcp_',s,'.rds',sep=''))
  syn_resids_prcp<-readRDS(paste(out,'/syn_resids_prcp_',s,'.rds',sep=''))
  for(d in 1:m){
    syn_f<-syn_forc_prcp[d,,]
    #ensure all empirical non-occurrence events are set to zero
    syn_f[which(knn_prcp_occ[d,,]==0)]<-0
    syn_forc_prcp[d,,]<-syn_f
    syn_resids_prcp[d,,]<-prcp_loess - syn_forc_prcp[d,,] #recalculate errors after occurrence corrections
  }
  #save the corrected result
  print(paste(s,Sys.time()))
  saveRDS(syn_resids_prcp,paste(out,'/syn_resids_prcp_',s,'.rds',sep=''))
  saveRDS(syn_forc_prcp,paste(out,'/syn_forc_prcp_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

#2b.TMAX
#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
sn<-list(cold,warm)
#n<-1
#m<-10

#load observations
obtmax_arr<-readRDS('data/basic/obs_mat_fit.rds')

#Create forecasts from errors
syn_forc_tmax<-array(0,c(m,length(ix2),64))

for(s in 1:n){
  #load synthetic errors
  syn_resids_tmax<-readRDS(paste(out,'/syn_resids_tmax_',s,'.rds',sep=''))
  for(d in 1:m){
    for(k in 1:2){
      #subtract forecast errors from observations to make forecasts
      syn_forc_tmax[d,sn[[k]],]<-obtmax_arr[sn[[k]],] - syn_resids_tmax[d,sn[[k]],]
    }
  }
  print(paste(s,Sys.time()))
  saveRDS(syn_forc_tmax,paste(out,'/syn_forc_tmax_',s,'.rds',sep=''))
}

#rm(list=ls());gc()

#2c.TMIN
#basic parameters
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)
cold<-sort(c(which(ixx2$mon>8),which(ixx2$mon<3)))
warm<-1:length(ixx2)
warm<-warm[-c(cold)]
sn<-list(cold,warm)
#n<-1
#m<-10

#load observations
obtmin_arr<-readRDS('data/basic/obs_mat_fit.rds')

#Create forecasts from errors
syn_forc_tmin<-array(0,c(m,length(ix2),64))

for(s in 1:n){
  #load synthetic errors
  syn_resids_tmin<-readRDS(paste(out,'/syn_resids_tmin_',s,'.rds',sep=''))
  for(d in 1:m){
    for(k in 1:2){
      #subtract forecast errors from observations to make forecasts
      syn_forc_tmin[d,sn[[k]],]<-obtmin_arr[sn[[k]],] - syn_resids_tmin[d,sn[[k]],]
    }
  }
  saveRDS(syn_forc_tmin,paste(out,'/syn_forc_tmin_',s,'.rds',sep=''))
}

rm(list=ls());gc()
print(paste('S2 complete',Sys.time()))
##########################################################END##############################################################