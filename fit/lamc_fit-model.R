args <- commandArgs()
print(args)
s <- as.numeric(args[6])

e_idx<-cbind(seq(1,51,10),seq(10,61,10))
e_idx[6,2]<-61

print(paste(s,'start',Sys.time()))

library(fGarch)
library(BigVAR)
#library(doParallel)
#library(abind)

#parallel::detectCores()

#n.cores <- parallel::detectCores()-2
#my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
#print(my.cluster)

#doParallel::registerDoParallel(cl = my.cluster)
#foreach::getDoParRegistered()

#ens_num <- 61
leads <- 14 #number of leads desired

#1) Read in observed and simulated inflows
inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]

obs[obs<0]<-0
obs_mat<-cbind(matrix(rep(obs[,1],leads),ncol=leads),matrix(rep(obs[,2],leads),ncol=leads),matrix(rep(obs[,3],leads),ncol=leads))

loess_mat<-readRDS('data/fit/loess_mat_loessv3-1_hopper_12z.rds')
norm_fit<-readRDS('data/fit/norm_fit_loessv3-1_hopper_12z.rds')
sd_arr<-readRDS('data/fit/sd_arr_loessv3-1_hopper_12z.rds')

#seasonal inference stuff
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

lin_mod<-function(intcpt,slope,x){
  y = intcpt + slope * x
  return(y)
}

param_lst<-vector('list',4)

#foreach(e = e_idx[s,1]:e_idx[s,2],.packages = c('fGarch','BigVAR')) %dopar% {
for(e in e_idx[s,1]:e_idx[s,2]){
#2) Normalize
rresids<-readRDS(paste('data/fit/rresids_loessv3-1_hopper_ens_',e,'_12z.rds',sep=''))
nresids<-vector('list',12)
  
for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  obmat<-loess_mat[seas,]
  rresid_mat<-rresids[[i]]
  nresid_mat<-array(NA,dim(rresid_mat))
  for(k in 1:(leads*3)){
    ob<-obmat[,k]
    res<-rresid_mat[,k]
    norm_vec<-lin_mod(norm_fit[[i]][[k]][1],norm_fit[[i]][[k]][2],ob)
    norm_vec[norm_vec<=0]<-sd_arr[i,k]
    nresid_mat[,k]<-res / norm_vec
  }
  nresids[[i]]<-nresid_mat
}
  
#saveRDS(nresids,paste('data/fit/nresids_loessv3_hopper_ens_',e,'_12z.rds',sep=''))
param_lst[[1]]<-nresids

#3) VAR model to create uncorrelated residuals
ar <- 3
uc_resid <- vector('list',12)
var_coefs <- array(NA,c(12,leads*3,(leads*3*ar+1)))
  
#Calculate uncorrelated matrices
#print(Sys.time())
for(i in 1:12){
  rresid_mat<-nresids[[i]]
  #m1 = constructModel(rresid_mat, p = ar, struct = "BGR", gran = c(25, 10),IC = F,
                        #verbose = F, VARX = list(),separate_lambdas = F,model.controls=list(intercept = F, MN=F))
  m1 = constructModel(rresid_mat, p = ar, struct = "Basic", gran = c(25, 10),IC = F,
                        verbose = F, VARX = list(),separate_lambdas = T, rolling_oos=T,model.controls=list(intercept = F, MN=F))
  m1_res = cv.BigVAR(m1)
    
  var_coefs[i,,]<-m1_res@betaPred
  resids<-m1_res@fitted #calculate decorr residuals from 
    
  mat<-matrix(0,dim(rresid_mat)[1],dim(rresid_mat)[2])
  mat[(ar+1):dim(mat)[1],]<-resids #resids set to 0 for first columns = max lag
    
  uc_mat<-rresid_mat - mat
    
  uc_resid[[i]]<-uc_mat
}
  
#saveRDS(var_coefs,paste('data/fit/var_coefs_loessv3_hopper_ens_',e,'_12z.rds',sep=''))

param_lst[[2]]<-var_coefs

print(paste(e,'varfit complete',Sys.time()))
#4) GL SGED Model
#monthly model fits for each lead time
gl_par_arr<-array(NA,c(12,leads*3,4))
at_lst<-vector('list',12)
  
for(i in 1:12){
  uc_sim<-uc_resid[[i]]
  seas<-which(ix2$mon==(i-1))
  at_arr<-array(NA,c(dim(uc_sim)[1],leads*3))
  for(j in 1:(leads*3)){
      
    gl_mle<-sgedFit(uc_sim[,j])
      
    at<-uc_sim[,j]
    gl_par_arr[i,j,]<-gl_mle$par
    at_arr[,j]<-at
  }
  at_lst[[i]]<-at_arr
}
  
#saveRDS(gl_par_arr,paste('data/fit/gl_par_arr_loessv3_hopper_ens_',e,'_12z.rds',sep=''))
#saveRDS(at_lst,paste('data/fit/at_lst_loessv3_hopper_ens_',e,'_12z.rds',sep=''))

param_lst[[3]]<-gl_par_arr
param_lst[[4]]<-at_lst

saveRDS(param_lst,paste('data/fit/param_lst_loessv3-1_hopper_ens_',e,'_12z.rds',sep=''))
}

print(paste(s,'end',Sys.time()))

rm(list=ls());gc()
###########################################END################################