#args <- commandArgs()
#print(args)
#m=as.numeric(args[6])

library(abind)
library(fGarch)
library(BigVAR)

print(paste('start',Sys.time()))

ens_num <- 61
leads<-14 #number of leads desired
#1) Read in observed and simulated inflows

inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]
obs_tot<-inf[which(inf$GMT=='10/1/1948 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]

obs[obs<0]<-0
obs_tot[obs_tot<0]<-0
obs_mat<-cbind(matrix(rep(obs[,1],leads),ncol=leads),matrix(rep(obs[,2],leads),ncol=leads),matrix(rep(obs[,3],leads),ncol=leads))
loess_mat<-array(NA,dim(obs_mat))
  
#seasonal inference stuff
ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)
ix3<-seq(as.Date('1985-10-01'),as.Date('2010-09-30'),'day')
hefs_idx<-which(ix3=='1985-10-15'):which(ix3=='2010-09-30')

lamc_hefs<-readRDS('data/lamc_hefs_ens_forc_act-meteo_12z.rds')
hopc_hefs<-readRDS('data/hopc_hefs_ens_forc_act-meteo_12z.rds')
ukac_hefs<-readRDS('data/ukac_hefs_ens_forc_act-meteo_12z.rds')
lamc_hefs<-lamc_hefs[,hefs_idx,]*1000;hopc_hefs<-hopc_hefs[,hefs_idx,]*1000;ukac_hefs<-ukac_hefs[,hefs_idx,]*1000

o_hopc<-array(rep(obs$HOPC1),c(length(hefs_idx),ens_num))
o_lamc<-array(rep(obs$LAMC1),c(length(hefs_idx),ens_num))
o_ukac<-array(rep(obs$UKAC1),c(length(hefs_idx),ens_num))

loess_fit<-vector('list',leads)
loess_fit_sub<-vector('list',12)
loess_fit_ssub<-vector('list',3)

for(i in 1:leads){
  loess_fit[[i]]<-loess_fit_sub
  for(j in 1:12){
  loess_fit[[i]][[j]]<-loess_fit_ssub
  seas<-which(ix2$mon==(j-1))
  loess_fit[[i]][[j]][[1]]<-loess(apply(hopc_hefs[,seas,i],2,mean)~o_hopc[seas,1],span=1,degree = 2,family='symmetric',control=loess.control(surface='interpolate'))
  loess_fit[[i]][[j]][[2]]<-loess(apply(lamc_hefs[,seas,i],2,median)~o_lamc[seas,1],span=1,degree = 2,family='symmetric',control=loess.control(surface='interpolate'))
  loess_fit[[i]][[j]][[3]]<-loess(apply(ukac_hefs[,seas,i],2,mean)~o_ukac[seas,1],span=1,degree = 2,family='symmetric',control=loess.control(surface='interpolate'))
  }
}

saveRDS(loess_fit,'data/fit/loess_fit_loessv3-1_hopper_12z.rds')

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  for(k in 1:3){
    obs_inf<-matrix(rep(obs[seas,k],leads),ncol=leads,byrow=F)
    loess_inf<-array(NA,dim(obs_inf))
    for(j in 1:leads){
      ob_pred<-obs_inf[,j]
      cmn<-predict(loess_fit[[j]][[i]][[k]],ob_pred)
      zero_idx<-which(cmn<0)
      cmn[zero_idx]<-ob_pred[zero_idx]
      loess_inf[,j]<-cmn
    }
    loess_mat[seas,((k-1)*leads+1):(k*leads)]<-loess_inf
  }
}

saveRDS(loess_mat,'data/fit/loess_mat_loessv3-1_hopper_12z.rds')

max_val<-array(NA,c(12,leads*3))

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  max_val[i,1:14]<-apply(apply(hopc_hefs,c(2,3),max)[seas,],2,function(x){max(x,na.rm=T)})
  max_val[i,15:28]<-apply(apply(lamc_hefs,c(2,3),max)[seas,],2,function(x){max(x,na.rm=T)})
  max_val[i,29:42]<-apply(apply(ukac_hefs,c(2,3),max)[seas,],2,function(x){max(x,na.rm=T)})
}

saveRDS(max_val,'data/fit/max_val_loessv3-1_hopper_12z.rds')

for(e in 1:ens_num){
#1a. Raw Resid Matrix
rresids<-vector('list',12)
locs <- list(hopc_hefs,lamc_hefs,ukac_hefs)

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  rresid_mat<-matrix(ncol=leads * 3,nrow=length(seas))
  for(k in 1:3){
    sim_inf<-locs[[k]][e,,]
    loess_inf<-loess_mat[seas,((k-1)*leads+1):(k*leads)]
    rresid_mat[,((k-1)*leads+1):(k*leads)] <- loess_inf - sim_inf[seas,]
  }
  rresids[[i]]<-rresid_mat
}
saveRDS(rresids,paste('data/fit/rresids_loessv3-1_hopper_ens_',e,'_12z.rds',sep=''))
}

conc_resids<-vector('list',12)

for(i in 1:12){
  res_mat<-readRDS('data/fit/rresids_loessv3-1_hopper_ens_1_12z.rds')[[i]]
  for(e in 2:ens_num){
    rresids<-readRDS(paste('data/fit/rresids_loessv3-1_hopper_ens_',e,'_12z.rds',sep=''))
    res_add<-rresids[[i]]
    res_mat<-abind(res_mat,res_add,along=3)
  }
  conc_resids[[i]]<-res_mat
}

###############

sd_arr<-array(NA,c(12,leads*3))

lin_mod<-function(intcpt,slope,x){
  y = intcpt + slope * x
  return(y)
}

lin_fit<-function(pars,x,y){
  pred = pars[1] + pars[2] * x
  err = sum((y - pred)^2)
  return(err)
}

lbs<-c(0,0) #lower bounds for intcpt, slope
ubs<-c(1000,10) #upper bounds for intcpt, slope
sts<-c(0,1) #starting parameters for intcpt, slope

norm_fit<-vector('list',12)
norm_fit_sub<-vector('list',leads*3)

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  obmat<-loess_mat[seas,]
  #obmat<-ob_mat[seas,]
  rresid_mat<-conc_resids[[i]]
  rresid_mat<-sqrt(rresid_mat^2)
  rresid_md_mat<-apply(rresid_mat,c(1,2),mean)
  norm_fit[[i]]<-norm_fit_sub
  for(k in 1:(leads*3)){
    ob<-obmat[,k]
    res<-rresid_md_mat[,k]
    res_sort<-sort(res)
    lbs[1]<-mean(res_sort[1:round(0.1*length(seas))])
    ubs[1]<-mean(res)
    sd_arr[i,k]<-mean(rresid_md_mat[,k])
    
    lfit<-optim(par=sts,lin_fit,x=ob,y=res,
                method = 'L-BFGS-B',lower = lbs,upper = ubs,
                control = list(fnscale=1,maxit=100000))
    
    norm_fit[[i]][[k]]<-lfit$par
  }
}

saveRDS(norm_fit,'data/fit/norm_fit_loessv3-1_hopper_12z.rds')
saveRDS(sd_arr,'data/fit/sd_arr_loessv3-1_hopper_12z.rds')

print(paste('end',Sys.time()))
###########################################END################################