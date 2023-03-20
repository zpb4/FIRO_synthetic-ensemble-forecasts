#7) Synthetic forecast generation
args <- commandArgs()
print(args)
z=as.numeric(args[6])
print(paste('start',Sys.time()))

library(fGarch)
library(doParallel)
library(abind)

parallel::detectCores()

n.cores <- parallel::detectCores()-2
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

ar <- 3
ens_num <- 61
leads <- 14
env_scale<-10
#n<-100

midx<-cbind(seq(1,100,10),seq(10,100,10))

msst<-midx[z,]

obs_nep_hist<-readRDS('data/obs_nep_hist_12z.rds')
obs_nep<-readRDS('data/obs_nep_12z.rds')

inf<-read.csv('data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]

obs[obs<0]<-0

ob_samp<-array(NA,dim(obs))

for(i in 1:length(ob_samp[,1])){
  idx<-i:min(length(ob_samp[,1]),(i+14))
  ob_samp[i,]<-apply(obs[idx,],2,sum)
}

obs_mat<-cbind(matrix(rep(obs[,1],leads),ncol=leads),matrix(rep(obs[,2],leads),ncol=leads),matrix(rep(obs[,3],leads),ncol=leads))

#Generate new out of sample data with KNN method
new_obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]
new_obs[new_obs<0]<-0

new_ob_samp<-array(NA,dim(new_obs))

for(i in 1:length(new_ob_samp[,1])){
  idx<-i:min(length(new_ob_samp[,1]),(i+14))
  new_ob_samp[i,]<-apply(new_obs[idx,],2,sum)
}

new_obs_mat<-cbind(matrix(rep(new_obs[,1],leads),ncol=leads),matrix(rep(new_obs[,2],leads),ncol=leads),matrix(rep(new_obs[,3],leads),ncol=leads))

#seasonal inference stuff
ixx<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix3<-as.POSIXlt(ixx)

ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)

yr_idx<-ix2$year
yr_idx_lst<-vector('list',12)

for(i in 1:12){
  seas3<-which(ix3$mon==(i-1))
  yr_idx_lst[[i]]<-yr_idx[seas3]
}

yr_seq<-c(rep(85,3),rep(86:109,each=12),rep(110,9))
mo_seq<-c(10:12,rep(1:12,length(86:109)),1:9)

knn<-30
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn)

norm_fit<-readRDS('data/fit/norm_fit_loessv3-1_hopper_12z.rds')
sd_arr<-readRDS('data/fit/sd_arr_loessv3-1_hopper_12z.rds')
cmean_mat<-readRDS('data/fit/loess_mat_loessv3-1_hopper_12z.rds')
#max_val<-readRDS('data/fit/max_val_loessv3-1_hopper_12z.rds')

lin_mod<-function(intcpt,slope,x){
  y = intcpt + slope * x
  return(y)
}

abnd<-function(x,y){abind(x,y,along=4)}

#1) Generate kNN indices
#syn_hefs_resid<-array(NA,c(leads*3,length(ix3),ens_num))
synflow_out <- foreach(m = msst[1]:msst[2],.combine='abnd',.packages='fGarch') %dopar% {

syn_hefs_flow<-array(NA,c(leads*3,length(ix3),ens_num))
  
knn_lst<-vector('list',12)
  
for(i in 1:12){
  knn_vec<-c()
  seas<-which(ix2$mon==(i-1))
  n_obs<-new_ob_samp[seas,]
  ob<-ob_samp[seas,]
  #ob0 <- apply(ob,1,sum)
  #id0 <- which(ob0 == 0)
  for(j in 1:length(n_obs[,1])){
    #if(sum(n_obs[j,]) == 0 & length(id0)>=1) {
      #s<-sample(id0,1); knn_vec[j]<-s} 
    #else {
      ob_smp<-n_obs[j,]
      y<-sqrt(apply((ob_smp - ob)^2,1,sum)) #find NEP closest by Euclidean distance
      x<-sort(y)
      x<-x[1:knn] #top 28 values are KNN
      s<-sample(x,1)#,prob=wts)
      id<-which(y==s)
      if(length(id)>1) {
        id <- sample(id,1)} #resample the sample for any duplicated values
      knn_vec[j]<-id#}
  }
  knn_lst[[i]]<-knn_vec
}
  
  
for(e in 1:ens_num){
  param_lst<-readRDS(paste('data/fit/param_lst_loessv3-1_hopper_ens_',e,'_12z.rds',sep=''))
  var_coefs<-param_lst[[2]]
  gl_par_arr<-param_lst[[3]]
  at_lst<-param_lst[[4]]
  

  #1) Create array of n Schaake Shuffled sequences for KNN
  syn_cop<-vector('list',12)

  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    mat<-at_lst[[i]]
    ecop<-apply(mat,2,function(x){rank(x,ties.method = 'random')})
    syn_ecop_mat<-array(NA,c(dim(mat)[1],(leads*3)))
    for(j in 1:(leads*3)){
      syn_at<-rsged(dim(mat)[1],mean=gl_par_arr[i,j,1],sd=gl_par_arr[i,j,2],nu=gl_par_arr[i,j,3],xi=gl_par_arr[i,j,4])
      r_syn_at<-rank(syn_at,ties.method = 'random')
      for(k in 1:length(r_syn_at)){
      syn_ecop_mat[k,j]<-syn_at[which(r_syn_at==ecop[k,j])]
      }
    }
    syn_cop[[i]]<-syn_ecop_mat
  }

#KNN process to generate array of 'ens_num' samples
  syn_cop_knn<-vector('list',12)
 
  for(i in 1:12){
    synecop_knn<-syn_cop[[i]]
    syncop<-array(NA,c(length(n_obs[,1]),leads*3))
    id<-knn_lst[[i]]
    syncop<-synecop_knn[id,]
    syn_cop_knn[[i]]<-syncop
  }
    
  app_mat<-syn_cop_knn[[mo_seq[1]]][1:3,]
    
  for(i in 1:length(yr_seq)){
    mat2<-syn_cop_knn[[mo_seq[i]]][which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
    coeff<-var_coefs[mo_seq[i],,]
    seas3<-which(ix3$mon==(mo_seq[i]-1) & ix3$year==yr_seq[i])
    seas<-which(ix3$mon==(mo_seq[i]-1))
    ob<-cmean_mat[seas3,]
    aob_seas<-obs_mat[seas,]
    max_vec<-apply(aob_seas,2,max)#/apply(aob_seas,2,mean)
    act_obs<-obs_mat[seas3,]
    syn_resid_mat<-matrix(0,ncol=leads*3,nrow=(dim(mat2)[1]+ar))
    syn_var_mat<-matrix(0,ncol=leads*3,nrow=(dim(mat2)[1]+ar))
    syn_var_mat[1:3,]<-app_mat
    for(j in 1:(leads*3)){
      for(k in (ar+1):(dim(mat2)[1]+ar)){
        bvar_res<-t(matrix(c(1,syn_var_mat[(k-1),],syn_var_mat[(k-2),],syn_var_mat[(k-3),]))) %*% matrix(coeff[j,]) + mat2[(k-ar),j]
        aobs<-act_obs[(k-ar),j]
        o<-ob[(k-ar),j]
        syn_var_mat[k,j]<-bvar_res[1,1]
        norm_val<-lin_mod(norm_fit[[mo_seq[i]]][[j]][1],norm_fit[[mo_seq[i]]][[j]][2],o)
        #ensure no normalization values less than zero (bad prediction from linear model)
        if(norm_val<=0){norm_val<-sd_arr[mo_seq[i],j]}
        if(o>max_vec[j]){o<-aobs}
        
        res<-bvar_res[1,1] * norm_val 
        
        #if residual will exceed reasonable limits (1.5X the max observed forecast)
        #try first to use a the more stable normalization constant
        if((o - res)>(env_scale*max_vec[j])){
          res<-bvar_res[1,1]*sd_arr[mo_seq[i],j]
        }
        
        #then try removing the VAR part to prevent runaway autocorrelation
        if((o - res)>(env_scale*max_vec[j])){
          res<-mat2[(k-ar),j]*sd_arr[mo_seq[i],j]
          syn_var_mat[k,j]<-res}
        
        #finally, try resampling from the distribution and pick from samples that don't exceed reasonable limits
        if((o - res)>(env_scale*max_vec[j])){
          new_at<-rsged(100,mean=gl_par_arr[mo_seq[i],j,1],sd=gl_par_arr[mo_seq[i],j,2],nu= gl_par_arr[mo_seq[i],j,3],xi=gl_par_arr[mo_seq[i],j,4])
          new_res<-new_at*sd_arr[mo_seq[i],j]
          samp_idx<-which((o-new_res)<(env_scale*max_vec[j]))
          res<-sample(new_res[samp_idx],1)
          syn_var_mat[k,j]<-res}
          
        #set to 0 if no work
        if((o - res)>(env_scale*max_vec[j])){
          res<-0
          syn_var_mat[k,j]<-res}
        
        
        syn_resid_mat[k,j]<-res
        
        #want to minimize 0 value forecasts
        #if residual from above will produce a negative value
        if(res > o){
          #generate a bunch of new samples from fitted distribution
          syn_at<-rsged(1000,mean=gl_par_arr[mo_seq[i],j,1],sd=gl_par_arr[mo_seq[i],j,2],nu= gl_par_arr[mo_seq[i],j,3],xi=gl_par_arr[mo_seq[i],j,4])
          #apply VAR model to all residuals
          var_res1<-t(matrix(c(1,syn_var_mat[(k-1),],syn_var_mat[(k-2),],syn_var_mat[(k-3),]))) %*% matrix(coeff[j,])
          #secondary set of residuals if can't get a suitable sample
          var_res2<-syn_at
          var_res<-(var_res1[1,1] + var_res2) * norm_val
          #determine positive samples which will produce a forecast flow value between the observation(cmean) and 0 flow
          var_idx<-which(var_res>=0 & var_res<=o)
          at_idx<-which(var_res2>=0 & var_res2<=o)
          #select first any residuals generated by full VAR model that will work
          if(length(at_idx)>0){syn_resid_mat[k,j]<-var_res2[sample(at_idx,1)]}
          #if no good samples there, try the raw a_t samples
          if(length(var_idx)>0){syn_resid_mat[k,j]<-var_res[sample(var_idx,1)]}
          #accept 0 as a last resort
          else{syn_resid_mat[k,j]<-o}
          #else{syn_resid_mat[k,j]<-(-min(abs(var_res2)))}
        }
        if((o - syn_resid_mat[k,j])>(env_scale*max_vec[j])){print(paste(o - syn_resid_mat[k,j],env_scale*max_vec[j],e,i,j,k,sep=','))}
      }
      #syn_hefs_resid[j,seas3,e]<-syn_resid_mat[(ar+1):k,j]
      syn_hefs_flow[j,seas3,e]<-cmean_mat[seas3,j] - syn_resid_mat[(ar+1):k,j]
    }
    app_mat<-syn_var_mat[(k-2):k,]
  }
    #print(paste(e,Sys.time()))
}

print(paste(m,Sys.time()))
return(syn_hefs_flow)
}

print(paste('end',Sys.time()))

saveRDS(synflow_out, paste('output/v10_sumrsamp30/syn_hefs_flow_8510_v10_sumrsamp30_ens_',z,'_12z.rds',sep=''))
#saveRDS(syn_hefs_resid, paste('output/syn_hefs_resid_8510_loessv3_hopper_ens_',m,'_12z.rds',sep=''))

rm(list=ls());gc()

###################################END######################################

