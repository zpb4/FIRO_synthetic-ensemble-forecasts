
#setwd('z:/syn_forecast_test/output/')
args <- commandArgs()
print(args)
z=as.numeric(args[6])

#z<-1

print(paste('start',Sys.time()))

n<-10
leads<-14
typ<-'v10-2_sumrsamp30'
env_scale<-3
init_chk<-F

syn_flow<-readRDS(paste(typ,'/syn_hefs_flow_8510_',typ,'_ens_',z,'_12z.rds',sep=''))
syn_flow_new<-syn_flow

inf<-read.csv('~/syn_forecast_test/data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/15/1985 12:00'):which(inf$GMT=='9/30/2010 12:00'),4:6]
obs[obs<0]<-0

ix<-seq(as.Date('1985-10-15'),as.Date('2010-09-30'),'day')
ix2<-as.POSIXlt(ix)
h<-c(2:n,1)
ens_num<-61

if(init_chk==T){
for(m in 1:n){
  syn_hefs_flow1<-syn_flow[,,,m]
  #error checking
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    ob_seas<-obs$LAMC1[seas]
    rvec<-syn_hefs_flow1[15:28,seas,]
    print(max(ob_seas))
    print(max(rvec))
    for(j in 1:leads){
      for(k in 1:ens_num){
        if(any(rvec[j,,k] > (env_scale*max(ob_seas)))==T) {print(paste('Bad member - month',i,'lead',j,'ens',k))} 
      }
    }
  }
  print(m)
}
}

for(m in 1:n){
  syn_hefs_flow1<-syn_flow[,,,m]
  bup_samp1<-1:10;bup_samp1<-bup_samp1[-c(m)];bsamp1<-sample(bup_samp1,1)
  syn_hefs_flow2<-syn_flow[,,,bsamp1]
  syn_hefs_new<-syn_hefs_flow1
  
  
  #error checking
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    rvec<-syn_hefs_flow1[15:28,seas,]
    samp_vec<-syn_hefs_flow2[15:28,seas,]
    syn_new<-syn_hefs_new[15:28,seas,]
    ob_seas<-obs$LAMC1[seas]
    for(j in 1:leads){
      bvec1<-c()
      bvec2<-c()
      s<-0
      b<-0
      d<-0
      while(s==0 & b<100){
        for(k in 1:ens_num){
          #med_vec<-apply(rvec[j,,],1,median)
          if(any(rvec[j,,k] > (env_scale*max(ob_seas)))==T) {bvec1<-c(bvec1,k)}
          if(any(samp_vec[j,,k] > (env_scale*max(ob_seas)))==T) {bvec2<-c(bvec2,k)}
        }
        if(length(bvec2)>(ens_num-5)){
          s<-0
          b<-b+1
          bup_samp1<-1:10;bup_samp1<-bup_samp1[-c(m)];bsamp1<-sample(bup_samp1,1)
          syn_hefs_flow2<-syn_flow[,,,bsamp1]
          samp_vec<-syn_hefs_flow2[15:28,seas,]}
        else{s<-1}
      }
      if(s==0 & b==100){
        bvec1<-c()
        bvec2<-c()
        t<-0
        l<-0
        while(t==0 & l<100){
          for(k in 1:ens_num){
            if(any(rvec[j,,k] > (env_scale*max(ob_seas)))==T) {bvec1<-c(bvec1,k)}
            if(j<leads){
              if(any(samp_vec[j+1,,k] > (env_scale*max(ob_seas)))==T) {bvec2<-c(bvec2,k)}
              d<-1
            }
            if(j==leads){
              if(any(samp_vec[j-1,,k] > (env_scale*max(ob_seas)))==T) {bvec2<-c(bvec2,k)}
              d<-(-1)
            }
          }
          if(length(bvec2)>(ens_num-5)){
            t<-0
            l<-l+1
            bup_samp1<-1:10;bup_samp1<-bup_samp1[-c(m)];bsamp1<-sample(bup_samp1,1)
            syn_hefs_flow2<-syn_flow[,,,bsamp1]
            samp_vec<-syn_hefs_flow2[15:28,seas,]}
          else{t<-1}
        }
      }
      if(length(bvec1)>0){
        gvec<-1:61
        if(length(bvec2)>0){gvec<-gvec[-c(bvec2)]}
        input<-sample(gvec,length(bvec1),replace=T)
        syn_new[j,,bvec1]<-samp_vec[j+d,,input]}
    }
    
    syn_hefs_new[15:28,seas,]<-syn_new
    
    rvec<-syn_hefs_new[15:28,seas,]
    for(j in 1:leads){
      for(k in 1:ens_num){
        if(any(rvec[j,,k] > (env_scale*max(ob_seas)))==T) {print(paste('Corr-ens_bad member - month',i,'lead',j,'ens',k))} 
      }
    }
  }
  syn_flow_new[15:28,,,m]<-syn_hefs_new[15:28,,]
  print(m)
}

saveRDS(syn_flow_new,paste('corr_',typ,'/syn_hefs_flow_8510_corr_',typ,'_ens_',z,'_12z.rds',sep=''))

print(paste('end',Sys.time()))


#################END######################



######################################END#############################################