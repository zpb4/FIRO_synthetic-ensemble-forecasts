library(feather)
library(doParallel)
library(abind)

parallel::detectCores()

n.cores <- parallel::detectCores()-2
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

print(paste('start',Sys.time()))

inf<-read.csv('~/syn_forecast_test/data/LAMC_local.csv')
obs<-inf[which(inf$GMT=='10/1/1948 12:00'):which(inf$GMT=='9/30/1985 12:00'),4:6]
ix<-seq(as.Date('1948-10-01'),as.Date('1985-09-30'),'day')

ens_num<-61
nsets<-10
nsamps<-10

samp_vec<-matrix(1:(nsets*nsamps),ncol=nsamps,byrow=T)

outfile_hefs<-'corr_v10_sumrsamp30'

if(dir.exists(paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc',sep=''))==T){
  unlink(paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc',sep=''),recursive = T)}
dir.create(paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc',sep=''))

if(dir.exists(paste('~/syn_forecast_test/output/',outfile_hefs,'/npz_phc',sep=''))==T){
  unlink(paste('~/syn_forecast_test/output/',outfile_hefs,'/npz_phc',sep=''),recursive = T)}
dir.create(paste('~/syn_forecast_test/output/',outfile_hefs,'/npz_phc',sep=''))

foreach(m = 1:nsets,.packages='feather') %dopar% {
  syn_forc_in<-readRDS(paste('~/syn_forecast_test/output/',outfile_hefs,'/syn_hefs_flow_4885_',outfile_hefs,'_ens_',m,'_12z.rds',sep=''))

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

    for(i in 1:14){
      syn_hop_out[,(i+1),]<-syn_hop[i,c((i+1):dim(syn_hop)[2],rep(dim(syn_hop)[2],(i))),1:ens_num]
      syn_lam_out[,(i+1),]<-syn_lam[i,c((i+1):dim(syn_hop)[2],rep(dim(syn_hop)[2],(i))),1:ens_num]
      syn_uka_out[,(i+1),]<-syn_uka[i,c((i+1):dim(syn_hop)[2],rep(dim(syn_hop)[2],(i))),1:ens_num]
    }

    if(dir.exists(paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc/samp',samp_vec[m,n],sep=''))==T){
      unlink(paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc/samp',samp_vec[m,n],sep=''),recursive = T)}
    dir.create(paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc/samp',samp_vec[m,n],sep=''))
    
    for(e in 1:ens_num){
      syn_hop_fthr<-as.data.frame(syn_hop_out[,,e])
      syn_lam_fthr<-as.data.frame(syn_lam_out[,,e])
      syn_uka_fthr<-as.data.frame(syn_uka_out[,,e])
      colnames(syn_hop_fthr)<-paste(0:14,'d')
      colnames(syn_lam_fthr)<-paste(0:14,'d')
      colnames(syn_uka_fthr)<-paste(0:14,'d')
      rownames(syn_hop_fthr)<-as.Date(ix)
      rownames(syn_lam_fthr)<-as.Date(ix)
      rownames(syn_uka_fthr)<-as.Date(ix)
      write_feather(syn_hop_fthr, paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc/samp',samp_vec[m,n],'/syn-hop_',e,'.feather',sep=''))
      write_feather(syn_lam_fthr, paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc/samp',samp_vec[m,n],'/syn-lam_',e,'.feather',sep=''))
      write_feather(syn_uka_fthr, paste('~/syn_forecast_test/output/',outfile_hefs,'/feather_phc/samp',samp_vec[m,n],'/syn-uka_',e,'.feather',sep=''))
    }
  }
}

print(paste('end',Sys.time()))

#############################################END############################################

