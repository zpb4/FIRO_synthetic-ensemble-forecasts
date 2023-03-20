setwd('h:/firo_lamc/temp_precip/')
ix2<-seq(as.Date('1985-01-01'),as.Date('2017-09-30'),'day') #same length
ixx2<-as.POSIXlt(ix2)

src<-'out/mag-samp_loess'
out<-'mag-samp_loess'

m<-1
s<-10

form<-read.table('fews_format.txt',sep=',',stringsAsFactors = F)

dt<-vector('list',length(ix2))

hrs<-rep(c('06:00:00','12:00:00','18:00:00','00:00:00'),16)

for(i in 1:(length(ix2)-16)){
  days<-c(rep(as.character(ix2[i]),3),rep(as.character(ix2[i+1]),4),rep(as.character(ix2[i+2]),4),rep(as.character(ix2[i+3]),4),
          rep(as.character(ix2[i+4]),4),rep(as.character(ix2[i+5]),4),rep(as.character(ix2[i+6]),4),rep(as.character(ix2[i+7]),4),
          rep(as.character(ix2[i+8]),4),rep(as.character(ix2[i+9]),4),rep(as.character(ix2[i+10]),4),rep(as.character(ix2[i+11]),4),
          rep(as.character(ix2[i+12]),4),rep(as.character(ix2[i+13]),4),rep(as.character(ix2[i+14]),4),rep(as.character(ix2[i+15]),4),as.character(ix2[i+16]))
  dt[[i]]<-paste(days,hrs)
}

saveRDS(dt,'hefs_gefs_date-time.rds')

syn_forc_prcp<-readRDS(paste(src,'/syn_forc_prcp_',m,'.rds',sep=''))
syn_forc_tmax<-readRDS(paste(src,'/syn_forc_tmax_',m,'.rds',sep=''))
syn_forc_tmin<-readRDS(paste(src,'/syn_forc_tmin_',m,'.rds',sep=''))

idx<-seq(4,(15*4),4)

n<-(length(ix2)-16)

gefs_arr<-array(NA,c(3,n,64))

for(j in 1:s){

for(i in 1:n){
  synf_prcp<-syn_forc_prcp[j,i,]
  synf_tmax<-syn_forc_tmax[j,i,]
  synf_tmin<-syn_forc_tmin[j,i,]
  for(k in 1:15){
    synf_prcp[idx[k]:(idx[k]+3)]<-syn_forc_prcp[j,i+k,idx[k]:(idx[k]+3)]
    synf_tmax[idx[k]:(idx[k]+3)]<-syn_forc_tmax[j,i+k,idx[k]:(idx[k]+3)]
    synf_tmin[idx[k]:(idx[k]+3)]<-syn_forc_tmin[j,i+k,idx[k]:(idx[k]+3)]
  }
  synf_prcp[64]<-syn_forc_prcp[j,i+16,64]
  synf_tmax[64]<-syn_forc_tmax[j,i+16,64]
  synf_tmin[64]<-syn_forc_tmin[j,i+16,64]
  
  gefs_arr[1,i,]<-synf_prcp
  gefs_arr[2,i,]<-synf_tmax
  gefs_arr[3,i,]<-synf_tmin
}

ix3<-seq(as.Date('1985-01-01'),as.Date('2017-09-14'),'day') #same length
ixx3<-as.POSIXlt(ix3)

library(stringr)

n<-length(ix3)

if(dir.exists(paste('csv/gefs_csv_',out,'_ens',j,sep=''))==T){
  unlink(paste('csv/gefs_csv_',out,'_ens',j,sep=''),recursive = T)}
dir.create(paste('csv/gefs_csv_',out,'_ens',j,sep=''))

for(i in 1:n){
  fname<-paste(ixx3[i]$year+1900,str_pad(ixx3[i]$mon+1,2,'left','0'),str_pad(ixx3[i]$mday,2,'left','0'),'12_GEFSv10.csv',sep='')
  form_out<-form[,-c(1)]
  form_out[4:67,1]<-dt[[i]]
  form_out[4:67,1]<-gefs_arr[1,i,]
  form_out[4:67,2]<-gefs_arr[2,i,]
  form_out[4:67,3]<-gefs_arr[3,i,]
  form_out<-data.frame(form_out)
  write.table(form_out,paste('csv/gefs_csv_',out,'_ens',j,'/',fname,sep=''),row.names = c('Location Names','Location Ids','Time',dt[[i]]),quote=F,col.names = F,sep=',')
}

print(paste(j,Sys.time()))
}