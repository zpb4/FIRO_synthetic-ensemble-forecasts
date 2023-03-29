setwd('h:/firo_lamc/plot/')
library(stringr)
samps<-100

results<-read.csv('results_v10_sumrsamp30_all.csv')
pfo<-read.csv('PfoResults_Final.csv')
pfo_store<-pfo$stor_mendo;pfo_store_comp<-pfo_store[which(pfo$date_time=='10/16/1985 12:00'):which(pfo$date_time=='9/30/2010 12:00')]
cnrfc<-read.csv('EfoResults_CNRFC_Final.csv')
cnrfc_store<-cnrfc$stor_mendo;cnrfc_store_comp<-cnrfc_store[which(cnrfc$date_time=='10/16/1985 12:00'):which(cnrfc$date_time=='9/30/2010 12:00')]

r1<-results$name_scenario[1]
r1_strip<-str_remove(r1,'_1')
len<-length(which(results$name_scenario==r1))
r1_date_st<-results$date_time[1]
date_start<-str_remove(r1_date_st,' 12:00')
r1_date_end<-results$date_time[len]
date_end<-str_remove(r1_date_end,' 12:00')

ix<-seq(as.Date('1985-10-16'),as.Date('2010-09-30'),'day') #same length
ixx<-as.POSIXlt(ix)

store_array<-array(NA,c(len,samps))
rownames(sim_array)<-as.character(ix)
spill_array<-array(NA,c(len,samps))
rownames(spill_array)<-as.character(ix)

for(n in 1:samps){
  store_array[,n]<-results$stor_mendo[which(results$name_scenario==paste(r1_strip,'_',n,sep=''))]
  spill_array[,n]<-results$spill[which(results$name_scenario==paste(r1_strip,'_',n,sep=''))]
}

day<-'2006-01-03'
idx<-which(ix==day)

length(which(spill_array[idx,]>0))
length(which(store_array[idx,]>116500))


#plot zoomed in spill events for SI
evt_idx<-which(ix=='1986-02-18')
evt<-ix[evt_idx]
png(paste('figs/spill_no-spill_comp_',evt,'.png',sep=''),width=768,height = 924)
par(mfrow=c(3,1),mar=c(4,5,3,1),cex.main=2,cex.lab=2,cex.axis=2)
start<-ix[evt_idx-15]
end<-ix[evt_idx+10]
idx2<-which(ix==start):which(ix==end)


syn_sset<-store_array[idx2,]

spill_idx<-apply(syn_sset,2,function(x){
  t<-any(x>116500)
  if(length(t[t==T])==0){y<-0}
  if(length(t[t==T])==1){y<-1}
  if(length(t[t==T])==2){y<-2}
  if(length(t[t==T])==3){y<-3}
  if(length(t[t==T])==4){y<-4}
  if(length(t[t==T])==5){y<-5}
  return(y)})

library(scales)
my_orng<-alpha('orange',alpha=0.5)

plot(1:length(idx2),pfo_store_comp[idx2]/1000,type='l',ylim=c(50,125),col='green',lwd=2,axes=F,xlab='Day',ylab='Storage (TAF)',main=paste(evt,'Spill'))
axis(1,at=seq(4,length(idx2),4),labels = start+seq(3,length(idx2),4))
axis(2,at=seq(50,125,25),labels=seq(50,125,25))
box(which='plot')
abline(h=116.5,col='black',lty=2,lwd=2)
abline(v=16,col='red',lty=3)
for(i in 1:100){
  if(spill_idx[i]==0){cl=my_orng}
  if(spill_idx[i]==1){cl='sienna4'}
  if(spill_idx[i]>1){cl='deeppink4'}
  lines(1:length(idx2),syn_sset[,i]/1000,col=cl)
}
lines(1:length(idx2),cnrfc_store_comp[idx2]/1000,col='blue',lwd=3)
lines(1:length(idx2),pfo_store_comp[idx2]/1000,col='green',lwd=3)
text(5,120,'Spillway: 116.5 TAF',cex=2)
legend('bottomright',c('HEFS','PFO','syn_GEFS spill','syn-GEFS no-spill'),col=c('blue','green','sienna4',my_orng),lwd=c(3,3,1,1),cex=2)

#non-spill only
plot(1:length(idx2),pfo_store_comp[idx2]/1000,type='l',ylim=c(50,125),col='green',lwd=2,axes=F,xlab='Day',ylab='Storage (TAF)',main='Non-spills only')
axis(1,at=seq(4,length(idx2),4),labels = start+seq(3,length(idx2),4))
axis(2,at=seq(50,125,25),labels=seq(50,125,25))
box(which='plot')
abline(h=116.5,col='black',lty=2,lwd=2)
for(i in 1:100){
  if(spill_idx[i]==0){lines(1:length(idx2),syn_sset[,i]/1000,col='orange')}
}
lines(1:length(idx2),cnrfc_store_comp[idx2]/1000,col='blue',lwd=3)
lines(1:length(idx2),pfo_store_comp[idx2]/1000,col='green',lwd=3)

#spill only
plot(1:length(idx2),pfo_store_comp[idx2]/1000,type='l',ylim=c(50,125),col='green',lwd=2,axes=F,xlab='Day',ylab='Storage (TAF)',main='Spills only')
axis(1,at=seq(4,length(idx2),4),labels = start+seq(3,length(idx2),4))
axis(2,at=seq(50,125,25),labels=seq(50,125,25))
box(which='plot')
abline(h=116.5,col='black',lty=2,lwd=2)
for(i in 1:100){
  if(spill_idx[i]>=1){lines(1:length(idx2),syn_sset[,i]/1000,col='sienna4')}
}
lines(1:length(idx2),cnrfc_store_comp[idx2]/1000,col='blue',lwd=3)
lines(1:length(idx2),pfo_store_comp[idx2]/1000,col='green',lwd=3)

dev.off()

########################################END#####################################################