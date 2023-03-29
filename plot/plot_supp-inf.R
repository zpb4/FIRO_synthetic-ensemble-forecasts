e_cons<-c(9.5,10.5,9.5,10,9.5,11,10,10.5,9.5)
sum(e_cons)

e_ovr<-c(18,16.5,14,12,9.5,7,6,4,3)
sum(e_ovr)

e_ovrd<-c(3,7,8,14,17,16,11,9,5)
sum(e_ovrd)

e_underd<-c(20,12,8,3,2,4,7,13,21)
sum(e_underd)

e_udr<-rev(e_ovr)

roll_sum<-function(x){out<-c();for(i in 1:length(x)){out[i]<-sum(x[1:i])};return(out)}

par(mfrow=c(3,3),mar=c(0,1,4,1),cex.main=1.5)

plot(NULL)

barplot(e_ovr,axes=F,ylim=c(0,25),main='Overforecast bias',space=0)
segments(0,10,9,10,lty=2,lwd=3)
rect(0,0,9,25)

plot(NULL)

barplot(e_ovrd,axes=F,ylim=c(0,25),main='Overdispersion',space=0)
segments(0,10,9,10,lty=2,lwd=3)
rect(0,0,9,25)

barplot(e_cons,axes=F,ylim=c(0,25),main='Rank uniformity',space=0)
segments(0,10,9,10,lty=2,lwd=3)
rect(0,0,9,25)

barplot(e_underd,axes=F,ylim=c(0,25),main='Underdispersion',space=0)
segments(0,10,9,10,lty=2,lwd=3)
rect(0,0,9,25)

plot(NULL)

barplot(e_udr,axes=F,ylim=c(0,25),main='Underforecast bias',space=0)
segments(0,10,9,10,lty=2,lwd=3)
rect(0,0,9,25)

plot(NULL)


#cumulative histogram

par(mfrow=c(3,3),mar=c(0,1,4,1),cex.main=1.5)

plot(NULL)

plot(0:9,c(0,roll_sum(e_ovr)/90),axes=F,ylim=c(0,1),type='l',main='Overforecast bias',col='grey',lwd=3)
abline(0,0.1111,lwd=3,lty=2)
box(which='plot')

plot(NULL)

plot(0:9,c(0,roll_sum(e_ovrd)/90),axes=F,ylim=c(0,1),type='l',main='Overdispersion',col='grey',lwd=3)
abline(0,0.1111,lwd=3,lty=2)
box(which='plot')

plot(0:9,c(0,roll_sum(e_cons)/90),axes=F,ylim=c(0,1),type='l',main='Rank uniformity',col='grey',lwd=3)
abline(0,0.1111,lwd=3,lty=2)
box(which='plot')

plot(0:9,c(0,roll_sum(e_underd)/90),axes=F,ylim=c(0,1),type='l',main='Underdispersion',col='grey',lwd=3)
abline(0,0.1111,lwd=3,lty=2)
box(which='plot')

plot(NULL)

plot(0:9,c(0,roll_sum(e_udr)/90),axes=F,ylim=c(0,1),type='l',main='Underforecast bias',col='grey',lwd=3)
abline(0,0.1111,lwd=3,lty=2)
box(which='plot')

plot(NULL)