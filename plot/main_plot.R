setwd('h:/firo_lamc/')

#----------------------------------------------------------------
source('plot/plot_ensembles_hopper.r')
source('plot/plot_cumul_ensembles_hopper.r')

nsets<-10
nsamps<-10
lds<-c(5)
site<-'lamc'
forc<-'hefs'
show_leads<-T
show_legend<-T
leg_pos<-'bottomright'
show_xlab<-T
show_ylab<-T
show_yaxslab<-T
plt_par<-T
plot_samp<-F
plot_median<-F
plot_mean<-T
bounds<-c(0.05,0.95)
sgefs_seed<-1
shefs_set_seed<-1
shefs_seed<-5
hefs_samp_seed<-1
sgefs_samp_seed<-1
shefs_samp_seed<-1
event<-'1986-02-18'
outfile<-'corr_v10_sumrsamp30'

png(filename=paste('plot/figs/fig4_',outfile,'.png',sep=''),width = 1024, height = 768)

par(mfrow=c(3,2),mar=c(2,4,0,1),mgp=c(2.5,0.75,0),tcl=-0.25,cex.axis=1.75,cex.lab=2,font.lab=2,oma=c(2.5,5,4,0.5))

plt_ensembles(nsamps,nsets,lds,site,forc,show_leads,show_xlab=F,show_ylab=T,show_yaxslab=T,show_legend=T,plt_par,plot_samp,
                        plot_median,plot_mean=T,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                        sgefs_samp_seed,shefs_samp_seed,event,outfile)

plt_cumul_ensembles(nsamps,nsets,lds,site,forc,show_leads,show_xlab=F,show_ylab=F,show_yaxslab=T,show_legend=F,plt_par,plot_samp,
              plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
              sgefs_samp_seed,shefs_samp_seed,event,outfile)

plt_ensembles(nsamps,nsets,lds,site,forc='shefs',show_leads=F,show_xlab=F,show_ylab=T,show_yaxslab=T,show_legend=T,plt_par,plot_samp,
              plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
              sgefs_samp_seed,shefs_samp_seed,event,outfile)

plt_cumul_ensembles(nsamps,nsets,lds,site,forc='shefs',show_leads=F,show_xlab=F,show_ylab=F,show_yaxslab=T,show_legend=F,plt_par,plot_samp,
                    plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                    sgefs_samp_seed,shefs_samp_seed,event,outfile)

plt_ensembles(nsamps,nsets,lds,site,forc='sgefs',show_leads=F,show_xlab=T,show_ylab=T,show_yaxslab=T,show_legend=T,plt_par,plot_samp,
              plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
              sgefs_samp_seed,shefs_samp_seed,event,outfile)

plt_cumul_ensembles(nsamps,nsets,lds,site,forc='sgefs',show_leads=F,show_xlab=T,show_ylab=F,show_yaxslab=T,show_legend=F,plt_par,plot_samp,
                    plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                    sgefs_samp_seed,shefs_samp_seed,event,outfile)

dev.off()




#------------------------------------------------------------------------
#plot figure 5
source('plot/plot_cumul_rank_hist_hopper.r')
source('plot/plot_ecrps_hopper.r')

nsamps<-10
nsets<-10
lds<-c(1,3,5,10)
lpcnt<-0
upcnt<-90
show_leads<-T
show_xlab<-T
plt_par<-F
#outfile<-'v10_wtsamp30'

png(filename=paste('plot/figs/fig5_',outfile,'.png',sep=''),width = 1024, height = 896)

par(mfrow=c(4,length(lds)),mar=c(3,4,1,0),mgp=c(2.5,0.75,0),tcl=-0.25,cex.axis=1.75,cex.lab=2,font.lab=2,oma=c(1,5,4,0.5))

plt_cumul_rank_hist(nsamps,nsets,lds,lpcnt,upcnt,show_leads=T,show_xlab=F,show_xaxs=F,show_legend=F,plt_par,outfile)
plt_cumul_rank_hist(nsamps,nsets,lds,lpcnt=90,upcnt=100,show_leads=F,show_xlab=T,show_xaxs=T,show_legend=T,plt_par,outfile)
plt_ecrps(nsamps,nsets,lds,lpcnt=0,upcnt=90,show_leads=F,show_xlab=F,plt_par,outfile,ylm=200,yinc=50)
plt_ecrps(nsamps,nsets,lds,lpcnt=90,upcnt=100,show_leads=F,show_xlab=T,plt_par,outfile,ylm=2000,yinc=1000)

dev.off()

#----------------------------------------------------------------
#plot figure 6
source('plot/plot_spread-skill_hopper.r')

nsamps<-10
nsets<-10
lds<-c(1,3,5,10)
lpcnt<-0
upcnt<-100
show_leads<-T
show_xlab<-F
plt_par<-T
plt_lim<-1500

png(filename=paste('plot/figs/fig6_',outfile,'-',plt_lim,'.png',sep=''),width = 1024, height = 256)

par(mfrow=c(1,length(lds)),mar=c(3,3,0,0),mgp=c(3,0.75,0),tcl=-0.25,cex.axis=2,cex.lab=2.25,font.lab=2,oma=c(0.5,3.5,3.5,0.5))

plt_sprd_skill(nsamps,nsets,lds,lpcnt,upcnt,show_leads,show_xlab,plt_par,outfile,plt_lim)

dev.off()

#zoomed in plot
plt_lim <- 500

png(filename=paste('plot/figs/fig5_',outfile,'-',plt_lim,'.png',sep=''),width = 1024, height = 256)

par(mfrow=c(1,length(lds)),mar=c(3,3,0,0),mgp=c(2,0.5,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3.5,3.5,0.5))

plt_sprd_skill(nsamps,nsets,lds,lpcnt,upcnt,show_leads,show_xlab,plt_par,outfile,plt_lim)

dev.off()

#rm(list=ls());gc()

#-------------------------------------------------------------------------
#plot figure 7
source('plot/plot_top-10-events_hopper.r')

nsamps<-10
nsets<-10
lds<-c(1,3,5,10)
show_xlab<-F
show_extremes<-T
show_mn_xlab<-F
show_var_xlab<-T
calc_shefs<-F
syn_forc_typ<-'shefs'
outfile_hefs<-outfile
outfile_gefs<-'v3'
mn_lim<-3
vr_lim<-4
max_idx<-c("2005-12-31", "1995-01-09", "1986-02-18", "1993-01-21", "1997-01-02",
            "2004-02-17", "1998-02-20","1995-03-10", "2002-12-21","2006-03-06")

png(filename=paste('plot/figs/fig7_',outfile,'.png',sep=''),width = 1024, height = 768)

par(mfrow=c(2,1),mar=c(2,3,0,1),mgp=c(2,0.75,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3,0.5,0.5))

plt_top10_events(nsamps,nsets,lds,show_xlab,show_extremes,show_mn_xlab,show_var_xlab,calc_shefs,outfile_gefs,outfile_hefs,syn_forc_typ,mn_lim,vr_lim,max_idx)

dev.off()


#plot for GEFS
png(filename=paste('plot/figs/figs3a_sgefs_',outfile_gefs,'.png',sep=''),width = 1024, height = 768)

par(mfrow=c(2,1),mar=c(2,3,0,1),mgp=c(2,0.75,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3,0.5,0.5))

plt_top10_events(nsamps,nsets=1,lds,show_xlab,show_extremes,show_mn_xlab,show_var_xlab,calc_shefs=T,outfile_gefs,outfile_hefs,syn_forc_typ='sgefs',mn_lim,vr_lim,max_idx)

dev.off()


#-------------------------------------------------------------------------
#plot figure S3
source('plot/plot_top-100-events_hopper.r')

nsamps<-10
nsets<-10
lds<-c(1,3,5,10)
show_xlab<-F
show_extremes<-F
show_mn_xlab<-F
show_var_xlab<-T
calc_shefs<-F
syn_forc_typ<-'shefs'
outfile_hefs<-outfile
outfile_gefs<-'v3'
mn_lim<-3
vr_lim<-4

png(filename=paste('plot/figs/figs3b_',outfile,'.png',sep=''),width = 1024, height = 768)

par(mfrow=c(2,1),mar=c(2,3,0,1),mgp=c(2,0.75,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3,0.5,0.5))

plt_top100_events(nsamps,nsets,lds,show_xlab,show_extremes,show_mn_xlab,show_var_xlab,calc_shefs,outfile_gefs,outfile_hefs,syn_forc_typ,mn_lim,vr_lim)

dev.off()


#plot for GEFS
png(filename=paste('plot/figs/figs3c_sgefs_',outfile_gefs,'.png',sep=''),width = 1024, height = 768)

par(mfrow=c(2,1),mar=c(2,3,0,1),mgp=c(2,0.75,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3,0.5,0.5))

plt_top100_events(nsamps,nsets=1,lds,show_xlab,show_extremes,show_mn_xlab,show_var_xlab,calc_shefs=T,outfile_gefs,outfile_hefs,syn_forc_typ='sgefs',mn_lim,vr_lim)

dev.off()

#------------------------------------------------------------------------------------------------
#Fig 10
source('plot/plot_ensembles_hopper.r')
source('plot/plot_cumul_ensembles_hopper.r')

nsets<-10
nsamps<-10
lds<-c(5)
site<-'lamc'
forc<-'hefs'
show_leads<-T
show_legend<-T
leg_pos<-'bottomright'
show_xlab<-T
show_ylab<-T
show_yaxslab<-T
plt_par<-T
plot_samp<-F
plot_median<-F
plot_mean<-T
bounds<-c(0.05,0.95)
sgefs_seed<-1
shefs_set_seed<-1
shefs_seed<-1
hefs_samp_seed<-1
sgefs_samp_seed<-1
shefs_samp_seed<-1
event<-'1986-02-18'
#outfile<-'v9_wtsamp30'

png(filename=paste('plot/figs/fig10_',outfile,'.png',sep=''),width = 1024, height = 512)

par(mfrow=c(2,3),mar=c(2,4,0,1),mgp=c(2.5,0.75,0),tcl=-0.25,cex.axis=1.75,cex.lab=2,font.lab=2,oma=c(2.5,5,4,0.5))

for(l in c(8,7,6)){
  lds<-l
  if(l==8){
    plt_ensembles(nsamps,nsets,lds,site,forc,show_leads,show_xlab=F,show_ylab=T,show_yaxslab=T,show_legend=F,plt_par,plot_samp,
              plot_median,plot_mean=T,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
              sgefs_samp_seed,shefs_samp_seed,event,outfile)
  }
  if(l!=8){
    plt_ensembles(nsamps,nsets,lds,site,forc,show_leads,show_xlab=F,show_ylab=F,show_yaxslab=F,show_legend=F,plt_par,plot_samp,
                  plot_median,plot_mean=T,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                  sgefs_samp_seed,shefs_samp_seed,event,outfile)
  }
}

for(l in c(8,7,6)){
  lds<-l
  if(l==8){
    plt_cumul_ensembles(nsamps,nsets,lds,site,forc,show_leads=F,show_xlab=T,show_ylab=T,show_yaxslab=T,show_legend=F,leg_pos,plt_par,plot_samp,
                    plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                    sgefs_samp_seed,shefs_samp_seed,event,outfile)
  }
  if(l==7){
    plt_cumul_ensembles(nsamps,nsets,lds,site,forc,show_leads=F,show_xlab=T,show_ylab=F,show_yaxslab=F,show_legend=F,leg_pos,plt_par,plot_samp,
                        plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                        sgefs_samp_seed,shefs_samp_seed,event,outfile)
  }
  if(l==6){
    plt_cumul_ensembles(nsamps,nsets,lds,site,forc,show_leads=F,show_xlab=T,show_ylab=F,show_yaxslab=F,show_legend=T,leg_pos,plt_par,plot_samp,
                        plot_median,plot_mean,bounds,sgefs_seed,shefs_set_seed,shefs_seed,hefs_samp_seed,
                        sgefs_samp_seed,shefs_samp_seed,event,outfile)
  }
}


dev.off()


#------------------------------------------------------------------------------------------------------
#rank histogram plots
par(mfrow=c(2,length(lds)),mar=c(3,4,0,1),mgp=c(2,0.5,0),tcl=-0.25,cex.axis=1.5,cex.lab=2,font.lab=2,oma=c(0.5,3.5,3.5,0.5))

plt_rank_hist(m,lds,lpcnt,upcnt,show_leads,syn_plot,shefs_typ,plt_par,show_xlab)
plt_rank_hist(m,lds,lpcnt,upcnt,show_leads=F,syn_plot='sHEFS',shefs_typ,plt_par)


############################################END##########################################

