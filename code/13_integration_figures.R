Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
a<-commandArgs(TRUE);S<-a[1];D<-if(length(a)>1)a[2] else a[1]
read<-function(n,root=S)read.csv(file.path(root,'tables',n),fileEncoding='UTF-8')

m<-read('models_repaired.csv')
makefigdata<-function(secondary=FALSE){rows<-list();for(type in if(secondary)c('telehealth','offline') else c('stacked','first_entry'))for(oc in c('composite','death')){
term<-if(!secondary)'offline断连' else if(type=='telehealth')'digi3远程使用者' else 'digi3完全断连'
x<-m[m$analysis==paste0(oc,':',if(secondary)'stacked' else type)&m$model=='M4'&m$weight=='w0'&m$term==term,]
rows[[length(rows)+1]]<-data.frame(label=paste(oc,'-',gsub('_',' ',type)),HR=x$HR,lo=x$lo,hi=x$hi,display=sprintf('%.2f (%.2f-%.2f)',x$HR,x$lo,x$hi))};do.call(rbind,rows)}
write.csv(makefigdata(),file.path(D,'tables/Figure2_data.csv'),row.names=FALSE)
write.csv(makefigdata(TRUE),file.path(D,'tables/Figure4_data.csv'),row.names=FALSE)

raw<-read('R2_11_raw_eligibility_flow.csv');flow<-read('R2_11_flow_by_entry.csv')
png(file.path(D,'figures/Figure1_flow.png'),width=1474,height=1606,res=220,type='cairo',family='Times New Roman')
par(mar=c(1,1,2,1),mfrow=c(1,3),family='Times New Roman')
for(i in 1:3){r<-raw[i,];f<-flow[i,];plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));title(paste('Round',r$cohort,'entry'),cex.main=1.15)
 labs<-c(sprintf('Sample-person records\nn = %s',r$SP_records),sprintf('Exclude no qualifying condition\nn = %s\nExclude non-community: %s',r$excluded_not_CV,r$excluded_noncommunity),sprintf('Eligible community dwelling\nn = %s',f$eligible),sprintf('Exposure unavailable\nn = %s',f$missing_exposure),sprintf('Exposure observed\nn = %s',f$known_exposure),sprintf('No follow-up: %s\nM4 covariate exclusions: %s',f$no_followup,f$M4_covariate_excluded),sprintf('M4 entry-weight analysis\nn = %s',f$M4_with_entry_weight))
 for(j in seq_along(labs)){y<-.94-(j-1)*.135;rect(.025,y-.048,.975,y+.048,border='grey45');text(.5,y,labs[j],cex=.875);if(j<7)arrows(.5,y-.05,.5,y-.083,length=.07)}
};dev.off()
forest<-function(n,title,out){z<-read(n,D);png(file.path(D,'figures',out),width=1474,height=800,res=220,type='cairo',family='Times New Roman');par(mar=c(4,10,3,1),family='Times New Roman');y<-rev(seq_len(nrow(z)));plot(z$HR,y,xlim=c(.45,3.5),ylim=c(.5,nrow(z)+.5),yaxt='n',ylab='',xlab='Hazard ratio (95% CI)',pch=16,bty='l',main=title,cex.main=.95,cex.axis=.875,cex.lab=.9);axis(2,at=y,labels=z$label,las=1,tick=FALSE,cex.axis=.875);abline(v=1,lty=2,col='grey50');arrows(z$lo,y,z$hi,y,angle=90,code=3,length=.04);text(2.55,y,z$display,adj=0,cex=.875);dev.off()}
forest('Figure2_data.csv','M4 associations with complete digital disconnection','Figure2_M4.png')
forest('Figure4_data.csv','Digital groups compared with connected non-users','Figure4_digital_groups.png')
# Publication-size panels use the same aggregate curve results, without reanalysis.
c<-read('R2_14_curves_by_entry.csv');nr<-read('R2_14_risk_tables_by_entry.csv');g<-unique(c$group);lab<-c('Connected non-user','Telehealth user','Offline');col<-c('black','grey45','grey65')
png(file.path(D,'figures/Figure3_all_entry_cohorts.png'),width=1474,height=1650,res=220,type='cairo',family='Times New Roman')
par(mfrow=c(3,1),mar=c(5.5,4.2,1.7,1),mgp=c(2.3,.6,0),family='Times New Roman',cex=.875)
for(b in 11:13){K<-14-b;plot(c(0,K+.08),c(0,.4),type='n',xlab='',ylab='Cumulative risk',xaxt='n',main=paste('Round',b,'entry'));axis(1,0:K)
 for(i in 1:3){x<-c[c$cohort==b & c$group==g[i],];xx<-rep(0:K,each=2)[-1];yl<-rep(c(0,x$lo),each=2)[-2*(K+1)];yh<-rep(c(0,x$hi),each=2)[-2*(K+1)];polygon(c(xx,rev(xx)),c(yl,rev(yh)),col=adjustcolor(col[i],.12),border=NA);lines(0:K,c(0,x$risk),type='s',lty=i,lwd=1.5,col=col[i]);arrows((1:K)+(i-2)*.012,x$lo,(1:K)+(i-2)*.012,x$hi,angle=90,code=3,length=.025,col=col[i]);ct<-nr[nr$cohort==b&nr$group==g[i],];mtext(paste(lab[i],paste(ct$n,collapse=' / ')),side=1,line=2+i*.8,cex=.77)}
 mtext(paste('Follow-up round; numbers at risk for rounds',paste(0:K,collapse='/')),side=1,line=1.7,cex=.77);legend('topleft',lab,col=col,lty=1:3,bty='n',cex=.8,horiz=TRUE)
};dev.off()
m<-read('models_repaired.csv')
png(file.path(D,'figures/FigureS1_adjustment.png'),width=1474,height=1000,res=220,type='cairo',family='Times New Roman')
par(mfrow=c(1,2),mar=c(4,3,2,1),family='Times New Roman',cex=.875)
for(oc in c('composite','death')){x<-m[m$analysis==paste0(oc,':stacked')&m$weight=='w0'&grepl('^offline',m$term),];x<-x[match(paste0('M',0:4),x$model),];plot(x$HR,5:1,xlim=c(.65,3.5),ylim=c(.5,5.5),yaxt='n',ylab='',xlab='Hazard ratio (95% CI)',pch=16,main=ifelse(oc=='composite','Composite outcome','Mortality'));axis(2,5:1,paste0('M',0:4),las=1);segments(x$lo,5:1,x$hi,5:1);abline(v=1,lty=2,col='grey50')};dev.off()
