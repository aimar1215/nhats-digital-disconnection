if(.Platform$OS.type=='windows')try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
args<-commandArgs(TRUE);source(file.path(args[1],'code/common_repair.R'),encoding='UTF-8')
res<-read.csv(file.path(ROOT,'tables/models_repaired.csv'),check.names=FALSE)
png(file.path(ROOT,'figures/FigureS1_adjustment.png'),width=2800,height=1700,res=300,type='cairo',family='Times New Roman')
par(mfrow=c(1,2),mar=c(5,6,3,1));for(oc in c('composite','death')){
 x<-res[res$analysis==paste0(oc,':stacked')&res$weight=='w_follow'&grepl('^offline',res$term),];x<-x[match(paste0('M',0:4),x$model),]
 plot(x$HR,5:1,xlim=c(.7,3.2),ylim=c(.5,5.5),yaxt='n',xlab='Hazard ratio (95% CI)',ylab='',pch=16,main=ifelse(oc=='composite','Composite outcome','Mortality'))
 axis(2,5:1,labels=paste0('M',0:4),las=1);segments(x$lo,5:1,x$hi,5:1,lwd=1.5);abline(v=1,lty=2,col='grey50')
}
dev.off()
# R2-14: each entry cohort displayed separately; no independent-cohort meta-analysis.
curve<-list();nr<-list();set.seed(20261028)
for(b in 11:13){
 z<-as.data.frame(filter(d,cohort==b));K<-14-b;des<-svydesign(ids=~vunit,strata=~vstrat,weights=~w0,data=z,nest=TRUE);rd<-as.svrepdesign(des,type='bootstrap',replicates=500)
 calc<-function(w,data){v<-c();for(g in levels(data$digi3)){ii<-data$digi3==g;s<-1;for(k in 1:K){den<-sum(w[ii&data$comp_time>=k]);num<-sum(w[ii&data$comp_time==k&data$evt_comp==1]);if(den>0)s<-s*(1-num/den);v<-c(v,1-s)}};v}
 rr<-withReplicates(rd,calc);ci<-qt(.975,degf(des));j<-0
 for(g in levels(z$digi3)){for(k in 1:K){j<-j+1;curve[[length(curve)+1]]<-data.frame(cohort=b,group=g,k=k,risk=coef(rr)[j],lo=max(0,coef(rr)[j]-ci*SE(rr)[j]),hi=min(1,coef(rr)[j]+ci*SE(rr)[j]))};nr[[length(nr)+1]]<-data.frame(cohort=b,group=g,k=0:K,n=sapply(0:K,function(k)sum(z$digi3==g&z$comp_time>=k)))}
}
curve<-bind_rows(curve);nr<-bind_rows(nr);put(curve,'R2_14_curves_by_entry');put(nr,'R2_14_risk_tables_by_entry')
png(file.path(ROOT,'figures/Figure3_all_entry_cohorts.png'),width=3600,height=1800,res=300,type='cairo',family='Times New Roman')
par(mfrow=c(1,3),mar=c(8,4.5,3,1));cols<-c('#111111','#666666','#999999');labs<-c('Connected non-user','Telehealth user','Offline')
for(b in 11:13){K<-14-b;plot(c(0,K+.08),c(0,.4),type='n',xlab='Follow-up round',ylab='Cumulative composite risk',xaxt='n',main=paste('R',b,'entry',sep=' '));axis(1,0:K)
 for(i in 1:3){x<-filter(curve,cohort==b,group==levels(d$digi3)[i]);xx<-rep(0:K,each=2)[-1];yylo<-rep(c(0,x$lo),each=2)[-2*(K+1)];yyhi<-rep(c(0,x$hi),each=2)[-2*(K+1)];polygon(c(xx,rev(xx)),c(yylo,rev(yyhi)),col=adjustcolor(cols[i],.12),border=NA);lines(0:K,c(0,x$risk),type='s',lty=i,lwd=2,col=cols[i]);arrows((1:K)+(i-2)*.015,x$lo,(1:K)+(i-2)*.015,x$hi,angle=90,code=3,length=.025,col=cols[i]);counts<-filter(nr,cohort==b,group==levels(d$digi3)[i]);mtext(paste(labs[i],paste(counts$n,collapse=' / ')),side=1,line=3.5+i,cex=.62)}
 legend('topleft',labs,col=cols,lty=1:3,lwd=2,bty='n',cex=.7)
}
dev.off();cat('FIGURES_DONE\n')
