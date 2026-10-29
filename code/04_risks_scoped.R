if(.Platform$OS.type=='windows')try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
args<-commandArgs(TRUE);source(file.path(args[1],'code/common_repair.R'),encoding='UTF-8')
# R2-7 and R2-14: R11 target; weighted discrete-round KM handles right censoring.
b<-first %>% filter(cohort==11);b<-as.data.frame(b)
des<-svydesign(ids=~vunit,strata=~vstrat,weights=~w0,data=b,nest=TRUE)
set.seed(20261028);repdes<-as.svrepdesign(des,type='bootstrap',replicates=500)
km<-function(w,time,event,k){s<-1;for(j in 1:k){risk<-sum(w[time>=j]);ev<-sum(w[time==j&event==1]);if(risk>0)s<-s*(1-ev/risk)};1-s}
stat<-function(w,data){
 r<-c()
 for(k in 1:3){
  for(g in levels(data$digi3)){ix<-data$digi3==g;r<-c(r,setNames(km(w[ix],data$comp_time[ix],data$evt_comp[ix],k),paste('group',g,k,sep=':')))}
  a<-data$digi3=='完全断连';off<-km(w[a],data$comp_time[a],data$evt_comp[a],k);con<-km(w[!a],data$comp_time[!a],data$evt_comp[!a],k)
  r<-c(r,setNames(c(off,con,off-con),paste(c('offline','all_connected','risk_difference'),k,sep=':')))
 };r
}
r<-withReplicates(repdes,stat);cc<-coef(r);ss<-SE(r);tt<-qt(.975,degf(des))
rr<-data.frame(estimand=names(cc),estimate=as.numeric(cc),se=as.numeric(ss),lo=as.numeric(cc-tt*ss),hi=as.numeric(cc+tt*ss))
put(rr,'R11_KM_risks_bootstrap');put(b %>% group_by(digi3) %>% summarise(N=n(),risk_round1=sum(comp_time>=1),risk_round2=sum(comp_time>=2),risk_round3=sum(comp_time>=3)),'R11_numbers_at_risk')
# R2-7: M4 conditional model standardization, same R11 covariate distribution for both exposure settings.
p<-filter(dat$pp_comp,cohort==11);z<-fit(p,'M4','w0');m<-attr(z,'fit');cv<-covs$M4
bb<-b[complete.cases(b[cv]),];bb$cohort_f<-droplevels(bb$cohort_f)
coefm<-coef(m);V<-vcov(m);if(anyNA(coefm))stop('Aliased coefficients: standardized risks require model review')
standard<-function(k,level){
 H<-rep(0,nrow(bb));grad<-matrix(0,nrow(bb),length(coefm));colnames(grad)<-names(coefm)
 for(j in 1:k){
  nd<-bb;nd$offline<-factor(level,levels=levels(p$offline));nd$tf<-factor(j,levels=levels(p$tf))
  X<-model.matrix(delete.response(terms(m)),nd,contrasts.arg=m$contrasts,xlev=m$xlevels);X<-X[,names(coefm),drop=FALSE]
  ee<-exp(drop(X%*%coefm));H<-H+ee;grad<-grad+X*ee
 }
 risk<-1-exp(-H);g<-colSums(grad*exp(-H)*bb$w0)/sum(bb$w0)
 list(risk=weighted.mean(risk,bb$w0),g=g)
}
rs<-list()
for(k in 1:3){a<-standard(k,levels(p$offline)[2]);c<-standard(k,levels(p$offline)[1]);for(n in c('offline','all_connected','risk_difference')){
 r<-switch(n,offline=a$risk,all_connected=c$risk,risk_difference=a$risk-c$risk);g<-switch(n,offline=a$g,all_connected=c$g,risk_difference=a$g-c$g);se<-sqrt(drop(t(g)%*%V%*%g));q<-qt(.975,m$df.residual)
 rs[[length(rs)+1]]<-data.frame(round=k,estimand=n,risk=r,se=se,lo=r-q*se,hi=r+q*se,N_standardized=nrow(bb),CI_method='design coefficient delta; conditional on observed target covariates')
}}
put(bind_rows(rs),'R11_M4_standardized_risks')
# Publication-style plot uses original study group labels in English.
png(file.path(ROOT,'figures/Figure3_R11_cumulative_risk.png'),width=2400,height=2000,res=300,type='cairo',family='Times New Roman')
par(mar=c(7,5,2,1));plot(0:3,rep(0,4),type='n',xlab='Follow-up round since R11 entry',ylab='Cumulative composite risk',ylim=c(0,.4),xaxt='n');axis(1,0:3)
cols<-c('#222222','#777777','#AAAAAA');labs<-c('Connected non-user','Telehealth user','Offline')
for(i in 1:3){g<-levels(b$digi3)[i];x<-rr[match(paste('group',g,1:3,sep=':'),rr$estimand),];polygon(c(0:3,3:0),c(0,pmax(0,x$lo),rev(c(0,pmin(1,x$hi)))) ,col=adjustcolor(cols[i],alpha.f=.15),border=NA);lines(0:3,c(0,x$estimate),type='s',col=cols[i],lwd=2,lty=i)
 nr<-sapply(0:3,function(k)sum(b$digi3==g&b$comp_time>=k));mtext(paste(labs[i],paste(nr,collapse=' / ')),side=1,line=3+i,cex=.7)
}
legend('topleft',labs,col=cols,lty=1:3,lwd=2,bty='n');dev.off()
cat('RISK_DONE\n');print(bind_rows(rs))
