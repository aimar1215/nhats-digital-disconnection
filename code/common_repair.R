if(.Platform$OS.type=='windows')try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
suppressPackageStartupMessages({library(haven);library(dplyr);library(tidyr);library(survey);library(survival)})
options(survey.lonely.psu='adjust')
ROOT<-commandArgs(TRUE)[1];dat<-readRDS(file.path(ROOT,'data/repaired_data.rds'));d<-dat$coh;first<-dat$first
put<-function(x,n)write.csv(x,file.path(ROOT,'tables',paste0(n,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
cv1<-c('agec','sex','race','educ3');cv2<-c(cv1,'medicaid_f','livealone','married','income_cat')
cv3<-c(cv2,'como','srh_f','dementia_f','estcvd_f');cv4<-c(cv3,'adl_cat','dep_f','mem_f','out3')
covs<-list(M0=character(),M1=cv1,M2=cv2,M3=cv3,M4=cv4)
fit<-function(p,mm='M4',weight='w0',expo='offline',extra=NULL,omit=NULL){
 cv<-setdiff(covs[[mm]],omit);tt<-c(expo,if(n_distinct(p$cohort)>1)'cohort_f',if(n_distinct(p$t)>1)'tf',cv,extra)
 p<-p[is.finite(p[[weight]])&p[[weight]]>0,]; p<-as.data.frame(p);rownames(p)<-NULL;p$ww<-p[[weight]]
 des<-svydesign(ids=~vunit+spid,strata=~vstrat,weights=~ww,data=p,nest=TRUE)
 m<-svyglm(as.formula(paste('event~',paste(tt,collapse='+'))),design=des,family=quasibinomial('cloglog'))
 co<-coef(m);se<-sqrt(diag(vcov(m)));ix<-grep(paste0('^',expo),names(co));df<-m$df.residual
 ans<-data.frame(term=names(co)[ix],beta=co[ix],se=se[ix],HR=exp(co[ix]),lo=exp(co[ix]-qt(.975,df)*se[ix]),hi=exp(co[ix]+qt(.975,df)*se[ix]),df=df,n_pp=nrow(model.frame(m)),n_person=n_distinct(p$spid[as.integer(rownames(model.frame(m)))]),model=mm,weight=weight,row.names=NULL)
 attr(ans,'fit')<-m;ans
}
# Add current-round observed assessment flag for contiguous-follow-up sensitivity.
addobs<-function(p){
 f<-dat$follow;j<-match(paste(p$spid,p$round_end),paste(f$spid,f$round));p$obs<-as.integer(f$res[j] %in% 1:6|p$event==1)
 p %>% arrange(spid,cohort,t) %>% group_by(spid,cohort) %>% mutate(contiguous=cumprod(obs)) %>% ungroup()
}
pp<-list(composite=addobs(dat$pp_comp),death=addobs(dat$pp_death))
