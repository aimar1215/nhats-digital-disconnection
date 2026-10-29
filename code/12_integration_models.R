Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
a<-commandArgs(TRUE);source(file.path(a[1],'code/common_repair.R'),encoding='UTF-8');OUT<-if(length(a)>1)a[2] else a[1]
put<-function(x,n)write.csv(x,file.path(OUT,'tables',paste0(n,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
s<-list()
for(oc in names(pp))for(type in c('stacked','first_entry')){
 p<-pp[[oc]];if(type=='first_entry')p<-semi_join(p,first,by=c('spid','cohort'))
 for(sub in c('exclude_first_interval','self_respondents','no_established_cvd','complete_M4')){
  ps<-switch(sub,exclude_first_interval=filter(p,t>=2),self_respondents=filter(p,proxy_f=='自答'),no_established_cvd=p,complete_M4=p)
  if(sub=='complete_M4')for(v in cv4)ps<-ps[!is.na(ps[[v]]) & !grepl('未知|Unknown',as.character(ps[[v]])),]
  r<-fit(ps,'M4','w0',omit=if(sub=='no_established_cvd')'estcvd_f' else NULL);attr(r,'fit')<-NULL;r$analysis<-paste(oc,type,sub,sep=':');s[[length(s)+1]]<-r
 }
 for(k in 1:3){r<-fit(filter(p,t==k),'M4','w0');attr(r,'fit')<-NULL;r$analysis<-paste(oc,type,'interval',k,sep=':');s[[length(s)+1]]<-r}
}
for(oc in names(pp))for(b in 11:13){r<-fit(filter(pp[[oc]],cohort==b),'M4','w0');attr(r,'fit')<-NULL;r$analysis<-paste(oc,'cohort',b,sep=':');s[[length(s)+1]]<-r}
put(bind_rows(s),'baseline_weight_sensitivities')
# Same residential endpoint classification, M4 and survey design as primary models.
z<-d %>% mutate(cr_status=case_when(!is.na(nh_round)&(is.na(death_round)|nh_round<=death_round)~1L,!is.na(death_round)&death_round>cohort&death_round<=14~2L,TRUE~0L),cr=factor(cr_status,0:2,c('censor','NH','death')))
np<-pp$composite %>% left_join(select(z,spid,cohort,cr_status),by=c('spid','cohort')) %>% mutate(event=event*as.integer(cr_status==1))
cs<-fit(np,'M4','w0');attr(cs,'fit')<-NULL;cs$analysis<-'residential:cause_specific:M4'
put(cs,'residential_M4_cause_specific')
# Retain actual discrete round times; expansion weights estimated within entry cohort.
crdat<-filter(z,fu_rounds>=1);crdat<-crdat[complete.cases(crdat[,c('offline',cv4)]),]
fg<-finegray(as.formula(paste('Surv(fu_rounds,cr)~',paste(c('offline',cv4,'cohort_f','spid','vunit','vstrat','w0','strata(cohort)'),collapse='+'))),data=crdat,etype='NH',weights=w0)
# finegray with weights=w0 returns fgwt including case weight; do not multiply twice.
fgdes<-svydesign(ids=~vunit+spid,strata=~vstrat,weights=~fgwt,data=fg,nest=TRUE)
# Compute the same Taylor covariance as survey::svycoxph without its omnibus
# Wald inversion (singular joint covariance does not require inversion for the
# prespecified single exposure contrast). No coefficients/covariates are removed.
fm<-coxph(as.formula(paste('Surv(fgstart,fgstop,fgstatus)~offline+cohort_f+',paste(cv4,collapse='+'))),data=fg,weights=fgwt/mean(fgwt),model=TRUE)
db<-as.matrix(residuals(fm,type='dfbeta',weighted=TRUE))
fm$var<-svyrecvar(db,fgdes$cluster,fgdes$strata,fgdes$fpc,postStrata=fgdes$postStrata)
stopifnot(is.finite(coef(fm)['offline断连']),all(is.finite(fm$var)))
put(data.frame(term=names(coef(fm)),beta=coef(fm),se=sqrt(diag(vcov(fm)))),'FineGray_coefficient_diagnostics')

b<-coef(fm)['offline断连'];se<-sqrt(vcov(fm)['offline断连','offline断连']);df<-degf(fgdes)-sum(!is.na(coef(fm)))+1
put(data.frame(HR=exp(b),lo=exp(b-qt(.975,df)*se),hi=exp(b+qt(.975,df)*se),beta=b,se=se,df=df,n_baseline=nrow(crdat),n_person=n_distinct(crdat$spid),method='M4 survey Fine-Gray; entry-stratified weighted expansion; round-time ties retained; SE conditional on expansion weights'),'residential_M4_FineGray')
capture.output(sessionInfo(),file=file.path(OUT,'documentation/R_session_info.txt'))
cat('INTEGRATION_ANALYSIS_DONE\n')
