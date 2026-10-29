if(.Platform$OS.type=="windows")try(Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8"),silent=TRUE)
args<-commandArgs(TRUE);source(file.path(args[1],'code/common_repair.R'),encoding='UTF-8')
suppressPackageStartupMessages(library(mice))
# One imputation record per participant: do not treat repeated entries as independent.
mi_base<-as.data.frame(first)
mi_base[]<-lapply(mi_base,function(x)if(inherits(x,"haven_labelled"))as.numeric(x) else x)
mi_base$logw<-log(mi_base$w0);mi_base$stratum<-factor(mi_base$vstrat)
mi_base$income_log<-log1p(rowMeans(mi_base[paste0('income',1:20)],na.rm=TRUE))
mi_base$income_log[!is.finite(mi_base$income_log)]<-NA_real_
# Item nonresponse versus structural skips retained as explicit Unknown levels.
a<-as.matrix(mi_base[c('sc_eat','sc_bath','sc_toil','sc_dres')])
mi_base$adl_cat[mi_base$adl_cat=='Unknown' & rowSums(a==-1,na.rm=TRUE)==0]<-NA
mi_base$dep_f[grepl('未知',mi_base$dep_f)& !(mi_base$dep1==-1|mi_base$dep2==-1)]<-NA
mi_base$mem_f[grepl('未知',mi_base$mem_f)&mi_base$ratemem!=-1]<-NA
mi_base$out3[grepl('未知',mi_base$out3)&mi_base$outoft!=-1]<-NA
covs<-lapply(covs,function(x)ifelse(x=='income_cat','income_log',x))
impvars<-unique(c(covs$M4,'offline','cohort_f','logw','stratum','evt_comp','death_event','comp_time','death_time','proxy_f'))
allres<-list();diagnostics<-list();completed<-list();counts<-list()
for(k in 1:20){
 cat('IMPUTATION',k,'of 20\n');flush.console()
 x<-mi_base[impvars];x$income_log<-log1p(mi_base[[paste0('income',k)]])
 meth<-mice::make.method(x);meth[]<-''
 for(v in names(x)[colSums(is.na(x))>0])meth[v]<-if(is.factor(x[[v]]))'cart' else 'pmm'
 pred<-mice::make.predictorMatrix(x);diag(pred)<-0
 im<-mice(x,m=1,maxit=10,method=meth,predictorMatrix=pred,seed=20261028+k,printFlag=FALSE)
 z<-complete(im,1);counts[[k]]<-data.frame(imputation=k,missing_remaining=sum(is.na(z)),logged_events=if(is.null(im$loggedEvents))0 else nrow(im$loggedEvents))
 saveRDS(im,file.path(ROOT,'data',paste0('mice_first_entry_',k,'.rds')))
 completed[[k]]<-cbind(first[c('spid','cohort')],z)
 for(oc in c('composite','death')){
  p<-if(oc=='composite')dat$pp_comp else dat$pp_death
  p<-semi_join(p,first,by=c('spid','cohort'));jj<-match(p$spid,first$spid)
  for(v in covs$M4)p[[v]]<-z[[v]][jj]
  for(w in c('w0','w_follow'))for(mm in c('M3','M4')){
   r<-fit(p,mm,w);attr(r,'fit')<-NULL;r$endpoint<-oc;r$imputation<-k;allres[[length(allres)+1]]<-r
  }
 }
}
r<-bind_rows(allres)
pool<-r %>% group_by(endpoint,model,weight,term) %>% summarise(m=n(),B=var(beta),W=mean(se^2),beta=mean(beta),dfcom=min(df),n_pp=min(n_pp),n_person=min(n_person),.groups='drop') %>% mutate(Tvar=W+(1+1/m)*B,lambda=(1+1/m)*B/Tvar,df_old=ifelse(lambda>0,(m-1)/lambda^2,Inf),df_obs=(dfcom+1)/(dfcom+3)*dfcom*(1-lambda),df=1/(1/df_old+1/df_obs),se=sqrt(Tvar),HR=exp(beta),lo=exp(beta-qt(.975,df)*se),hi=exp(beta+qt(.975,df)*se))
put(r,'MI_individual_estimates');put(pool,'MI_first_entry_pooled');put(bind_rows(counts),'MI_diagnostics');saveRDS(completed,file.path(ROOT,'data/MI_completed_first_entry.rds'))
cat('MI_DONE\n');print(pool)
