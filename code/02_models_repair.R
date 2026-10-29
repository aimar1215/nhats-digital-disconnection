if(.Platform$OS.type=="windows")try(Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8"),silent=TRUE)
args<-commandArgs(TRUE);source(file.path(args[1],'code/common_repair.R'),encoding='UTF-8')
res<-list();sens<-list();errors<-list()
safe<-function(p,label,mm='M4',weight='w0',expo='offline',omit=NULL){
 tryCatch({r<-fit(p,mm,weight,expo,omit=omit);attr(r,'fit')<-NULL;r$analysis<-label;r},error=function(e){errors[[length(errors)+1]]<<-data.frame(analysis=label,error=conditionMessage(e));NULL})
}
for(oc in names(pp)){
 for(type in c('stacked','first_entry')){
  p<-pp[[oc]];if(type=='first_entry')p<-semi_join(p,first,by=c('spid','cohort'))
  for(w in c('w0','w_follow'))for(mm in names(covs))for(expo in c('offline','digi3')){
   cat('FIT',oc,type,w,mm,expo,'\n');flush.console()
   r<-safe(p,paste(oc,type,sep=':'),mm,w,expo);res[[length(res)+1]]<-r
  }
  for(sub in c('exclude_first_interval','self_respondents','contiguous_followup','no_established_cvd','complete_M4')){
   ps<-switch(sub,exclude_first_interval=filter(p,t>=2),self_respondents=filter(p,proxy_f=='自答'),contiguous_followup=filter(p,contiguous==1),no_established_cvd=p,complete_M4=p)
   if(sub=='complete_M4')for(v in cv4)ps<-ps[!is.na(ps[[v]]) & !grepl('未知|Unknown',as.character(ps[[v]])),]
   r<-safe(ps,paste(oc,type,sub,sep=':'),'M4','w_follow',omit=if(sub=='no_established_cvd')'estcvd_f' else NULL);sens[[length(sens)+1]]<-r
  }
  for(k in 1:3){r<-safe(filter(p,t==k),paste(oc,type,'interval',k,sep=':'),'M4','w_follow');sens[[length(sens)+1]]<-r}
 }
 for(b in 11:13){r<-safe(filter(pp[[oc]],cohort==b),paste(oc,'cohort',b,sep=':'),'M4','w_follow');sens[[length(sens)+1]]<-r}
}
put(bind_rows(res),'models_repaired');put(bind_rows(sens),'sensitivity_repaired');put(bind_rows(errors),'model_errors')
# Variable-by-variable missingness, by entry, including unknown categories.
vars<-unique(c('digi3',cv4,'income_obs'))
missing<-bind_rows(lapply(11:13,function(b)bind_rows(lapply(vars,function(v){x<-d[d$cohort==b,v];data.frame(cohort=b,variable=v,N=length(x),NA_n=sum(is.na(x)),unknown_n=sum(grepl('未知|Unknown',as.character(x)),na.rm=TRUE),missing_pct=100*mean(is.na(x)|grepl('未知|Unknown',as.character(x))))}))))
put(missing,'missingness_per_variable')
# Current vs entry weights coverage is reported, not silently dropped.
put(bind_rows(lapply(names(pp),function(oc) pp[[oc]] %>% group_by(cohort) %>% summarise(endpoint=oc,N=n(),nonpositive_or_missing=sum(!is.finite(w_follow)|w_follow<=0),events_without_positive_weight=sum(event==1&(!is.finite(w_follow)|w_follow<=0))))),'followup_weight_coverage')
# Table 1: earliest eligible entry, explicit missing rows; weighted percentages within group.
rows<-list()
for(v in cv4){for(g in levels(first$digi3)){
 z<-first[first$digi3==g,];x<-as.character(z[[v]]);x[is.na(x)]<-'Missing'
 if(v=='como'){rows[[length(rows)+1]]<-data.frame(variable=v,level='Mean',group=g,n=sum(!is.na(z[[v]])),weighted_percent_or_mean=weighted.mean(z[[v]],z$w0,na.rm=TRUE));next}
 for(l in sort(unique(x)))rows[[length(rows)+1]]<-data.frame(variable=v,level=l,group=g,n=sum(x==l),weighted_percent_or_mean=100*sum(z$w0[x==l])/sum(z$w0))
}}
put(bind_rows(rows),'baseline_first_entry')
cat('MODEL_BATCH_DONE\n')
