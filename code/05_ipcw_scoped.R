if(.Platform$OS.type=='windows')try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
args<-commandArgs(TRUE);source(file.path(args[1],'code/common_repair.R'),encoding='UTF-8')
# R2-15: IPCW sensitivity, first eligible record per person. Baseline outcome definitions retained.
b<-as.data.frame(first)
# Explicit missing levels are confined to this censoring-weight sensitivity model.
for(v in covs$M4)if(is.factor(b[[v]])){x<-as.character(b[[v]]);x[is.na(x)]<-'Missing';b[[v]]<-factor(x)}
risk<-list()
for(k in 1:3){q<-b[b$cohort+k<=14 & b$comp_time>=k-1 & !(b$evt_comp==1 & b$comp_time==k-1),];q$t<-k;q$tf<-factor(k,levels=1:3);q$observed<-as.integer(q$comp_time>=k);risk[[k]]<-q}
risk<-bind_rows(risk)
des<-svydesign(ids=~vunit,strata=~vstrat,weights=~w0,data=risk,nest=TRUE)
num<-svyglm(observed~offline+cohort_f+tf,design=des,family=quasibinomial())
den<-svyglm(as.formula(paste('observed~offline+cohort_f+tf+',paste(covs$M4,collapse='+'))),design=des,family=quasibinomial())
risk$pn<-as.numeric(predict(num,type='response'));risk$pd<-as.numeric(predict(den,type='response'))
rawpd<-risk$pd;risk$pn<-pmax(.001,pmin(.999,risk$pn));risk$pd<-pmax(.001,pmin(.999,risk$pd))
risk<-risk %>% arrange(spid,t) %>% group_by(spid) %>% mutate(sw=cumprod(pn/pd)) %>% ungroup()
p<-filter(risk,observed==1);p$event<-as.integer(p$t==p$comp_time & p$evt_comp==1)
cut<-quantile(p$sw,c(.01,.99),na.rm=TRUE);p$ipcw_raw<-p$w0*p$sw;p$ipcw_trim<-p$w0*pmax(cut[1],pmin(cut[2],p$sw))
res<-list()
for(w in c('w0','ipcw_raw','ipcw_trim')){z<-fit(p,'M4',w);attr(z,'fit')<-NULL;z$analysis<-'R2-15:first_entry:composite';res[[w]]<-z}
put(bind_rows(res),'IPCW_sensitivity')
put(data.frame(metric=c('baseline_persons','zero_followup_baselines','censor_model_intervals','minimum_unclipped_observation_probability','clipped_observation_probabilities','SW_p01','SW_p99','SW_max','outcome_intervals','effective_n_weighted_intervals'),value=c(nrow(b),sum(b$comp_time==0),nrow(risk),min(rawpd),sum(rawpd<.001|rawpd>.999),cut,max(p$sw),nrow(p),sum(p$ipcw_trim)^2/sum(p$ipcw_trim^2))),'IPCW_diagnostics')
saveRDS(list(numerator=num,denominator=den,intervals=risk),file.path(ROOT,'data/IPCW_diagnostics.rds'))
cat('IPCW_DONE\n');print(bind_rows(res))
