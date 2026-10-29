Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8');library(haven);library(dplyr)
ROOT<-commandArgs(TRUE)[1];raw<-Sys.getenv('NHATS_RAW_DIR')
old<-readRDS(file.path(ROOT,'baseline/data_derived/analytic_cohorts_v2.rds'));rows<-list();why<-list()
for(b in 11:13){
 f<-list.files(raw,pattern=paste0('Round_',b,'_SP_File'),full.names=TRUE,recursive=TRUE)[1];x<-read_dta(f)
 g<-function(n)as.numeric(x[[n]])
 cv<-g(paste0('hc',b,'disescn2')) %in% c(1,7)|g(paste0('hc',b,'disescn3')) %in% c(1,7)|g(paste0('hc',b,'disescn6')) %in% c(1,7)
 res<-g(paste0('r',b,'dresid'));rows[[b]]<-data.frame(cohort=b,SP_records=nrow(x),excluded_not_CV=sum(!cv),excluded_noncommunity=sum(cv&!res %in% 1),eligible_community=sum(cv&res==1))
 z<-old[old$cohort==b & is.na(old$digi3),];on<-as.numeric(z$online);te<-as.numeric(z$teleh)
 why[[b]]<-data.frame(cohort=b,reason=c('Online response unavailable/other','Online yes, telehealth response unavailable','Other uncategorized exposure'),n=c(sum(!on %in% c(1,2)),sum(on==1,na.rm=TRUE),sum(on==2,na.rm=TRUE)))
}
write.csv(bind_rows(rows),file.path(ROOT,'tables/R2_11_raw_eligibility_flow.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(bind_rows(why),file.path(ROOT,'tables/R2_11_exposure_exclusion_reasons.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(data.frame(entry=rep(11:13,c(3,2,1)),followup=c(12:14,13:14,14),baseline_weight=rep(c('w11anfinwgt0','w12anfinwgt0','w13anfinwgt0'),c(3,2,1)),followup_weight=c('w12an2015wgt0','w13an2015wgt0','w14an2015wgt0','w13an2022wgt0','w14an2022wgt0','w14anfinwgt0')),file.path(ROOT,'tables/R2_12_weight_variable_map.csv'),row.names=FALSE)
