if(.Platform$OS.type=='windows')try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
suppressPackageStartupMessages({library(haven);library(dplyr)})
ROOT<-commandArgs(TRUE)[1]; SRC<-dirname(ROOT);RAW<-Sys.getenv('NHATS_RAW_DIR')
d<-readRDS(file.path(ROOT,'baseline/data_derived/analytic_cohorts_v2.rds'));p<-readRDS(file.path(ROOT,'baseline/data_derived/analytic_personperiod_v2.rds'))
d<-as.data.frame(d);d$income_cat<-d$inc3;d$income_obs<-ifelse(d$totinc>=0,d$totinc,NA_real_)
d$comp_time<-d$fu_rounds;d$death_time<-d$fu_rounds;d$death_event<-d$evt_death
follow<-list();sources<-list()
for(b in 11:14){
 fn<-list.files(RAW,pattern=paste0('Round_',b,'_SP_File'),full.names=TRUE)[1]
 x<-read_dta(fn);g<-function(n)if(n %in% names(x))as.numeric(x[[n]]) else rep(NA_real_,nrow(x))
 follow[[as.character(b)]]<-data.frame(spid=g('spid'),round=b,res=g(paste0('r',b,'dresid')),w_fin=g(paste0('w',b,'anfinwgt0')),w_2015=g(paste0('w',b,'an2015wgt0')),w_2022=g(paste0('w',b,'an2022wgt0')))
 if(b<=13){
  ff<-list.files(RAW,pattern=paste0('(R',b,'_Int_Inc_Imp_File|Round_',b,'_Int_Inc_Imp_File).*dta$'),recursive=TRUE,full.names=TRUE)[1]
  inc<-read_dta(ff);ii<-which(d$cohort==b);jj<-match(d$spid[ii],inc$spid)
  for(k in 1:20){v<-as.numeric(inc[[paste0('ia',b,'dtoincimi',k)]][jj]);d[ii,paste0('income',k)]<-ifelse(v>=0,v,NA_real_)}
  sources[[b]]<-data.frame(cohort=b,income_source=ff)
 }
}
f<-bind_rows(follow);d<-filter(d,!is.na(digi3));p<-filter(p,!is.na(digi3))
p$income_cat<-d$income_cat[match(paste(p$spid,p$cohort),paste(d$spid,d$cohort))];p$round_end<-p$cohort+p$t
fw<-f[match(paste(p$spid,p$round_end),paste(f$spid,f$round)),];p$w_follow<-ifelse(p$cohort==11,fw$w_2015,ifelse(p$cohort==12,fw$w_2022,fw$w_fin))
p1<-p;p1$event<-p1$evt_comp;p1$endpoint<-'composite';p2<-p;p2$event<-p2$evt_death;p2$endpoint<-'death'
first<-d %>% arrange(spid,cohort) %>% distinct(spid,.keep_all=TRUE)
saveRDS(list(coh=d,first=first,pp_comp=p1,pp_death=p2,follow=f),file.path(ROOT,'data/repaired_data.rds'))
write.csv(bind_rows(sources),file.path(ROOT,'documentation/income_sources.csv'),row.names=FALSE,fileEncoding='UTF-8')
cat('REVIEWER_SCOPED_DATA',nrow(d),nrow(first),nrow(p),'\n')
