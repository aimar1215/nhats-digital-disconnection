if(.Platform$OS.type=='windows')try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
args<-commandArgs(TRUE);source(file.path(args[1],'code/common_repair.R'),encoding='UTF-8')
old<-readRDS(file.path(ROOT,'baseline/data_derived/analytic_cohorts_v2.rds'))
put(d %>% group_by(spid) %>% summarise(entries=n(),groups=n_distinct(digi3),death=any(evt_death==1),residential=any(evt_nh==1),composite=any(evt_comp==1),.groups='drop') %>% count(entries,groups),'R2_1_repeat_and_exposure_change')
put(first %>% group_by(digi3) %>% summarise(N=n(),death=sum(evt_death),residential=sum(evt_nh),composite=sum(evt_comp)),'R2_1_unique_events')
put(d %>% group_by(cohort) %>% summarise(person_baselines=n(),unweighted_offline=mean(digi3=='完全断连')*100,weighted_offline=weighted.mean(digi3=='完全断连',w0)*100),'R2_12_prevalence_by_entry')
f<-old %>% group_by(cohort) %>% summarise(eligible=n(),missing_exposure=sum(is.na(digi3)),known_exposure=sum(!is.na(digi3)),no_followup=sum(!is.na(digi3)&fu_rounds==0))
f$eligible_after_followup<-f$known_exposure-f$no_followup
for(b in 11:13){p<-filter(dat$pp_comp,cohort==b);key<-paste(p$spid,p$cohort);cv<-covs$M4;cc<-complete.cases(p[cv]);f$M4_covariate_excluded[f$cohort==b]<-length(unique(key))-length(unique(key[cc]));f$M4_with_entry_weight[f$cohort==b]<-length(unique(key[cc]));f$M4_positive_followup_weight[f$cohort==b]<-length(unique(key[cc&is.finite(p$w_follow)&p$w_follow>0]))}
put(f,'R2_11_flow_by_entry')
put(data.frame(variable=c('Exposure','Age/sex/race/education/Medicaid/marital/SRH','Income','PHQ2/memory/going outside','ADL'),original=c('641 excluded; not imputed','complete-case in primary models','unknown category','unknown category for unavailable response','original count coding retained'),supplement=c('unchanged exclusion','first-entry MICE for item nonresponse','20 official draws plus residual missing imputation','structural skips retained; ordinary unknown imputed in MI','original definition retained under scope instruction')),'R2_4_missing_handling_rules')
put(data.frame(round=11:14,field_start=c('2021-06','2022-03','2023-05','2024-05'),field_end=c('2021-11','2022-12','2024-04','2025-03'),source='NHATS R1-R14 User Guide February 2026 p16'),'R2_5_field_periods')
cat('DESCRIPTIVES_DONE\n')
