if(.Platform$OS.type=="windows")try(Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8"),silent=TRUE)
ROOT<-commandArgs(TRUE)[1]
RAW<-Sys.getenv("NHATS_RAW_DIR")
## ============================================================
## NHATS 第四篇 V2 — 04_main：主分析（审稿修订版）
## #4 stacked 主模型(队列固定效应+spid聚类) —— 主结果
##    per-cohort + meta 合并 —— 一致性检验
## #9 M4 老年脆弱模型(M3 + ADL/抑郁/认知/收入/外出)
## + Table 1(含 M4 维度) + E-value + 敏感性
## ============================================================
suppressPackageStartupMessages({
  library(haven);library(dplyr);library(tidyr);library(survey)})
options(survey.lonely.psu="adjust")
PROJ <- file.path(ROOT,"baseline")
TAB  <- file.path(PROJ,"tables"); dir.create(TAB,showWarnings=FALSE)

coh <- readRDS(file.path(PROJ,"data_derived","analytic_cohorts.rds"))
pp  <- readRDS(file.path(PROJ,"data_derived","analytic_personperiod.rds"))

## ============================================================
## 0. 派生暴露与 M4 协变量
## ============================================================
## 收入三分位切点（基于 person-baseline 正值收入）
inc_cut <- quantile(coh$totinc[coh$totinc>0], c(1/3,2/3), na.rm=TRUE)
cat("收入三分位切点：", round(inc_cut), "\n")

derive <- function(d){
  d %>% mutate(
    offline = factor(ifelse(digi3=="完全断连","断连","非断连"),
                     levels=c("非断连","断连")),
    ## M4：ADL 求助计数（吃/洗澡/如厕/穿衣，值1=有人协助）
    adl_n   = (sc_eat==1)+(sc_bath==1)+(sc_toil==1)+(sc_dres==1),
    adl_cat = factor(pmin(adl_n,2), levels=0:2,
                     labels=c("无ADL求助","1项ADL求助","≥2项ADL求助")),
    ## M4：抑郁 PHQ-2（depresan1/2，1-4 → 0-3；和≥3 为阳性）
    phq2    = ifelse(dep1 %in% 1:4 & dep2 %in% 1:4, (dep1-1)+(dep2-1), NA_real_),
    dep_f   = factor(ifelse(is.na(phq2),"抑郁筛查未知",
                     ifelse(phq2>=3,"抑郁筛查阳性","抑郁筛查阴性")),
                     levels=c("抑郁筛查阴性","抑郁筛查阳性","抑郁筛查未知")),
    ## M4：自评记忆
    mem_f   = factor(ifelse(ratemem %in% 1:2,"记忆好",
                     ifelse(ratemem==3,"记忆一般",
                     ifelse(ratemem %in% 4:5,"记忆差","记忆未知"))),
                     levels=c("记忆好","记忆一般","记忆差","记忆未知")),
    ## M4：外出频率（居家受限代理）
    out3    = factor(ifelse(outoft %in% 1:2,"经常外出",
                     ifelse(outoft==3,"有时外出",
                     ifelse(outoft %in% 4:5,"罕见/不外出","外出未知"))),
                     levels=c("经常外出","有时外出","罕见/不外出","外出未知")),
    ## M4：家庭收入三分位
    inc3    = factor(ifelse(is.na(totinc)|totinc<0,"收入未知",
                     ifelse(totinc<=inc_cut[1],"收入低",
                     ifelse(totinc<=inc_cut[2],"收入中","收入高"))),
                     levels=c("收入中","收入低","收入高","收入未知")),
    tf      = factor(t %||% NA))
}
`%||%` <- function(a,b) if(is.null(a)) b else a
coh <- derive(coh %>% mutate(t=NA))
pp  <- derive(pp)
pp_m  <- pp  %>% filter(!is.na(digi3))
coh_m <- coh %>% filter(!is.na(digi3))
cat("分析集：person-period", nrow(pp_m), "行；person-baseline", nrow(coh_m), "\n")

## ============================================================
## 1. 协变量分块（新增 M4）
## ============================================================
covs <- list(
  M0=character(0),
  M1=c("agec","sex","race","educ3"),
  M2=c("agec","sex","race","educ3","medicaid_f","livealone","married"),
  M3=c("agec","sex","race","educ3","medicaid_f","livealone","married",
       "como","srh_f","dementia_f","estcvd_f"),
  M4=c("agec","sex","race","educ3","medicaid_f","livealone","married",
       "como","srh_f","dementia_f","estcvd_f",
       "adl_cat","dep_f","mem_f","out3","inc3"))

fmtHR <- function(beta,se){
  sprintf("%.2f (%.2f-%.2f)",exp(beta),exp(beta-1.96*se),exp(beta+1.96*se))}

## ============================================================
## 2. STACKED 主模型（队列固定效应 + spid 聚类稳健方差）
## ============================================================
fit_stacked <- function(outcome, expo, model){
  cv <- covs[[model]]
  terms <- c(expo, "cohort_f", "tf", cv)
  fml <- as.formula(paste(outcome,"~",paste(terms,collapse="+")))
  des <- svydesign(ids=~spid, weights=~w0, data=pp_m)
  m <- tryCatch(svyglm(fml,design=des,family=quasibinomial(link="cloglog")),
                error=function(e){cat("stacked 失败:",model,expo,conditionMessage(e),"\n");NULL})
  if(is.null(m)) return(NULL)
  s <- summary(m)$coefficients
  keep <- grep(paste0("^",expo), rownames(s))
  data.frame(term=rownames(s)[keep], beta=s[keep,1], se=s[keep,2],
             n_pp=nrow(model.frame(m)), method="stacked", model=model,
             stringsAsFactors=FALSE)
}

## ============================================================
## 3. per-cohort 模型（NHATS 复杂抽样设计，一致性检验）
## ============================================================
fit_cohort <- function(dat, outcome, expo, model, drop_t=FALSE){
  cv <- covs[[model]]
  terms <- c(expo, if(!drop_t) "tf", cv)
  fml <- as.formula(paste(outcome,"~",paste(terms,collapse="+")))
  des <- svydesign(ids=~vunit,strata=~vstrat,weights=~w0,data=dat,nest=TRUE)
  m <- tryCatch(svyglm(fml,design=des,family=quasibinomial(link="cloglog")),
                error=function(e) NULL)
  if(is.null(m)) return(NULL)
  s <- summary(m)$coefficients
  keep <- grep(paste0("^",expo), rownames(s))
  if(length(keep)==0) return(NULL)
  data.frame(term=rownames(s)[keep], beta=s[keep,1], se=s[keep,2],
             n_pp=nrow(model.frame(m)), stringsAsFactors=FALSE)
}
run_cohorts <- function(outcome, expo){
  res <- list()
  for(b in c(11,12,13)){
    d <- pp_m %>% filter(cohort==b)
    dropt <- (length(unique(d$t))<2)
    for(mm in names(covs)){
      r <- fit_cohort(d,outcome,expo,mm,drop_t=dropt)
      if(!is.null(r)){ r$cohort<-as.character(b); r$model<-mm; res[[paste(b,mm)]]<-r }
    }
  }
  bind_rows(res)
}
meta_pool <- function(df){
  df %>% group_by(term,model) %>%
    summarise(k=n(), beta_p=sum(beta/se^2)/sum(1/se^2),
              se_p=sqrt(1/sum(1/se^2)),
              Q=sum((beta-sum(beta/se^2)/sum(1/se^2))^2/se^2),
              .groups="drop") %>%
    mutate(I2=pmax(0,(Q-(k-1))/Q)*100, cohort="meta", beta=beta_p, se=se_p) %>%
    select(term,model,cohort,beta,se,I2,k)
}

## ============================================================
## 4. 汇总主表（stacked + per-cohort + meta），两结局 × 三对比
## ============================================================
build_all <- function(outcome,label){
  ## stacked
  st <- bind_rows(lapply(names(covs), function(mm)
          bind_rows(fit_stacked(outcome,"digi3",mm),
                    fit_stacked(outcome,"offline",mm))))
  st <- st %>% transmute(结局=label,term,model,cohort="stacked",
                         HR=fmtHR(beta,se),I2=NA_real_)
  ## per-cohort + meta
  raw <- bind_rows(run_cohorts(outcome,"digi3"), run_cohorts(outcome,"offline"))
  mp  <- meta_pool(raw)
  pc  <- bind_rows(
           raw %>% transmute(结局=label,term,model,cohort,
                             HR=fmtHR(beta,se),I2=NA_real_),
           mp  %>% transmute(结局=label,term,model,cohort,
                             HR=fmtHR(beta,se),I2))
  bind_rows(st,pc)
}
T_comp  <- build_all("evt_comp","复合硬终点(死亡或养老院)")
T_death <- build_all("evt_death","全因死亡")
T_main  <- bind_rows(T_comp,T_death) %>%
  mutate(term=recode(term, digi3远程使用者="远程使用者 vs 联网非使用者",
                     digi3完全断连="完全断连 vs 联网非使用者",
                     offline断连="完全断连 vs 非断连")) %>%
  arrange(结局,term,model,match(cohort,c("stacked","11","12","13","meta")))
write.csv(T_main, file.path(TAB,"T3_main_HR_v2.csv"),
          row.names=FALSE, fileEncoding="UTF-8")
cat("\n主表 T3 完成，", nrow(T_main), "行\n")
cat("\n=== STACKED 主结果（完全断连 vs 非断连）===\n")
print(T_main %>% filter(term=="完全断连 vs 非断连", cohort=="stacked"))
cat("\n=== STACKED 主结果（远程 vs 联网非使用）===\n")
print(T_main %>% filter(term=="远程使用者 vs 联网非使用者", cohort=="stacked"))

## ============================================================
## 5. E-value（stacked M3 与 M4，完全断连 vs 非断连）
## ============================================================
evalue <- function(est){ e<-ifelse(est<1,1/est,est); e+sqrt(e*(e-1)) }
ev_rows <- list()
for(oc in c("evt_comp","evt_death")){
  for(mm in c("M3","M4")){
    r <- fit_stacked(oc,"offline",mm)
    hr <- exp(r$beta); lo <- exp(r$beta-1.96*r$se)
    ev_rows[[paste(oc,mm)]] <- data.frame(
      结局=ifelse(oc=="evt_comp","复合硬终点","全因死亡"), 模型=mm,
      HR=sprintf("%.2f",hr),
      `HR_E.value`=sprintf("%.2f",evalue(hr)),
      `CI下界_E.value`=sprintf("%.2f",evalue(lo)))
  }
}
T_ev <- bind_rows(ev_rows)
write.csv(T_ev, file.path(TAB,"T5_evalue_v2.csv"), row.names=FALSE, fileEncoding="UTF-8")
cat("\nE-value 完成\n"); print(T_ev)

## ============================================================
## 6. 敏感性分析（stacked：排除首区间 / 限自答者）
## ============================================================
sens <- list()
## (a) 排除入组后第一区间事件
des_a <- svydesign(ids=~spid,weights=~w0,data=pp_m %>% filter(t>=2))
m_a <- svyglm(as.formula(paste("evt_comp ~ offline+cohort_f+tf+",
        paste(covs$M3,collapse="+"))),design=des_a,family=quasibinomial("cloglog"))
sa <- summary(m_a)$coefficients["offline断连",]
sens[["排除首区间事件(stacked)"]] <- data.frame(分析="排除首区间事件(stacked)",
  HR=fmtHR(sa[1],sa[2]), n_pp=nrow(model.frame(m_a)))
## (b) 限自答者
des_b <- svydesign(ids=~spid,weights=~w0,data=pp_m %>% filter(proxy_f=="自答"))
m_b <- svyglm(as.formula(paste("evt_comp ~ offline+cohort_f+tf+",
        paste(covs$M3,collapse="+"))),design=des_b,family=quasibinomial("cloglog"))
sb <- summary(m_b)$coefficients["offline断连",]
sens[["限自答者(stacked)"]] <- data.frame(分析="限自答者(stacked)",
  HR=fmtHR(sb[1],sb[2]), n_pp=nrow(model.frame(m_b)))
## (c) per-cohort 敏感性（保留，与原文一致）
for(b in c(11,12)){
  d <- pp_m %>% filter(cohort==b, t>=2); if(nrow(d)<50) next
  r <- fit_cohort(d,"evt_comp","offline","M3",drop_t=(length(unique(d$t))<2))
  if(!is.null(r)) sens[[paste0("排除首区间_R",b)]] <-
    data.frame(分析=paste0("排除首区间事件_R",b),HR=fmtHR(r$beta,r$se),n_pp=r$n_pp)
}
for(b in c(11,12,13)){
  d <- pp_m %>% filter(cohort==b, proxy_f=="自答"); if(nrow(d)<50) next
  r <- fit_cohort(d,"evt_comp","offline","M3",drop_t=(length(unique(d$t))<2))
  if(!is.null(r)) sens[[paste0("限自答_R",b)]] <-
    data.frame(分析=paste0("限自答者_R",b),HR=fmtHR(r$beta,r$se),n_pp=r$n_pp)
}
T_sens <- bind_rows(sens)
write.csv(T_sens, file.path(TAB,"T6_sensitivity_v2.csv"), row.names=FALSE, fileEncoding="UTF-8")
cat("\n敏感性分析完成\n"); print(T_sens)

## ============================================================
## 7. 保存
## ============================================================
saveRDS(list(T_main=T_main,T_ev=T_ev,T_sens=T_sens,inc_cut=inc_cut),
        file.path(PROJ,"data_derived","model_results_v2.rds"))
saveRDS(coh, file.path(PROJ,"data_derived","analytic_cohorts_v2.rds"))
saveRDS(pp,  file.path(PROJ,"data_derived","analytic_personperiod_v2.rds"))
cat("\n04_main 完成。\n")
