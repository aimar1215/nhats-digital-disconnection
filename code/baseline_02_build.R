if(.Platform$OS.type=="windows")try(Sys.setlocale("LC_CTYPE","Chinese (Simplified)_China.utf8"),silent=TRUE)
ROOT<-commandArgs(TRUE)[1]
RAW<-Sys.getenv("NHATS_RAW_DIR")
## ============================================================
## NHATS 第四篇 V2 — 02_build：错开入组分析数据集构建（扩展版）
## 在基线复现版基础上，新增 M4 老年脆弱协变量原始题项 + 失访指标
## 队列构建逻辑与原 02_build.R 完全一致（12,139 person-baseline 不变）
## 新增列：ADL 求助 / 移动协助 / 抑郁 PHQ-2 / 自评记忆 / 收入 / 外出频率 / ltfu
## ============================================================
suppressPackageStartupMessages({library(haven); library(dplyr); library(tidyr)})

DATA <- RAW
TRK <- file.path(RAW,"NHATS_Round_14_Tracker_File.dta")
PROJ <- file.path(ROOT,"baseline")
spf  <- c("10"="NHATS_Round_10B_SP_File.dta","11"="NHATS_Round_11_SP_File_V2.dta",
          "12"="NHATS_Round_12_SP_File.dta","13"="NHATS_Round_13_SP_File.dta",
          "14"="NHATS_Round_14_SP_File.dta")

neg2na <- function(x){ x <- as.numeric(x); x[x < 0] <- NA; x }

## ---- 读取四轮 SP 文件（R11–R14 分析用）----
read_round <- function(rd){
  d <- read_dta(file.path(DATA, spf[rd]))
  gv <- function(stem) if (stem %in% names(d)) d[[stem]] else rep(NA_real_, nrow(d))
  dis <- sapply(1:10, function(k) gv(paste0("hc",rd,"disescn",k)))
  como <- rowSums(matrix(dis %in% c(1,7), ncol=10), na.rm=TRUE)
  tibble(
    spid = as.numeric(d[["spid"]]), round = as.integer(rd),
    hd  = gv(paste0("hc",rd,"disescn2")), htn = gv(paste0("hc",rd,"disescn3")),
    dm  = gv(paste0("hc",rd,"disescn6")), mi  = gv(paste0("hc",rd,"disescn1")),
    strk= gv(paste0("hc",rd,"disescn8")), dementia = gv(paste0("hc",rd,"disescn9")),
    como = como, srh = neg2na(gv(paste0("hc",rd,"health"))),
    online = gv(paste0("te",rd,"online")), teleh = gv(paste0("te",rd,"intrntmd4")),
    mctele = gv(paste0("mc",rd,"telehealth")), regdoc = neg2na(gv(paste0("mc",rd,"havregdoc"))),
    dresid = gv(paste0("r",rd,"dresid")),
    age = neg2na(gv(paste0("r",rd,"d2intvrage"))), sex = neg2na(gv(paste0("r",rd,"dgender"))),
    marital = neg2na(gv(paste0("hh",rd,"dmarstat"))), hhnum = neg2na(gv(paste0("hh",rd,"dhshldnum"))),
    medicaid = neg2na(gv(paste0("ip",rd,"cmedicaid"))), proxyid = gv(paste0("is",rd,"dproxyid")),
    w0 = gv(paste0("w",rd,"anfinwgt0")), vstrat = gv(paste0("w",rd,"varstrat")),
    vunit = gv(paste0("w",rd,"varunit")),
    ## ---- M4 老年脆弱协变量原始题项（入组轮取值）----
    sc_eat  = gv(paste0("sc",rd,"eathlp")),  sc_bath = gv(paste0("sc",rd,"bathhlp")),
    sc_toil = gv(paste0("sc",rd,"toilhlp")), sc_dres = gv(paste0("sc",rd,"dreshlp")),
    mobhelp = gv(paste0("dm",rd,"dmohlp")),
    dep1    = gv(paste0("hc",rd,"depresan1")), dep2 = gv(paste0("hc",rd,"depresan2")),
    ratemem = gv(paste0("cg",rd,"ratememry")),
    totinc  = gv(paste0("ia",rd,"totinc")),
    outoft  = gv(paste0("mo",rd,"outoft"))
  )
}
sp <- bind_rows(lapply(c("11","12","13","14"), read_round))
cat("各轮 SP 行数：\n"); print(table(sp$round))

## ---- 时不变人口学(种族/教育)查找表：从所有含该变量的轮次取首个非缺失 ----
read_demog <- function(rd){
  d <- read_dta(file.path(DATA, spf[rd]))
  gv <- function(stem) if (stem %in% names(d)) neg2na(d[[stem]]) else rep(NA_real_, nrow(d))
  tibble(spid=as.numeric(d[["spid"]]),
         race=gv(paste0("rl",rd,"dracehisp")), educ=gv(paste0("el",rd,"dhigstschl")),
         gender=gv(paste0("r",rd,"dgender")))
}
f1 <- function(x){ x <- na.omit(x); if(length(x)==0) NA_real_ else x[1] }
demog <- bind_rows(lapply(c("10","12","13"), read_demog)) %>%
  group_by(spid) %>%
  summarise(race=f1(race), educ=f1(educ), gender=f1(gender), .groups="drop")
cat("时不变人口学查找表 N =", nrow(demog), "；race 覆盖", sum(!is.na(demog$race)),
    "；educ 覆盖", sum(!is.na(demog$educ)), "；gender 覆盖", sum(!is.na(demog$gender)), "\n")

## ---- 关键变量缺失体检 ----
chk <- sp %>% group_by(round) %>%
  summarise(across(c(age,sex,marital,hhnum,medicaid,regdoc,srh,dementia,vstrat,w0),
                   ~round(mean(is.na(.))*100,1), .names="{.col}"))
cat("\n关键变量缺失率(%)：\n"); print(as.data.frame(chk))

## ---- 死亡：Tracker ----
tk <- read_dta(TRK, col_select=c("spid", paste0("r",11:14,"status")))
tk$spid <- as.numeric(tk$spid)
deadcode <- c(62,81,86)
death_round <- apply(tk[,paste0("r",11:14,"status")], 1, function(s){
  hit <- which(s %in% deadcode); if(length(hit)==0) NA_integer_ else (10+min(hit))})
trk <- tibble(spid=tk$spid, death_round=death_round)

## ---- 暴露三梯度 ----
sp <- sp %>% mutate(
  digi3 = case_when(
    teleh==1 ~ "1_telehealth_user",
    online==1 & teleh==2 ~ "2_connected_nonuser",
    online==2 ~ "3_offline", TRUE ~ NA_character_))

## ---- 宽表：各轮 dresid / 观测 ----
res_wide <- sp %>% select(spid,round,dresid) %>%
  pivot_wider(names_from=round, values_from=dresid, names_prefix="res")
obs_wide <- sp %>% mutate(obs=1L) %>% select(spid,round,obs) %>%
  pivot_wider(names_from=round, values_from=obs, names_prefix="obs")

## ---- 构建错开入组队列 ----
build_cohort <- function(b){
  base <- sp %>% filter(round==b)
  cv <- (base$hd %in% c(1,7)) | (base$htn %in% c(1,7)) | (base$dm %in% c(1,7))
  hist <- sp %>% filter(round<=b) %>% group_by(spid) %>%
    summarise(estcvd = any(mi==1,na.rm=TRUE) | any(strk==1,na.rm=TRUE) |
                       any(hd %in% c(1,7),na.rm=TRUE), .groups="drop")
  d <- base %>% filter(cv) %>%
    left_join(hist,by="spid") %>% left_join(trk,by="spid") %>%
    left_join(res_wide,by="spid") %>% left_join(obs_wide,by="spid") %>%
    mutate(cohort=b) %>% filter(dresid==1 | is.na(dresid))
  rescols <- paste0("res",(b+1):14); obscols <- paste0("obs",(b+1):14)
  resm <- as.matrix(d[,rescols,drop=FALSE]); obsm <- as.matrix(d[,obscols,drop=FALSE])
  d$nh_round <- apply(resm,1,function(r){
    hit <- which(r %in% c(2,3,4,5)); if(length(hit)==0) NA_integer_ else (b+min(hit))})
  d$last_obs <- apply(obsm,1,function(o){
    hit <- which(o==1); if(length(hit)==0) b else (b+max(hit))})
  ev <- pmin(d$death_round, d$nh_round, na.rm=TRUE)
  ev[is.infinite(ev)] <- NA; ev[!is.na(ev) & (ev<=b | ev>14)] <- NA
  d$event_round <- ev
  d$end_round <- ifelse(!is.na(d$event_round), d$event_round, pmin(pmax(d$last_obs,b),14))
  d$evt_comp  <- as.integer(!is.na(d$event_round))
  d$evt_death <- as.integer(!is.na(d$death_round) & d$death_round<=14 & d$death_round>b)
  d$evt_nh    <- as.integer(!is.na(d$nh_round))
  d$fu_rounds <- pmax(d$end_round - b, 0L)
  d
}
coh <- bind_rows(lapply(c(11,12,13), build_cohort))

## ---- 合并时不变人口学 + 协变量重编码 ----
coh <- coh %>% left_join(demog, by="spid") %>% mutate(
  digi3 = factor(digi3, levels=c("2_connected_nonuser","1_telehealth_user","3_offline"),
                 labels=c("联网非使用者","远程使用者","完全断连")),
  sex   = factor(coalesce(sex, gender), levels=c(1,2), labels=c("男","女")),
  race  = factor(case_when(race==1~"非西语裔白人", race==2~"非西语裔黑人",
                           race==4~"西语裔", is.na(race)~NA_character_, TRUE~"其他/多种"),
                 levels=c("非西语裔白人","非西语裔黑人","西语裔","其他/多种")),
  educ3 = factor(case_when(educ %in% 1:2 ~"高中以下", educ==3~"高中/同等",
                           educ>=4 ~"高中以上", TRUE~NA_character_),
                 levels=c("高中以上","高中/同等","高中以下")),
  livealone = factor(ifelse(hhnum==1,"独居","非独居"), levels=c("非独居","独居")),
  married   = factor(ifelse(marital==1,"在婚","非在婚"), levels=c("在婚","非在婚")),
  medicaid_f= factor(ifelse(medicaid==1,"是","否"), levels=c("否","是")),
  dementia_f= factor(ifelse(dementia %in% c(1,7),"是","否"), levels=c("否","是")),
  srh_f     = factor(srh, levels=1:5, labels=c("极好","很好","好","一般","差")),
  srh_bad   = factor(ifelse(srh>=4,"一般/差","极好~好"), levels=c("极好~好","一般/差")),
  estcvd_f  = factor(ifelse(estcvd,"是","否"), levels=c("否","是")),
  proxy_f   = factor(ifelse(!is.na(proxyid) & proxyid>0,"代答","自答"), levels=c("自答","代答")),
  agec      = factor(age),
  cohort_f  = factor(cohort, levels=c(11,12,13), labels=c("R11入组","R12入组","R13入组")),
  ## ---- 失访指标：无复合事件且末次观测早于 R14 ----
  ltfu      = as.integer(evt_comp==0 & end_round<14)
)

## ---- person-period 长表（离散时间生存）----
make_pp <- function(d){
  d <- d %>% filter(fu_rounds>=1)
  bind_rows(lapply(seq_len(nrow(d)), function(i){
    r <- d[i,]; T <- r$fu_rounds
    tibble(spid=r$spid, cohort=r$cohort, t=1:T,
           evt_comp =ifelse(1:T==T & r$evt_comp==1,1L,0L),
           evt_death=ifelse(1:T==T & r$evt_death==1,1L,0L))
  }))
}
pp <- make_pp(coh) %>% left_join(
  coh %>% select(spid,cohort,digi3,sex,race,educ3,agec,livealone,married,medicaid_f,
                 dementia_f,srh,srh_f,srh_bad,como,estcvd_f,proxy_f,regdoc,mctele,
                 cohort_f,w0,vstrat,vunit,
                 sc_eat,sc_bath,sc_toil,sc_dres,mobhelp,dep1,dep2,ratemem,totinc,outoft),
  by=c("spid","cohort"))

## ---- 保存 ----
dir.create(file.path(PROJ,"data_derived"), showWarnings=FALSE)
saveRDS(coh, file.path(PROJ,"data_derived","analytic_cohorts.rds"))
saveRDS(pp,  file.path(PROJ,"data_derived","analytic_personperiod.rds"))

## ---- 核验摘要 ----
cat("\n========== 队列摘要 ==========\n")
print(as.data.frame(coh %>% group_by(cohort_f, digi3) %>%
  summarise(N=n(), 死亡=sum(evt_death), 养老院=sum(evt_nh),
            复合事件=sum(evt_comp), 中位随访轮=median(fu_rounds), .groups="drop")))
cat("\n建模协变量缺失率(%)，三队列合并：\n")
print(round(sapply(coh[,c("race","educ3","agec","sex","srh_f","medicaid_f","livealone")],
                   function(x) mean(is.na(x))*100),1))
cat("\nperson-period 行数：", nrow(pp), "；暴露缺失：", sum(is.na(coh$digi3)),"/",nrow(coh),"\n")

## ---- M4 协变量原始题项取值核查（用于确定派生编码）----
cat("\n========== M4 协变量原始取值分布 ==========\n")
for(v in c("sc_eat","sc_bath","sc_toil","sc_dres","mobhelp","dep1","dep2",
           "ratemem","outoft")){
  cat("--", v, "--\n"); print(table(coh[[v]], useNA="ifany"))
}
cat("-- totinc（家庭总收入）摘要 --\n")
print(summary(coh$totinc))
cat("totinc 为负/缺失比例：", round(mean(coh$totinc<0 | is.na(coh$totinc),na.rm=TRUE)*100,1),"%\n")
cat("失访(ltfu) × 暴露组：\n")
print(table(coh$digi3, coh$ltfu, useNA="ifany"))
cat("02_build 完成。\n")
