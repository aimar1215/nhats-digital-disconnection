Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8');library(dplyr)
ROOT<-commandArgs(TRUE)[1]
a<-read.csv(file.path(ROOT,'tables/MI_individual_estimates.csv'));b<-read.csv(file.path(ROOT,'tables/MI_first_entry_pooled.csv'));dg<-read.csv(file.path(ROOT,'tables/MI_diagnostics.csv'))
stopifnot(nrow(a)==160,nrow(b)==8,nrow(dg)==20,all(dg$missing_remaining==0),all(b$m==20),all(is.finite(b$lo)),all(is.finite(b$hi)),all(b$B>=0))
for(i in seq_len(nrow(b))){z<-a[a$endpoint==b$endpoint[i]&a$model==b$model[i]&a$weight==b$weight[i]&a$term==b$term[i],];stopifnot(abs(var(z$beta)-b$B[i])<1e-12,abs(mean(z$se^2)+1.05*var(z$beta)-b$Tvar[i])<1e-12)}
cat('PASS: 20 completed imputations, 160 estimates, 8 pooled rows, independent Rubin variance check.\n')
