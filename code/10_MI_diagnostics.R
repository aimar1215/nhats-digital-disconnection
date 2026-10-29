Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
suppressPackageStartupMessages({library(mice);library(dplyr)})
ROOT<-commandArgs(TRUE)[1]
ims<-lapply(1:20,function(k)readRDS(file.path(ROOT,'data',paste0('mice_first_entry_',k,'.rds'))))
logs<-bind_rows(lapply(1:20,function(k){x<-ims[[k]]$loggedEvents;if(is.null(x))return(NULL);x$imputation<-k;x}))
write.csv(logs,file.path(ROOT,'tables/MI_logged_events.csv'),row.names=FALSE,fileEncoding='UTF-8')
methods<-bind_rows(lapply(1:20,function(k)data.frame(imputation=k,variable=names(ims[[k]]$method),method=ims[[k]]$method)))
write.csv(methods,file.path(ROOT,'tables/MI_methods.csv'),row.names=FALSE,fileEncoding='UTF-8')
vars<-names(ims[[1]]$method)[ims[[1]]$method!='']
png(file.path(ROOT,'figures/MI_chain_mean_traces.png'),width=2700,height=2500,res=300,type='cairo',family='Times New Roman')
par(mfrow=c(3,3),mar=c(4,4,2,1));for(v in vars){mat<-sapply(ims,function(im)im$chainMean[v,,1]);matplot(1:10,mat,type='l',lty=1,col=adjustcolor('black',.25),xlab='Iteration',ylab='Mean imputed code',main=v)};dev.off()
cat('COMPLETE_DATA_MISSING',sapply(ims,function(x)sum(is.na(complete(x)))),'\n')
cat('Methods',unique(methods$method),'\n');cat('Logged events',nrow(logs),'primarily automatic collinearity exclusions\n')
