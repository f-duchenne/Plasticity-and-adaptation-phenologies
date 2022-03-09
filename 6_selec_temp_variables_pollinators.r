.libPaths(c("/home/duchenne/R/x86_64-pc-linux-gnu-library/3.4/",.libPaths()))
library(data.table)
library(dplyr)
library(glmmTMB)
#library(doMPI)
#library(igraph)
library(doParallel)
library(foreach)
library(parallel)
# cl<-startMPIcluster()
# registerDoMPI(cl)
#cl<-makeCluster(48)
#registerDoParallel(cl)
setwd(dir="D:/cluster/plast")
donnf2=fread("dataf_avec_elev.txt",sep="\t")
donnf2$latitude=round(donnf2$latitude,5)
donnf2$longitude=round(donnf2$longitude,5)
tabtemp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")
tabtemp=tabtemp[,-which(names(tabtemp)=="temp_360_85"),with=F]
donnf2=merge(donnf2,tabtemp,by=c("latitude","longitude","Annee"))
liste=as.data.frame(unique(donnf2[,c("Speciesgen","Species","FAMILLE","ORDRE")]))
liste[,names(tabtemp)[4:ncol(tabtemp)]]=NA
liste$rsq=NA
liste$nb_pre1990=NA
liste$nb_post1990=NA
liste$nb=NA
liste$inter=NA
liste$sing=NA
donnf2=subset(donnf2,!is.na(elev))
subs=matrix(TRUE,ncol(tabtemp)-3,ncol(tabtemp)-3)
for (i in 1:ncol(subs)){
if(i<ncol(subs)){subs[i+1,i]=FALSE
subs[i,i+1]=FALSE}
}
colnames(subs)=c("temp_0_90","temp_30_120","temp_60_150","temp_90_180","temp_120_210","temp_150_240","temp_180_270",
"temp_210_300","temp_240_330","temp_270_360","temp_300_25","temp_330_55") 
rownames(subs)=c("temp_0_90","temp_30_120","temp_60_150","temp_90_180","temp_120_210","temp_150_240","temp_180_270",
"temp_210_300","temp_240_330","temp_270_360","temp_300_25","temp_330_55") 
subs[nrow(subs),1]=FALSE
subs[1,nrow(subs)]=FALSE
subs=as.data.frame(subs)
subs$vari=rownames(subs)
subs=melt(subs,id.var="vari")
subs=subset(subs,value=="FALSE")

bid1=as.data.frame(t(combn(c("temp_0_90","temp_30_120","temp_60_150","temp_90_180","temp_120_210","temp_150_240","temp_180_270",
"temp_210_300","temp_240_330","temp_270_360","temp_300_25","temp_330_55"),m=1)))
names(bid1)=c("var1")
bid1$var2=NA
bid1$var3=NA
bid2=as.data.frame(t(combn(c("temp_0_90","temp_30_120","temp_60_150","temp_90_180","temp_120_210","temp_150_240","temp_180_270",
"temp_210_300","temp_240_330","temp_270_360","temp_300_25","temp_330_55"),m=2)))
names(bid2)=c("var1","var2")
bid2$var3=NA
bid3=as.data.frame(t(combn(c("temp_0_90","temp_30_120","temp_60_150","temp_90_180","temp_120_210","temp_150_240","temp_180_270",
"temp_210_300","temp_240_330","temp_270_360","temp_300_25","temp_330_55"),m=3)))
names(bid3)=c("var1","var2","var3")
bidf=rbind(bid1,bid2,bid3)
bidf$autor=NA
for(i in 1:nrow(bidf)){
bidf$autor[i]=max(apply(subs,1,function(x)length(which(t(bidf[i,]) %in% t(x))))) 
}
bidf=subset(bidf,autor<2)

rm(i)
resultat=foreach(i=1:nrow(liste),.combine=rbind)%dopar%{
library(data.table)
library(dplyr)
library(lme4)
library(spaMM)
bidon=subset(donnf2,Speciesgen==liste$Speciesgen[i] & !is.na(temp_0_90))
bidon$maille=as.factor(bidon$maille)
if(length(unique(bidon$maille))<nrow(bidon) & nrow(bidon)>50){
a=0
bidon$elev=(bidon$elev-mean(bidon$elev,na.rm=T))/1000
dredgi=bidf
dredgi$aicc=NA
if(length(unique(bidon$maille))>1){
for(j in 1:nrow(dredgi)){
varia=t(dredgi)[1:3,j][which(!is.na(t(dredgi)[1:3,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1+",paste(varia,collapse="+"),"|maille)+(1|Annee)"))
model=fitme(form,data=bidon)
dredgi$aicc[j]=extractAIC(model)[2]
}
dredgi=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:3,1][which(!is.na(t(dredgi)[1:3,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1+",paste(varia,collapse="+"),"|maille)+(1|Annee)"))
model=fitme(form,data=bidon,REML=T)
}else{
for(j in 1:nrow(dredgi)){
varia=t(dredgi)[1:3,j][which(!is.na(t(dredgi)[1:3,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1|Annee)"))
model=fitme(form,data=bidon)
dredgi$aicc[j]=extractAIC(model)[2]
}
dredgie=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:3,1][which(!is.na(t(dredgi)[1:3,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1|Annee)"))
model=fitme(form,data=bidon)
	}
bidon$resid=residuals(model)
liste$rsq[i]=NA
liste[i,names(model$fixef[-1])]=model$fixef[-1]
liste$nb_pre1990[i]=nrow(subset(bidon,Annee<1990))
liste$nb_post1990[i]=nrow(subset(bidon,Annee>=1990))
liste$inter[i]=model$fixef[1,1]
if(length(unique(bidon$maille))>1){liste$sing[i]=isSingular(model)}else{liste$sing[i]=NA}
}else{
liste$rsq[i]=NA
liste[i,names(tabtemp)[4:ncol(tabtemp)]]=NA
liste$nb[i]=NA
liste$inter[i]=NA
liste$sing[i]=NA
}
# resultat=rbind(resultat,liste[i,])
return(liste[i,])
}

fwrite(resultat,"selec_temp_var_beta.txt",sep="\t")

closeCluster(cl)
mpi.quit()



