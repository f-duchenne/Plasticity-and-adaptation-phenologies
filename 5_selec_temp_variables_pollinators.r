library(data.table)
library(dplyr)
#library(doMPI)
#library(igraph)
library(doParallel)
library(foreach)
library(parallel)


setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
donnf2=fread("dataf_avec_elev.txt",sep="\t")
donnf2$latitude=round(donnf2$latitude,5)
donnf2$longitude=round(donnf2$longitude,5)
tabtemp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")
tabtemp$latitude=round(tabtemp$latitude,5)
tabtemp$longitude=round(tabtemp$longitude,5)
donnf2=merge(donnf2,tabtemp,by=c("latitude","longitude","Annee"))
dim(donnf2)
liste=as.data.frame(unique(donnf2[,c("Speciesgen","Species","FAMILLE","ORDRE")]))
liste[,names(tabtemp)[4:ncol(tabtemp)]]=NA
liste$rsq=NA
liste$nb_pre1990=NA
liste$nb_post1990=NA
liste$nb=NA
liste$inter=NA
liste$sing=NA
donnf2=subset(donnf2,!is.na(elev)  & !is.na(temp_0_90))

#defining temperature indice which are overlapping too much
couples=data.frame(deb=rep(seq(0,365,15),3),fin= c(seq(0,365,15)+30,seq(0,365,15)+60,seq(0,365,15)+90))
couples=subset(couples,deb<360)
subs=matrix(TRUE,nrow(couples),nrow(couples))
for (i in 1:nrow(couples)){
for (j in i:nrow(couples)){
vec1=c(couples$deb[i]:couples$fin[i]) 
overlap=length(vec1[vec1 %in% c(couples$deb[j]:couples$fin[j])])
if(overlap>0){
subs[i,j]=FALSE
subs[j,i]=FALSE}
}
}
couples$fin[which(couples$fin>365)]=couples$fin[which(couples$fin>365)]-365
colnames(subs)=paste("temp",couples$deb,couples$fin,sep="_")
rownames(subs)=paste("temp",couples$deb,couples$fin,sep="_")
subs[nrow(subs),1]=FALSE
subs[1,nrow(subs)]=FALSE
subs=as.data.frame(subs)
subs$vari=rownames(subs)
subs=melt(subs,id.var="vari")
subs=subset(subs,value=="FALSE")

#create all combinations of 1 and 2 climatic variables
bid1=as.data.frame(t(combn(paste("temp",couples$deb,couples$fin,sep="_"),m=1)))
names(bid1)=c("var1")
bid1$var2=NA
bid1$var3=NA
bid2=as.data.frame(t(combn(paste("temp",couples$deb,couples$fin,sep="_"),m=2)))
names(bid2)=c("var1","var2")
bid2$var3=NA
# bid3=as.data.frame(t(combn(paste("temp",couples$deb,couples$fin,sep="_"),m=3)))
# names(bid3)=c("var1","var2","var3")
bidf=rbind(bid1,bid2)


#removing combination including temperature indices that are overlapping too much
bidf$autor=NA
func=function(y){max(apply(subs,1,function(x)length(which(t(y) %in% t(x)))))}
bidf$autor=apply(bidf,1,func)
bidf=subset(bidf,autor<2)

fwrite(bidf,"combination_for_dredge.txt")
fwrite(donnf2,"data_final.txt")
#################################################
library(data.table)
library(dplyr)
library(glmmTMB)
#library(doMPI)
#library(igraph)
library(doParallel)
library(foreach)
library(parallel)
library(lme4)
library(spaMM)

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])
i=1

nvariables=2


setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")

bidf=fread("combination_for_dredge.txt")
donnf2=fread("data_final.txt")
liste=as.data.frame(unique(donnf2[,c("Speciesgen","Species","FAMILLE","ORDRE")]))
liste[,names(donnf2)[grep("temp",names(donnf2))]]=NA
liste$rsq=NA
liste$nb_pre1990=NA
liste$nb_post1990=NA
liste$nb=NA
liste$inter=NA
liste$sing=NA

bidon=subset(donnf2,Speciesgen==liste$Speciesgen[i])
bidon$maille=as.factor(bidon$maille)
if(length(unique(bidon$maille))<nrow(bidon) & nrow(bidon)>50){
a=0
bidon$elev=(bidon$elev-mean(bidon$elev,na.rm=T))/1000
dredgi=bidf
dredgi$aicc=NA
if(length(unique(bidon$maille))>1){
dredgi=foreach(j=1:2,.combine=rbind)%dopar%{
varia=t(bidf)[1:nvariables,j][which(!is.na(t(bidf)[1:nvariables,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"elev+(",paste(varia,collapse="+"),"1|maille)+(1|Annee)"))
model=fitme(form,data=bidon)
return(cbind(bidf[j,],data.frame(aicc=extractAIC(model)[2])))
}
dredgi=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:nvariables,1][which(!is.na(t(dredgi)[1:nvariables,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"elev+(",paste(varia,collapse="+"),"1|maille)+(1|Annee)"))
dredgi=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:nvariables,1][which(!is.na(t(dredgi)[1:nvariables,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"elev+(",paste(varia,collapse="+"),"1|maille)+(1|Annee)"))
model=fitme(form,data=bidon)
}else{
dredgi=foreach(j=1:nrow(bidf),.combine=rbind)%dopar%{
varia=t(bidf)[1:nvariables,j][which(!is.na(t(bidf)[1:nvariables,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"elev+(1|Annee)"))
model=fitme(form,data=bidon)
return(cbind(bidf[j,],data.frame(aicc=extractAIC(model)[2])))
}
dredgie=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:nvariables,1][which(!is.na(t(dredgi)[1:nvariables,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"elev+(1|Annee)"))
model=fitme(form,data=bidon)
}
bidon$resid=residuals(model)
liste$rsq[i]=NA
liste[i,names(model$fixef[-1])]=model$fixef[-1]
liste$nb_pre1990[i]=nrow(subset(bidon,Annee<1990))
liste$nb_post1990[i]=nrow(subset(bidon,Annee>=1990))
liste$inter[i]=model$fixef[1]
#if(length(unique(bidon$maille))>1){liste$sing[i]=isSingular(model)}else{liste$sing[i]=NA}
}else{
liste$rsq[i]=NA
liste[i,names(tabtemp)[4:ncol(tabtemp)]]=NA
liste$nb[i]=NA
liste$inter[i]=NA
liste$sing[i]=NA
}

fwrite(liste[i,],paste0(liste$Speciesgen[i],"_selec_temp.txt"),sep="\t")


