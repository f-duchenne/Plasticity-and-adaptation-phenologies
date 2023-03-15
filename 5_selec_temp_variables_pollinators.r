library(data.table)
library(dplyr)
#library(doMPI)
#library(igraph)
library(doParallel)
library(foreach)
library(parallel)


setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
donnf2=fread("dataf_avec_elev.txt",sep="\t")
donnf2$Latitude_10km=round(donnf2$Latitude_10km,5)
donnf2$Longitude_10km=round(donnf2$Longitude_10km,5)
tabtemp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")
tabtemp$Latitude_10km=round(tabtemp$Latitude_10km,5)
tabtemp$Longitude_10km=round(tabtemp$Longitude_10km,5)
dim(donnf2)
donnf2=merge(donnf2,tabtemp,by=c("Latitude_10km","Longitude_10km","Annee","maille"))
dim(donnf2)
donnf2=subset(donnf2,!is.na(elev)  & !is.na(temp_0_90))
dim(donnf2)
liste=as.data.frame(unique(donnf2[,c("Speciesgen","Species","FAMILLE","ORDRE")]))
liste[,names(tabtemp)[4:ncol(tabtemp)]]=NA
liste$rsq=NA
liste$nb_pre1990=NA
liste$nb_post1990=NA
liste$nb=NA
liste$inter=NA
liste$sing=NA


#defining temperature indice which are overlapping too much
couples=data.frame(deb=rep(seq(0,365,15),2),fin= c(seq(0,365,15)+30,seq(0,365,15)+90))
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


liste=donnf2 %>% dplyr::group_by(Speciesgen,Species,FAMILLE,ORDRE) %>% dplyr::summarise(nb_pre1990=length(Latitude[Annee<1990]),nb_post1990=length(Latitude[Annee>=1990]))
liste=subset(liste,nb_pre1990>=50 & nb_post1990>=500)
liste[,names(donnf2)[grep("temp",names(donnf2))]]=NA
liste$rsq=NA
liste$inter=NA

fwrite(liste,"liste_species_to_study.txt")
fwrite(bidf,"combination_for_dredge.txt")
#fwrite(donnf2,"data_final.txt")

for(i in 1:nrow(liste)){
bidon=subset(donnf2,Speciesgen==liste$Speciesgen[i])
fwrite(bidon,paste0("data_per_species/",liste$Speciesgen[i],".txt"))
}

#################################################
library(data.table)
library(dplyr)
library(glmmTMB)
library(doParallel)
library(foreach)
library(parallel)
library(lme4)
# setup parallel backend to use all available CPUs
ncpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
registerDoParallel(ncpus)

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])
print(i)

nvariables=2

setwd(dir="/home/duchenne/plast/")
temp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")
liste=as.data.frame(fread("liste_species_to_study.txt"))
bidf=fread("combination_for_dredge.txt")
#bidf=bidf[1:2,]
bidon=fread(paste0("data_per_species/",liste$Speciesgen[i],".txt"))
bidon$maille=as.factor(bidon$maille)


if(length(unique(bidon$maille))<nrow(bidon) & nrow(bidon)>50){
a=0
bidon$Altitude=(bidon$Altitude-mean(bidon$Altitude,na.rm=T))/1000
dredgi=bidf
dredgi$aicc=NA

if(length(unique(bidon$maille))>1){
dredgi=foreach(j=1:nrow(bidf),.combine=rbind)%dopar%{
varia=t(bidf)[1:nvariables,j][which(!is.na(t(bidf)[1:nvariables,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+Altitude+(",paste(varia,collapse="+"),"+1|maille)+(1|Annee)"))
model=lmer(form,data=bidon,REML=F)
return(cbind(bidf[j,],data.frame(aicc=extractAIC(model)[2])))
}
dredgi=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:nvariables,1][which(!is.na(t(dredgi)[1:nvariables,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+Altitude+(",paste(varia,collapse="+"),"+1|maille)+(1|Annee)"))
model=lmer(form,data=bidon)
}else{
dredgi=foreach(j=1:nrow(bidf),.combine=rbind)%dopar%{
varia=t(bidf)[1:nvariables,j][which(!is.na(t(bidf)[1:nvariables,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+Altitude+(1|Annee)"))
model=lmer(form,data=bidon,REML=F)
return(cbind(bidf[j,],data.frame(aicc=extractAIC(model)[2])))
}

dredgie=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:nvariables,1][which(!is.na(t(dredgi)[1:nvariables,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+Altitude+(1|Annee)"))
model=lmer(form,data=bidon)
}

res=cbind(liste[i,c("Speciesgen","Species","FAMILLE","ORDRE","nb_pre1990","nb_post1990")],data.frame(varia=varia))
res$cor=cor(bidon[,c("Annee",varia),with=F])[1,varia]
res$trend=NA
res$trend.se=NA
for(j in varia){
formu=formula(paste0(j,"~Annee+(1+Annee|maille)"))
data=subset(temp,Annee>=min(bidon$Annee) & Annee<=max(bidon$Annee) & maille %in% bidon$maille)
mod <- lmer(formu,data=data)
res$trend[varia==j]=fixef(mod)[2]
res$trend.se[varia==j]=summary(mod)$coefficients["Annee",2]
}
#if(length(unique(bidon$maille))>1){liste$sing[i]=isSingular(model)}else{liste$sing[i]=NA}
}else{
res=cbind(liste[i,c("Speciesgen","Species","FAMILLE","ORDRE","nb_pre1990","nb_post1990")],data.frame(varia=NA))
res$cor=NA
res$trend=NA
res$trend.se=NA
}

setwd(dir="/home/duchenne/plast/species_temp_sel")
fwrite(res,paste0(liste$Speciesgen[i],"_selec_temp.txt"),sep="\t")

stopCluster(ncpus)


