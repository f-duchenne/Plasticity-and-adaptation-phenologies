library(mgcv)
library(numDeriv)
library(data.table)
library(dplyr)
library(sf)
library(glmmTMB)
library(lme4)
library(e1071)
library(ade4)

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])
print(i)

setwd(dir="/home/duchenne/plast/")
shp <- st_read(".", "Europe_coastline")

temp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")

dat=fread("selec_temp_var_beta.txt",sep="\t")
liste=data.frame(Speciesgen=unique(dat$Speciesgen),nb_pre1990=dat$nb_pre1990[!duplicated(dat$Speciesgen)],nb_post1990=dat$nb_post1990[!duplicated(dat$Speciesgen)])
liste$nb=liste$nb_pre1990+liste$nb_post1990

resf=NULL

bidon=fread(paste0("data_per_species/",liste$Speciesgen[i],".txt"))
bidon=subset(bidon, !is.na(elev) & !is.na(temp_0_90))

bidon$Altitude2=bidon$Altitude
bidon$Altitude=bidon$Altitude
bidon$Annee2=bidon$Annee-mean(bidon$Annee)
bidon$Annee=bidon$Annee-1960
bidon$maille_slope=bidon$maille
varia=subset(dat,Speciesgen==liste$Speciesgen[i])$varia
bidon$loc=paste(bidon$Longitude_10km,bidon$Latitude_10km,sep="_")
bidon$period="pre1990"
bidon$period[(bidon$Annee+1960)>=1990]="post1990"


f.s=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+scale(Altitude)+Annee2+(1+",paste(varia,collapse="+"),"+Annee2|maille)+(1|Annee)"))
model <-lmer(f.s, data = bidon)


setwd(dir="/home/duchenne/plast/resultats")
save(model,file=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),".RData"))

#	
bidon$resid=residuals(model)
res=as.data.frame(summary(model)$coefficients)
res$varia=rownames(res)
res$core=NA
res$core[res$varia %in% varia]=t(cor(bidon[,c("Annee",varia),with=F])[-1,1])
res$STD_maille=NA
res$STD_maille=attr(summary(model)$varcor$maille,"stddev")[res$varia]
res$Speciesgen=liste$Speciesgen[i]
res$latmean=mean(bidon$Y)
res$longmean=mean(bidon$X)
res$Altitudemean=mean(bidon$Altitude2,na.rm=T)
res$mfd=mean(bidon$Jour.de.collecte)
res$nb_annee=length(unique(bidon$Annee))
res$nbdata=nrow(bidon)
bidon$Annee2=plyr::round_any(bidon$Annee,2)
res$skew_pre1990=skewness(bidon$resid[bidon$period=="pre1990"])
res$skew_post1990=skewness(bidon$resid[bidon$period=="post1990"])
res$kurt_pre1990=kurtosis(bidon$resid[bidon$period=="pre1990"])
res$kurt_post1990=kurtosis(bidon$resid[bidon$period=="post1990"])


resqual=as.data.frame(bidon %>% dplyr::group_by(Speciesgen,Species,FAMILLE,ORDRE,maille,ctr_lon,ctr_lat) %>%
dplyr::summarise(n=length(unique(Annee))))
randoms=as.data.frame(ranef(model)$maille)
resqual=cbind(resqual,randoms+t(replicate(nrow(resqual),res[names(randoms),"Estimate"])))
biche=as.data.frame(t(sapply(asplit(attr(ranef(model)$maille,"postVar"), 3),diag)))
names(biche)=sort(names(randoms))
biche=sqrt(biche+t(replicate(nrow(resqual),diag(vcov(model)[names(randoms),names(randoms)]))))
names(biche)=paste(names(biche),"sde",sep="_")
resqual=cbind(resqual,biche)


res$trend=NA
for(j in varia){
formu=formula(paste0(j,"~Annee+(1+Annee|maille)"))
data=subset(temp,Annee>=min(bidon$Annee2+1990) & Annee<=max(bidon$Annee2+1990) & maille %in% bidon$maille)
mod <- lmer(formu,data=data)
res$trend[res$varia==j]=fixef(mod)[2]
resqual[,paste0(j,"_trend")]=ranef(mod)$maille$Annee
}


resf=rbind(resf,res)
resqual=subset(resqual,n>=5)

setwd(dir="/home/duchenne/plast/resultats")
fwrite(resf,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),".txt"),sep="\t")
fwrite(resqual,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),"_qual.txt"),sep="\t")