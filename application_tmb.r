library(mgcv)
library(numDeriv)
library(data.table)
library(dplyr)
library(sf)
library(glmmTMB)
library(lme4)
library(e1071)

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
bidon$Altitude=(bidon$Altitude-mean(bidon$Altitude,na.rm=T))/1000
bidon$Annee2=bidon$Annee-1990
bidon$Annee=numFactor(bidon$Annee)
bidon$group=1
bidon$maille_slope=bidon$maille
bidon$pos <- numFactor(bidon$Longitude, bidon$Latitude)
varia=subset(dat,Speciesgen==liste$Speciesgen[i])$varia


f.s=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+Altitude+Annee+(1+",paste(varia,collapse="+"),"+Annee|maille)+exp(pos + 0 | group)+ou(Annee + 0 | group)"))
model <- glmmTMB(f.s, family = "gaussian", data = bidon)


setwd(dir="/home/duchenne/plast/resultats")
save(model,file=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),".RData"))

summary(model)
#	
#bidon$resid=bidon$Jour.de.collecte-model$summary.fitted.values$mean[grep("APredictor",rownames(model$summary.fitted.values),fixed=T)]
#res=as.data.frame(model$summary.fixed)
#res$varia=rownames(res)
#res$core=NA
#res$core[res$varia %in% varia]=t(cor(bidon[,c("Annee",varia),with=F])[-1,1])
#res$Var=NA
#res$Var[res$varia %in% rownames(model$summary.hyperpar)]=model$summary.hyperpar$mean
#res$Speciesgen=liste$Speciesgen[i]
#res$latmean=mean(bidon$Y)
#res$longmean=mean(bidon$X)
#res$Altitudemean=mean(bidon$Altitude2,na.rm=T)
#res$mfd=mean(bidon$Jour.de.collecte)
#res$nb_annee=length(unique(bidon$Annee))
#res$nbdata=nrow(bidon)
#bidon$Annee2=plyr::round_any(bidon$Annee,2)
#b=bidon %>% dplyr::group_by(Speciesgen,Annee2) %>% dplyr::summarise(skew=skewness(resid),nbd=length(resid),
#kurt=kurtosis(resid))
#modelskew=lm(skew~Annee2,data=b,weights=sqrt(nbd))
#res$trend_skew=NA
#res$trend_skew=modelskew$coeff[2]
#res$trend_skew_err=summary(modelskew)$coeff[2,2]
#modelkurt=lm(kurt~Annee2,data=b,weights=sqrt(nbd))
#res$trend_kurt=modelkurt$coeff[2]
#res$trend_kurt_err=summary(modelkurt)$coeff[2,2]
#
#resqual=as.data.frame(bidon %>% dplyr::group_by(Speciesgen,Species,FAMILLE,ORDRE,maille,ctr_lon,ctr_lat) %>%
#dplyr::summarise(n=length(unique(Annee))))
#resqual=cbind(resqual,data.frame(Annee=model$summary.random$maille_slope$mean+model$summary.fixed["Annee","mean"],Annee_sde=sqrt(model$summary.random$maille_slope$sd^2+model$summary.fixed["Annee","sd"]^2)))
#biche=data.frame(model$summary.random$maille_slope1$mean+model$summary.fixed[varia[1],"mean"],sqrt(model$summary.random$maille_slope1$sd^2+model$summary.fixed[varia[1],"sd"]^2),
#model$summary.random$maille_slope2$mean+model$summary.fixed[varia[2],"mean"],sqrt(model$summary.random$maille_slope2$sd^2+model$summary.fixed[varia[2],"sd"]^2))
#names(biche)=paste0(rep(varia,each=2),c("","_sde"))
#biche=biche[,!is.na(names(biche))]
#resqual=cbind(resqual,biche)
#
#
#res$trend=NA
#for(j in varia){
#formu=formula(paste0(j,"~Annee+(1+Annee|maille)"))
#data=subset(temp,Annee>=min(bidon$Annee2+1990) & Annee<=max(bidon$Annee2+1990) & maille %in% bidon$maille)
#mod <- lmer(formu,data=data)
#res$trend[res$varia==j]=fixef(mod)[2]
#resqual[,paste0(j,"_trend")]=ranef(mod)$maille$Annee
#}
#
#
#resf=rbind(resf,res)
#resqual=subset(resqual,n>=5)
#
#setwd(dir="/home/duchenne/plast/resultats")
#fwrite(resf,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),".txt"),sep="\t")
#fwrite(resqual,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),"_qual.txt"),sep="\t")