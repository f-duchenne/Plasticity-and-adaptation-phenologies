library(mgcv)
library(numDeriv)
library(data.table)
library(dplyr)
library(sf)
library(glmmTMB)
library(spaMM)
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
"+scale(Altitude)+Annee2+(1+",paste(varia,collapse="+"),"+Annee2|maille)+AR1(1|Annee)+ Matern(1|Longitude+Latitude)"))
model <-fitme(f.s, data = bidon,family = gaussian(),method="REML")


setwd(dir="/home/duchenne/plast/resultats")
save(model,file=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),".RData"))

#	
bidon$resid=residuals(model)
res=as.data.frame(summary(model)$beta_table)
res$varia=rownames(res)
res$Speciesgen=liste$Speciesgen[i]

resf=res

setwd(dir="/home/duchenne/plast/resultats")
fwrite(resf,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),".txt"),sep="\t")
fwrite(resqual,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),"_qual.txt"),sep="\t")