library(mgcv)
library(numDeriv)
library(data.table)
library(dplyr)
library(sf)
library(INLA)
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
bidon$Altitude=bidon$Altitude
bidon$Annee2=bidon$Annee-mean(bidon$Annee)
bidon$Annee=bidon$Annee-1960
bidon$maille_slope=bidon$maille
varia=subset(dat,Speciesgen==liste$Speciesgen[i])$varia
bidon$loc=paste(bidon$Longitude_10km,bidon$Latitude_10km,sep="_")
bidon$period="pre1990"
bidon$period[(bidon$Annee+1960)>=1990]="post1990"

listevar=list(beta0=rep(1,nrow(bidon)),Annee = bidon$Annee,Annee2=bidon$Annee2,maille=bidon$maille,
maille_slope=bidon$maille_slope,maille_slope1=bidon$maille_slope,maille_slope2=bidon$maille_slope,
maille_slope3=bidon$maille_slope,
Altitude = bidon$Altitude)
b=length(listevar)
for(j in 1:length(varia)){
listevar[j+b]=c(bidon[,paste(varia[j]),with=F])
names(listevar)[j+b]=varia[j]
}


sfpts=sf::st_as_sf(bidon, coords = c("X","Y"),crs=st_crs(shp))
spts=as(sfpts, "Spatial")

# C-LM bounday as SpatialPolygons
prdomain <- inla.nonconvex.hull(as.matrix(bidon[,c("X","Y")]), -0.03, -0.05, resolution = c(100, 100))
if(nrow(bidon)<20000){clm.mesh <- inla.mesh.2d(loc=spts,boundary=prdomain, max.edge = c(100000, 100000),cutoff=1)
}else{
clm.mesh <- inla.mesh.2d(loc=spts,boundary=prdomain, max.edge = c(100000, 100000),cutoff=1000)
}

#par(mar = c(0, 0, 0, 0))
#plot(clm.mesh, asp = 1, main = "")
#xmin=-11
#xmax=2.3
#ymin=49.9
#ymax=59.6
#e <- as(raster::extent(xmin,xmax,ymin,ymax), "SpatialPolygons") %>% 
#  st_as_sf() %>% st_set_crs(4326) %>% st_transform(st_crs(shp))
#grd_lrg <- st_make_grid(e, cellsize = c(50000,50000))
#plot(grd_lrg,add=T)

  
A <- inla.spde.make.A(clm.mesh, loc = as.matrix(bidon[,c("X","Y")]))

spde <- inla.spde2.matern(clm.mesh, alpha = 2)
mesh.index <- inla.spde.make.index(name = "field", n.spde = spde$n.spde)
stk.dat <- inla.stack(data = list(y = bidon$Jour.de.collecte), A = list(A, 1), tag = "est", effects = list(list(field=mesh.index$field),
listevar))
f.s=formula(paste0("y~0+beta0+",paste(varia,collapse="+"),
"+Altitude+Annee2+f(maille,model = 'iid')+f(maille_slope,Annee2,model = 'iid')+",paste(paste0("f(maille_slope",c(1:length(varia)),",",varia,",model='iid')"),collapse="+"),
"+f(field,model=spde)+f(Annee,model='ar1')"))
model <- inla(f.s, family = "gaussian", data = inla.stack.data(stk.dat), verbose = FALSE, 
    control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE), control.fixed=list(mean=list(beta0=180,default=0),prec=list(Annee2=5,beta0=0.0005,default=0.001)))


setwd(dir="/home/duchenne/plast/resultats")
#save(model,file=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),".RData"))
	
bidon$resid=bidon$Jour.de.collecte-model$summary.fitted.values$mean[grep("APredictor",rownames(model$summary.fitted.values),fixed=T)]
res_inla=as.data.frame(model$summary.fixed)
res_inla$varia=rownames(res_inla)
res_inla$Speciesgen=liste$Speciesgen[i]


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


res=cbind(res,res_inla)

resf=rbind(resf,res)
resqual=subset(resqual,n>=5)

setwd(dir="/home/duchenne/plast/resultats")
fwrite(resf,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),".txt"),sep="\t")
fwrite(resqual,paste0(gsub(" ","_",liste$Speciesgen[i],fixed=T),"_qual.txt"),sep="\t")