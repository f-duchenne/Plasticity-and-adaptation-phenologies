pkgs <- c("ggplot2", "mgcv", "MASS","car","gplots","doBy","dplyr","lme4","mgcv","numDeriv","nleqslv","lubridate","reshape2","Hmisc",
"fitdistrplus","geosphere","plyr","pastecs","SDMtools","jtools","gridExtra","stats","ggmap")
lapply(pkgs, require, character.only = TRUE,quietly=T)
setwd(dir="C:/Users/Francois/Documents/papier 1 - données collection/data")
library(mgcv)
library(numDeriv)
library(nleqslv)
library(lubridate)
library(Hmisc)
library(stats)
library(plyr)
library(pastecs)
library(SDMTools)
library(geosphere)
library(fitdistrplus)
library(reshape2)
library(data.table)
library(glmmTMB)
xmin=-11
xmax=2.3
ymin=49.9
ymax=59.6
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
donnf2=fread("data_final_avec_gen2.txt",sep="\t",header=T)
donnf2$Annee=as.numeric(as.character(donnf2$Annee))
donnf2=donnf2 %>% filter(Latitude>=ymin & Latitude<=ymax & Longitude>=xmin & Longitude<=xmax)
donnf2=donnf2 %>% dplyr::filter(Latitude>50.53 | Longitude<(-0.7))
donnf2=donnf2 %>% dplyr::filter(!(Latitude<51.04 & Longitude>1.2829))
donnf2=donnf2 %>% filter(Annee>=1960)
#plot(unique(donnf2[,c("Longitude2","Latitude2")]))
donnf2$categorie="A"
donnf2$categorie[which(donnf2$Annee>=1990)]="B"
names(donnf2)[4]="Species"
donnf2$Speciesgen=paste(donnf2$Species,donnf2$gen,sep="_")
b=donnf2 %>% dplyr::group_by(ORDRE,Speciesgen,categorie) %>% dplyr::summarise(n=n_distinct(Reference))
b=dcast(b,Speciesgen+ORDRE~categorie,value.var="n")
b$C=b$A+b$B
b2=subset(b,A>=50 & C>=1000)
summary(as.factor(b2$ORDRE))
donnf2=donnf2[donnf2$Speciesgen %in% b2$Speciesgen,]


library(raster)
library(sf)
library(rgdal)
setwd(dir="C:/Users/Francois/Documents/land use change/European_costlines")
shp <- st_read(".", "Europe_coastline")
e <- as(raster::extent(xmin,xmax,ymin,ymax), "SpatialPolygons") %>% 
  st_as_sf() %>% st_set_crs(4326) %>% st_transform(st_crs(shp))
grd_lrg <- st_make_grid(e, cellsize = c(50000,50000))
pts <- do.call(rbind,st_geometry(st_centroid(grd_lrg)))
pts=as.data.frame(pts)
names(pts)=c("ctr_lon","ctr_lat")
pts$maille=1:nrow(pts)
plot(grd_lrg)
plot(shp,add=T)
donnf=st_as_sf(donnf2,coords=c("Longitude","Latitude"))%>%
  st_set_crs(4326) %>% 
  st_transform(st_crs(shp))
random_joined = st_intersects(donnf,grd_lrg)
is.integer0 <- function(x){is.integer(x) && length(x) == 0L}
random_joined[unlist(lapply(random_joined ,is.integer0))] <- NA 
donnf2$maille=NA
donnf2$maille=sapply(random_joined,function(x) x[1])
donnf2=cbind(donnf2,as.data.table(st_coordinates(donnf)))
bidon=subset(donnf2,is.na(maille))
bidon=st_as_sf(bidon,coords=c("Longitude","Latitude"))%>%
  st_set_crs(4326) %>% 
  st_transform(st_crs(shp))
plot(bidon,add=T,pch=19,col="black")
donnf2=subset(donnf2,!is.na(maille))

bidon=data.frame(latitude=plyr::round_any(donnf2$Y,10000),
longitude=plyr::round_any(donnf2$X,10000))
bidon=sf::st_transform(sf::st_as_sf(bidon,coords = c('longitude','latitude'),crs=st_crs(shp)),crs=4326)
bidon2=as.data.frame(st_coordinates(bidon))
donnf2$latitude=bidon2[,2]
donnf2$longitude=bidon2[,1]
donnf2=merge(donnf2,pts,by="maille")
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
fwrite(donnf2,"dataf.txt",sep="\t")
#############################################################################################
####EXTRACTION DES TEMPERATURE###############################
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
#charger le package:
library(climateExtract)
library(plyr)

library(nlme)
library(raster)
library(data.table)
#on ne récupère que les données européennes selon un maillage de 0.3
#SI on veut toute l'Europe chargement du maillage de 0.3:
# a=getData('worldclim', var='tmin',res=10)
# b=rasterToPoints(a,spatial=TRUE)
# coord=as.data.frame(coordinates(b))
# coord=subset(coord,coord$x>(-15) & coord$x<34 & coord$y>32 & coord$y<73) 
# names(coord)=c("longitude","latitude")
# Points=coord
# Points$longitude=round_any(Points$longitude,0.3)
# Points$latitude=round_any(Points$latitude,0.3)
# Points=Points[!duplicated(Points[,c("longitude","latitude")]),]
memory.limit(size = 1e9)
#Que pour notre jeux de données:
donnf2=fread("dataf.txt",sep="\t")
Points=unique(donnf2[,c("longitude","latitude")])
Points$site_id=1:nrow(Points)
Points=Points[,c("site_id","longitude","latitude")]
memory.limit(size=5000000)
years=c(1960,1970,1980,1990,2000,2016)
climate_data <- extract_nc_value(1950,1964,local_file=T,clim_variable='mean temp',grid_size=0.1)
temp=point_grid_extract(climate_data,Points)
climate_data <- extract_nc_value(1965,1980,local_file=T,clim_variable='mean temp',grid_size=0.1)
temp2=point_grid_extract(climate_data,Points)
climate_data <- extract_nc_value(1981,1990,local_file=T,clim_variable='mean temp',grid_size=0.1)
temp3=point_grid_extract(climate_data,Points)
climate_data <- extract_nc_value(1991,2002,local_file=T,clim_variable='mean temp',grid_size=0.1)
temp4=point_grid_extract(climate_data,Points)
climate_data <- extract_nc_value(2003,2016,local_file=T,clim_variable='mean temp',grid_size=0.1)
temp5=point_grid_extract(climate_data,Points)

tempf=rbind(temp,temp2,temp3,temp4,temp5)
fwrite(tempf,"temp_1950_2016.txt",sep="\t",row.names=F)
fwrite(Points,"points_temp_1950_2016.txt",sep="\t",row.names=F)

####CALCUL DE L'INDICE ANNUEL DE TEMPERATURE###############################
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
library(lubridate)
library(doBy)
library(data.table)
library(dplyr)
#charger les tableaux un à un pour éviter de faire péter l'ordi
tab=fread("temp_1950_2016.txt",sep="\t",header=T)
Points=fread("points_temp_1950_2016.txt",sep="\t",header=T)
Points$latitude=round(Points$latitude,5)
Points$longitude=round(Points$longitude,5)
#définir julian day and year:
tab$jour=yday(tab$date_extract)
tab$y=year(tab$date_extract)
couples=data.frame(deb=rep(seq(0,365,30),3),fin= c(seq(0,365,30)+30,seq(0,365,30)+60,seq(0,365,30)+90))
couples$fin[which(couples$fin>365)]=couples$fin[which(couples$fin>365)]-365
couples=data.frame(deb=seq(0,365,30),fin= c(seq(0,365,30)+90))
couples$fin[which(couples$fin>365)]=couples$fin[which(couples$fin>365)]-365

for(i in 1:nrow(couples)){
#subseter la période voulue
if(couples$fin[i]<=couples$deb[i]){tab2=subset(tab,jour>couples$deb[i] | jour<=couples$fin[i])
}else{tab2=subset(tab,jour>couples$deb[i] & jour<=couples$fin[i])}
if(couples$fin[i]<=couples$deb[i]){tab2$y[which(tab2$jour>couples$deb[i])]=tab2$y[which(tab2$jour>couples$deb[i])]+1}
func=function(x){length(x[!is.na(x)])}
func2=function(x){mean(x,na.rm=T)}
obj=tab2[,which(names(tab) %in% c(paste(Points$site_id),"y")),with=F] %>%
dplyr::group_by(y) %>% dplyr::summarise_all(list(func2))
bidon=melt(obj,id.var="y")
bidon=merge(bidon,Points,by.x="variable",by.y="site_id")
bidon$indice=paste("temp",couples$deb[i],couples$fin[i],sep="_")
if(i==1){tabf=bidon}else{tabf=rbind(tabf,bidon)}
print(i)}
names(tabf)[1:2]=c("point","Annee")
tabf=tabf %>% dplyr::group_by(point,indice) %>% dplyr::mutate(moy=mean(value,na.rm=T))
tabf=tabf %>% dplyr::group_by(Annee,indice) %>% filter (!duplicated(latitude,longitude))
tabf$value=tabf$value-tabf$moy
tabf2=reshape2::dcast(tabf,latitude+longitude+Annee~indice,value.var="value")
write.table(tabf2,"annual_mean_by_indice_and_by_point_mat.txt",sep="\t",row.names=F)

######################################
###ADD ELEVATION###
pkgs <- c("ggplot2", "mgcv", "MASS","car","gplots","doBy","dplyr","lme4","mgcv","numDeriv","nleqslv","lubridate","reshape2","Hmisc",
"fitdistrplus","geosphere","plyr","pastecs","SDMtools","jtools","gridExtra","stats","ggmap")
lapply(pkgs, require, character.only = TRUE,quietly=T)
setwd(dir="C:/Users/Francois/Documents/papier 1 - données collection/data")
library(mgcv)
library(numDeriv)
library(nleqslv)
library(lubridate)
library(Hmisc)
library(stats)
library(plyr)
library(pastecs)
library(SDMTools)
library(geosphere)
library(fitdistrplus)
library(reshape2)
library(data.table)
library(glmmTMB)
library(spaMM)
library(raster)
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
donnf2=fread("dataf.txt",sep="\t")

all=as.data.frame(donnf2 %>% dplyr::group_by(X,Y) %>% dplyr::count())
step<-10000
data_size <- length(all$X)
#add altitude, remove it from the following loop to run faster:
setwd(dir="C:/Users/Francois/Documents/land use change/elevation")
elev=raster("elevation1x1_new.tif")
all[,"elev"]=NA

for(i in seq(1,data_size,step)){
  
  i1 <- i
  i2 <- i + (step-1)
  
  if (i2 > data_size) {i2 <- data_size}
  
  print(c(i1,i2))
  
  pts <- all[i1:i2,c('X','Y')]
  
  sf_pts <- sf::st_transform(sf::st_as_sf(pts, coords = c('X','Y'),crs=projection(elev)), raster::projection(elev)) # match the same projection to species obs points and raster
  all[i1:i2,"elev"]=raster::extract(elev, as(sf_pts, 'Spatial'))
 
}

donnf2=merge(as.data.table(donnf2),as.data.table(all),by=c('X','Y')) 
nrow(restot)-nrow(res)
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
fwrite(donnf2,"dataf_avec_elev.txt",sep="\t",row.names=F)
###################################################################################
.libPaths(c("/home/duchenne/R/x86_64-pc-linux-gnu-library/3.4/",.libPaths()))
library(data.table)
library(dplyr)
library(glmmTMB)
library(doMPI)
library(igraph)
library(doParallel)
library(foreach)
library(parallel)
# cl<-startMPIcluster()
# registerDoMPI(cl)
cl<-makeCluster(6)
registerDoParallel(cl)
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
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
# resultat=NULL
#for(i in 1:nrow(liste)){
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
# bidon$pos=numFactor(bidon$X,bidon$Y)
# bidon$group <- factor(rep(1, nrow(bidon)))
# form=formula(paste0("Jour.de.collecte~",paste(names(tabtemp)[4:ncol(tabtemp)],collapse="+"),
# "+(1|maille)+exp(pos+0|group)"))
# model=glmmTMB(form, data=bidon)
# form=formula(paste0("Jour.de.collecte~",paste(names(tabtemp)[4:ncol(tabtemp)],collapse="+"),
# "+elev+(1|maille)")) #+Matern(1 | X + Y)
# model=fitme(form, bidon)
if(length(unique(bidon$maille))>1){
for(j in 1:nrow(dredgi)){
varia=t(dredgi)[1:3,j][which(!is.na(t(dredgi)[1:3,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1+",paste(varia,collapse="+"),"|maille)+(1|Annee)"))
#model=fitme(form,data=bidon)
model=lmer(form,data=bidon,REML=F)
dredgi$aicc[j]=extractAIC(model)[2]
}
dredgi=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:3,1][which(!is.na(t(dredgi)[1:3,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1+",paste(varia,collapse="+"),"|maille)+(1|Annee)"))
#model=fitme(form,data=bidon)
model=lmer(form,data=bidon,REML=T)
}else{
for(j in 1:nrow(dredgi)){
varia=t(dredgi)[1:3,j][which(!is.na(t(dredgi)[1:3,j]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1|Annee)"))
#model=fitme(form,data=bidon)
model=lmer(form,data=bidon,REML=F)
dredgi$aicc[j]=extractAIC(model)[2]
}
dredgie=dredgi[order(dredgi$aicc),]
varia=t(dredgi)[1:3,1][which(!is.na(t(dredgi)[1:3,1]))]
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+(1|Annee)"))
#model=fitme(form,data=bidon)
model=lmer(form,data=bidon,REML=T)
	}
bidon$resid=residuals(model)
# pl1=ggplot(data=bidon,aes(x=X,y=Y,color=sqrt(abs(resid))*sign(resid)))+
# geom_point(size=2)+scale_color_viridis()
suma=summary(model)$coefficients
liste$rsq[i]=NA
liste[i,rownames(suma[-1,])]=suma[-1,1]
liste$nb_pre1990[i]=nrow(subset(bidon,Annee<1990))
liste$nb_post1990[i]=nrow(subset(bidon,Annee>=1990))
liste$inter[i]=suma[1,1]
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

fwrite(resultat,"temperature_select_var.txt",sep="\t")

closeCluster(cl)
#mpi.quit()



##############################################
###################################################################################
library(data.table)
library(dplyr)
library(glmmTMB)
library(igraph)
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
donnf2=fread("dataf_avec_elev.txt",sep="\t")
donnf2$latitude=round(donnf2$latitude,5)
donnf2$longitude=round(donnf2$longitude,5)
tabtemp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")
tabtemp=tabtemp[,-which(names(tabtemp)=="temp_360_85"),with=F]
donnf2=merge(donnf2,tabtemp,by=c("latitude","longitude","Annee"))

donnf2=subset(donnf2,!is.na(elev))
donnf2$cate=as.character(donnf2$cate)
donnf2$cate="A"
donnf2$cate[donnf2$Annee>=1990]="B"
b=donnf2 %>% group_by(Speciesgen,Species,FAMILLE,ORDRE,cate) %>% count()
b=dcast(b,Speciesgen+Species+FAMILLE+ORDRE~cate,value.var="n")
b=subset(b,A>=50 & B>=1000)
summary(as.factor(b$ORDRE))

fwrite(b,"nb_data.txt",sep="\t")



