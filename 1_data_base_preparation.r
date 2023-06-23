library(lubridate)
library(plyr)
library(ggmap)
library(SDMTools)
library(geosphere)
library(reshape2)
library(data.table)
library(raster)
library(sf)
library(rgdal)
xmin=-11
xmax=2.3
ymin=49.9
ymax=59.6
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
donnf2=fread("data_final_avec_gen2.txt",sep="\t",header=T)
donnf2=donnf2[,-c(1:2)]
names(donnf2)[c(2,22)]=c("Species","Speciesgen")
donnf2$Annee=as.numeric(as.character(donnf2$Annee))
donnf2=donnf2 %>% filter(Latitude>=ymin & Latitude<=ymax & Longitude>=xmin & Longitude<=xmax)
donnf2=donnf2 %>% dplyr::filter(Latitude>50.53 | Longitude<(-0.7))
donnf2=donnf2 %>% dplyr::filter(!(Latitude<51.04 & Longitude>1.2829))
donnf2=donnf2 %>% filter(Annee>=1960)
#plot(unique(donnf2[,c("Longitude2","Latitude2")]))
donnf2$categorie="A"
donnf2$categorie[which(donnf2$Annee>=1990)]="B"
donnf2$Speciesgen=paste(donnf2$Species,donnf2$gen,sep="_")
b=donnf2 %>% dplyr::group_by(ORDRE,Speciesgen,categorie) %>% dplyr::summarise(n=n_distinct(Reference))
b=dcast(b,Speciesgen+ORDRE~categorie,value.var="n")
b$C=b$A+b$B
b2=subset(b,A>=50 & C>=1000)
summary(as.factor(b2$ORDRE))
donnf2=donnf2[donnf2$Speciesgen %in% b2$Speciesgen,]



setwd(dir="D:/land use change/European_costlines")
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

bidon=data.frame(latitude_10km=plyr::round_any(donnf2$Y,10000),
longitude_10km=plyr::round_any(donnf2$X,10000))
bidon=sf::st_transform(sf::st_as_sf(bidon,coords = c('longitude_10km','latitude_10km'),crs=st_crs(shp)),crs=4326)
bidon2=as.data.frame(st_coordinates(bidon))
donnf2$Latitude_10km=bidon2[,2]
donnf2$Longitude_10km=bidon2[,1]
donnf2=merge(donnf2,pts,by="maille")
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
fwrite(donnf2,"dataf.txt",sep="\t")
#############################################################################################
####EXTRACTION DES TEMPERATURE###############################
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
#charger le package:
library(climateExtract)
library(plyr)
library(raster)
library(data.table)
memory.limit(size = 1e9)
#Que pour notre jeux de données:
donnf2=fread("dataf.txt",sep="\t")
Points=unique(donnf2[,c("maille","Longitude_10km","Latitude_10km")])
Points$site_id=1:nrow(Points)
Points=Points[,c("site_id","maille","Longitude_10km","Latitude_10km")]
space= sf::st_as_sf(Points,coords=c('Longitude_10km','Latitude_10km'),crs="+proj=longlat +datum=WGS84 +no_defs")
memory.limit(size=5000000)
years=c(1960,1970,1980,1990,2000,2016)
tempf=NULL
for(i in 1950:2016){
climate_data <- extract_nc_value(i,i,local_file=T,clim_variable='mean temp',grid_size=0.1,file_path="D:/land use change/tg_ens_mean_0.1deg_reg_v23.1e.nc",
spatial_extent = space)
aggrega=temporal_aggregate(x = climate_data,agg_function = "mean",variable_name = "average temp",time_step = "window",
win_length = 1)
bidon=cbind(as.data.frame(Points),extract(aggrega,space,df=TRUE))
bidon[is.na(bidon[,paste0("X",i,".01.01")]),grep("X",names(bidon))]=cbind(extract(aggrega,space[is.na(bidon[,paste0("X",i,".01.01")]),],
buffer=20000,fun=mean,na.rm=T))
bidon2=melt(bidon,id.vars=c("site_id","maille","Longitude_10km","Latitude_10km","ID"))
names(bidon2)[names(bidon2)=="variable"]="date_extract"
bidon2$date_extract=as.Date(gsub("X","",bidon2$date_extract),format="%Y.%m.%d")
tempf=rbind(tempf,bidon2)
print(i)
}

fwrite(tempf,"temp_1950_2016.txt",sep="\t",row.names=F)
fwrite(Points,"points_temp_1950_2016.txt",sep="\t",row.names=F)

####CALCUL DE L'INDICE ANNUEL DE TEMPERATURE###############################
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
library(lubridate)
library(data.table)
library(dplyr)
memory.limit(size=5000000)
#charger les tableaux un à un pour éviter de faire péter l'ordi
tab=fread("temp_1950_2016.txt",sep="\t",header=T)
Points=fread("points_temp_1950_2016.txt",sep="\t",header=T)
Points$Latitude_10km=round(Points$Latitude_10km,5)
Points$Longitude_10km=round(Points$Longitude_10km,5)
#définir julian day and year:
tab$jour=yday(tab$date_extract)
tab$y=year(tab$date_extract)
couples=data.frame(deb=rep(seq(0,365,15),3),fin= c(seq(0,365,15)+30,seq(0,365,15)+60,seq(0,365,15)+90))
couples$fin[which(couples$fin>365)]=couples$fin[which(couples$fin>365)]-365
couples=subset(couples,deb<360)

for(i in 1:nrow(couples)){
#subseter la période voulue
if(couples$fin[i]<=couples$deb[i]){
tab2=subset(tab,jour>couples$deb[i] | jour<=couples$fin[i])
}else{
tab2=subset(tab,jour>couples$deb[i] & jour<=couples$fin[i])}

if(couples$fin[i]<=couples$deb[i]){
tab2$y[which(tab2$jour>couples$deb[i])]=tab2$y[which(tab2$jour>couples$deb[i])]+1
}

bidon=tab2 %>% dplyr::group_by(site_id,y,maille,Longitude_10km,Latitude_10km) %>% dplyr::summarise(value=mean(value,na.rm=T))
bidon$indice=paste("temp",couples$deb[i],couples$fin[i],sep="_")
if(i==1){tabf=bidon}else{tabf=rbind(tabf,bidon)}
print(i)}
names(tabf)[1:2]=c("point","Annee")
tabf=tabf %>% dplyr::group_by(point,indice,maille) %>% dplyr::mutate(moy=mean(value,na.rm=T))
tabf=tabf %>% dplyr::group_by(Annee,indice,maille) %>% filter (!duplicated(Latitude_10km,Longitude_10km))
tabf$value=tabf$value-tabf$moy
tabf2=reshape2::dcast(tabf,Latitude_10km+Longitude_10km+maille+Annee~indice,value.var="value")
write.table(tabf2,"annual_mean_by_indice_and_by_point_mat.txt",sep="\t",row.names=F)

######################################
###ADD ELEVATION###
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
library(mgcv)
library(lubridate)
library(plyr)
library(SDMTools)
library(geosphere)
library(reshape2)
library(data.table)
library(dplyr)
library(raster)
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
donnf2=fread("dataf.txt",sep="\t")
donnf2$Latitude=round(donnf2$Latitude,5)
donnf2$Longitude=round(donnf2$Longitude,5)

all=as.data.frame(donnf2 %>% dplyr::group_by(Longitude,Latitude) %>% dplyr::count())
step<-10000
data_size <- nrow(all)
#add altitude, remove it from the following loop to run faster:
setwd(dir="D:/land use change/elevation")
elev=raster("elevation1x1_new.tif")
all[,"elev"]=NA

for(i in seq(1,data_size,by=step)){
  
  i1 <- i
  i2 <- i + (step-1)
  
  if (i2 > data_size) {i2 <- data_size}
  
  print(c(i1,i2))
  
  pts <- all[i1:i2,c('Longitude','Latitude')]
  
  sf_pts <- sf::st_transform(sf::st_as_sf(pts, coords = c('Longitude','Latitude'),crs="+proj=longlat +datum=WGS84 +no_defs"), raster::projection(elev)) # match the same projection to species obs points and raster
  all[i1:i2,"elev"]=raster::extract(elev, as(sf_pts, 'Spatial'))
 
}

nrow(donnf2)
donnf2=merge(as.data.table(donnf2),as.data.table(all),by=c('Longitude','Latitude')) 
nrow(donnf2)
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
fwrite(donnf2,"dataf_avec_elev.txt",sep="\t",row.names=F)
###################################################################################

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



