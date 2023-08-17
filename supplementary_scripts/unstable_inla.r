library(mgcv)
library(numDeriv)
library(data.table)
library(dplyr)
library(sf)
library(INLA)
library(lme4)
library(e1071)


setwd(dir="D:/land use change/European_costlines")
shp <- st_read(".", "Europe_coastline")

temp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")
nrep=3

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")

dat=fread("selec_temp_var_beta.txt",sep="\t")
liste=data.frame(Speciesgen=unique(dat$Speciesgen),nb_pre1990=dat$nb_pre1990[!duplicated(dat$Speciesgen)],nb_post1990=dat$nb_post1990[!duplicated(dat$Speciesgen)])
liste$nb=liste$nb_pre1990+liste$nb_post1990
liste=subset(liste,nb<=2000)

resf=NULL
for(i in 1:nrow(liste)){
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

f.s2=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+scale(Altitude)+Annee2+(1+",paste(varia,collapse="+"),"+Annee2|maille)+(1|Annee)"))

for(j in 1:nrep){
model <- inla(f.s, family = "gaussian", data = inla.stack.data(stk.dat), verbose = FALSE, 
    control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE), control.fixed=list(mean=list(beta0=180,default=0),prec=list(Annee2=5,beta0=0.0005,default=0.001)))

	
bidon$resid=bidon$Jour.de.collecte-model$summary.fitted.values$mean[grep("APredictor",rownames(model$summary.fitted.values),fixed=T)]
res_inla=as.data.frame(model$summary.fixed)
res_inla$varia=rownames(res_inla)
res_inla$Speciesgen=liste$Speciesgen[i]
res_inla$essai=j


model <-lmer(f.s2, data = bidon)
res=as.data.frame(summary(model)$coefficients)

resf=rbind(resf,cbind(res_inla,res))

}

}



