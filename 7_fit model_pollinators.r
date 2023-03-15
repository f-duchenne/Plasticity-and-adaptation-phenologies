library(mgcv)
library(numDeriv)
library(nleqslv)
library(data.table)
library(dplyr)
library(glmmTMB)
library(doParallel)
library(foreach)
library(parallel)
library(lme4)
library(spaMM)
library(sf)
library(INLA)
library(inlabru)
library(e1071)
setwd(dir="D:/land use change/European_costlines")
shp <- st_read(".", "Europe_coastline")
# setup parallel backend to use all available CPUs
# ncpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
# registerDoParallel(ncpus)

# # Collect command arguments
# args <- commandArgs(trailingOnly = TRUE)
# args_contents <- strsplit(args, ' ')
# # Get first argument
# i <- as.numeric(args_contents[[1]])
# print(i)

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
dat=fread("selec_temp_var_beta.txt",sep="\t")
liste=data.frame(Speciesgen=unique(dat$Speciesgen),nb_pre1990=dat$nb_pre1990[!duplicated(dat$Speciesgen)],nb_post1990=dat$nb_post1990[!duplicated(dat$Speciesgen)])
liste$nb=liste$nb_pre1990+liste$nb_post1990

resf=NULL
for(i in 1:nrow(liste)){
i=which(liste$Speciesgen=="Acasis viretata_2")
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
bidon=fread(paste0("data_per_species/",liste$Speciesgen[i],".txt"))
bidon=subset(bidon, !is.na(elev) & !is.na(temp_0_90))

bidon$Altitude2=bidon$Altitude
bidon$Altitude=(bidon$Altitude-mean(bidon$Altitude,na.rm=T))/1000
bidon$Annee2=bidon$Annee-1990
bidon$Annee=bidon$Annee-1960
bidon$maille_slope=bidon$maille
varia=subset(dat,Speciesgen==liste$Speciesgen[i])$varia

listevar=list(Annee = bidon$Annee,Annee2=bidon$Annee2,maille=bidon$maille,
maille_slope=bidon$maille_slope,maille_slope1=bidon$maille_slope,maille_slope2=bidon$maille_slope,
maille_slope3=bidon$maille_slope,
Altitude = bidon$Altitude,beta0=rep(1,nrow(bidon)))
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
"+Altitude+Annee+f(maille,model = 'iid')+f(maille_slope,Annee,model = 'iid')+",paste(paste0("f(maille_slope",c(1:length(varia)),",",varia,",model='iid')"),collapse="+"),
"+f(field,model=spde)+f(Annee2,model='rw1')"))
model <- inla(f.s, family = "gaussian", data = inla.stack.data(stk.dat), verbose = FALSE, 
    control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE))
print(model$summary.fixed["Annee",])

f.s2=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+Altitude+Annee2+",paste0("(1+Annee2+",paste0(varia,collapse="+"),"|maille)"),
"+(1|Annee)"))
model2=lmer(f.s2,data=bidon)


model2 <- lme(formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),"+Altitude+Annee2"))  , random = formula(paste0("~(1+Annee2+",paste0(varia,collapse="+"),")|maille")), 
                     data = bidon,correlation = corExp(form = ~ Longitude + Latitude))

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
save(model,file=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),".RData"))
	
bidon$resid=bidon$Jour.de.collecte-model$summary.fitted.values$mean[grep("APredictor",rownames(model$summary.fitted.values),fixed=T)]
res=as.data.frame(model$summary.fixed)
res$varia=rownames(res)
res$core=NA
res$core[res$varia %in% varia]=t(cor(bidon[,c("Annee",varia),with=F])[-1,1])
res$Var=NA
res$Var[res$varia %in% rownames(model$summary.hyperpar)]=model$summary.hyperpar$mean
res$Speciesgen=liste$Speciesgen[i]
res$latmean=mean(bidon$Y)
res$longmean=mean(bidon$X)
res$Altitudemean=mean(bidon$Altitude2,na.rm=T)
res$mfd=mean(bidon$Jour.de.collecte)
res$nb_annee=length(unique(bidon$Annee))
res$nbdata=nrow(bidon)
bidon$Annee2=plyr::round_any(bidon$Annee,2)
b=bidon %>% dplyr::group_by(Speciesgen,Annee2) %>% dplyr::summarise(skew=skewness(resid),nbd=length(resid),
kurt=kurtosis(resid))
modelskew=lm(skew~Annee2,data=b,weights=sqrt(nbd))
res$trend_skew=NA
res$trend_skew=modelskew$coeff[2]
res$trend_skew_err=summary(modelskew)$coeff[2,2]
modelkurt=lm(kurt~Annee2,data=b,weights=sqrt(nbd))
res$trend_kurt=modelkurt$coeff[2]
res$trend_kurt_err=summary(modelkurt)$coeff[2,2]

resqual=as.data.frame(bidon %>% dplyr::group_by(Speciesgen,Species,FAMILLE,ORDRE,maille,ctr_lon,ctr_lat) %>%
dplyr::summarise(n=length(unique(Annee))))

resqual=cbind(resqual,data.frame(Annee=model$summary.random$maille_slope$mean+model$summary.fixed["Annee","mean"],Annee_sde=sqrt(model$summary.random$maille_slope$sd^2+model$summary.fixed["Annee","sd"]^2)))
biche=data.frame(model$summary.random$maille_slope1$mean+model$summary.fixed[varia[1],"mean"],sqrt(model$summary.random$maille_slope1$sd^2+model$summary.fixed[varia[1],"sd"]^2),
model$summary.random$maille_slope2$mean+model$summary.fixed[varia[2],"mean"],sqrt(model$summary.random$maille_slope2$sd^2+model$summary.fixed[varia[2],"sd"]^2))
names(biche)=paste0(rep(varia,each=2),c("","_sde"))
biche=biche[,!is.na(names(biche))]
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
}

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
fwrite(resf,"plast_adapt.txt",sep="\t")

######################################
pkgs <- c("ggplot2", "mgcv", "MASS","car","gplots","doBy","dplyr","lme4","mgcv","numDeriv","nleqslv","lubridate","reshape2","Hmisc",
"fitdistrplus","geosphere","pastecs","SDMtools","jtools","gridExtra","stats","ggmap")
lapply(pkgs, require, character.only = TRUE,quietly=T)
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
library(nleqslv)
library(geosphere)
library(fitdistrplus)
library(data.table)
library(glmmTMB)
library(spaMM)
library(viridis)
library(MuMIn)
library(sf)
library(ggeffects)
library(RColorBrewer)
library(ggraph)
library(igraph)


setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
liste=fread("selec_temp_var_beta.txt",sep="\t")
liste=subset(liste,nb_pre1990>=50 & nb_post1990>=500)
vec=names(liste)[grep("temp",names(liste),fixed=T)]
nb=liste[,c("Speciesgen","Species","FAMILLE","ORDRE","nb_pre1990","nb_post1990")]

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
resf=NULL
lili=list.files()
lili=lili[grep("_qual.txt",lili,invert=T)]
for(i in lili){
bidon=fread(i,sep="\t")
resf=rbind(resf,bidon)
}

plot(mean~Estimate,data=resf[resf$varia %in% c("Annee2"),])
abline(0,1)
cor(resf[(resf$varia %in% c("Annee2")),c("mean","Estimate")])

plot(mean~Estimate,data=resf[grep("temp_",resf$varia),])
cor(resf[grep("temp_",resf$varia),c("mean","Estimate")])


resf$type=NA
resf$type[grep("temp",resf$varia)]="plast"
resf$type[grep("Annee",resf$varia)]="adapt"
names(resf)[3:5]=c("lwr","median","upr")
resf$Est.signi=">0.05"
resf$Est.signi[resf$lwr>0]="<0.05"
resf$Est.signi[resf$upr<0]="<0.05"

coucou=c("magenta2",brewer.pal(5,"Set2")[5],"red","dodgerblue4")

####### ADAPT ######
resf$Est.signi=factor(resf$Est.signi,c(">0.05","<0.05"))
b=subset(resf,type=="adapt" & Est.signi=="<0.05") %>% group_by(Speciesgen) %>% dplyr::count()
p=ggplot(data=subset(resf,type=="adapt"),aes(x=mean,fill=Est.signi))+geom_histogram(col="black")+
geom_vline(xintercept=0,col="red")+
scale_fill_manual(values=c("white","black"))+theme_bw()+
theme(panel.grid.minor=element_blank(),strip.background = element_blank(),
legend.title = element_blank(),legend.position="none",plot.title=element_text(size=14,face="bold"))+
xlab("Year effect on mean flight date (day/year)")+ylab("Number of species")+ggtitle("a")

bidon=subset(resf,type=="adapt")
model=lmer(Estimate~longmean2+latmean2+elev+mfd+(1|ORDRE/FAMILLE/Species),
data=bidon,weights=sqrt(1/Cond.SE))
b=ggeffect(model,c("mfd"))
pl1=ggplot(b,aes(x=x,y=predicted))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_point(data=bidon,aes(x=mfd,y=Estimate,col=ORDRE),alpha=0.4)+geom_line()+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("b")+xlab("Average mean flight date (day of the year)")+
ylab("Year effect on mean flight date (day/year)")+labs(color="")+
scale_color_manual(values=coucou)

grid.arrange(p,pl1,widths=c(1,1.5),ncol=2)

png(paste0("fig_adapt_season.png"),width=1400,height=600,res=160)
grid.arrange(p,pl1,widths=c(1,1.5),ncol=2)
dev.off();


####### CONTRIBUTIONS ######
resf$contrib=NA
resf[,contrib:=Estimate*trend]
resf[type=="adapt",contrib:=Estimate]
b=subset(resf,!is.na(type)) %>% dplyr::group_by(Species,ORDRE,FAMILLE,type,elev,mfd,
longmean2,latmean2) %>%
dplyr::summarise(contrib=sum(contrib))
b=b %>% dplyr::group_by(Species) %>% dplyr::mutate(contrib_tot=sum(abs(contrib)),resul=sum(contrib))
b$contrib_perc=b$contrib/b$contrib_tot
b=b[order(b$mfd),]
b$Species=factor(b$Species,unique(b$Species))
bmaster=b


pep=ggplot(data=b,aes(x=type,y=contrib,col=ORDRE))+geom_hline(yintercept=0)+geom_boxplot()+
scale_color_brewer(palette="Set1")+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.title.x=element_blank())+scale_x_discrete(labels=c("Year effect","Phenotypic plasticity"))+
ylab("Contribution to mean flight date shift (day/year)")+labs(color="")

bidon=subset(bmaster,type=="adapt")
bidon %>% dplyr::group_by(ORDRE) %>% dplyr::summarise(moy=mean(abs(contrib_perc)),
sde=sd(abs(contrib_perc))/sqrt(length(contrib_perc)))
model=glmer(abs(contrib_perc)~longmean2+latmean2+elev+mfd+(1|ORDRE/FAMILLE/Species),
data=bidon,family=binomial)
b=ggeffect(model,c("mfd"))
pl2=ggplot(b,aes(x=x,y=predicted))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_point(data=bidon,aes(x=mfd,y=abs(contrib_perc),col=ORDRE),alpha=0.4)+geom_line()+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("c")+xlab("Average mean flight date\n")+
ylab("Relative contribution of year effect to MFD shifts")+labs(color="")+
scale_color_manual(values=coucou)+scale_y_continuous(labels=scales::percent)

b=ggeffect(model,c("elev"))
pl3=ggplot(b,aes(x=x,y=predicted))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_point(data=bidon,aes(x=elev,y=contrib_perc,col=ORDRE),alpha=0.4)+geom_line()+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("d")+xlab("Average mean flight date\n")+
ylab("Mean growth rate")+labs(color="")+
scale_color_brewer(palette="Set1")



pl1=ggplot(data=bmaster,aes(x=Species,y=contrib,fill=type))+geom_bar(stat="identity")+
#geom_point(aes(y=resul))+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_blank(),panel.grid=element_blank(),axis.title=element_text(size=18),axis.text=element_text(size=16),legend.text=element_text(size=16))+
ggtitle("a")+xlab("Phenologies (ordered by mean flight date)")+
ylab("Contribution to mean flight date shifts\n(day/year)")+labs(color="")+
scale_fill_manual(values=c("dodgerblue4","gold3"),labels=c("Year effect","Plasticity"))+
guides(fill = guide_legend(override.aes = list(shape = NA),title=""))

pdf(paste0("fig_diapo.pdf"),width=11,height=5)
ggplot(data=bmaster,aes(x=Species,y=contrib,fill=type))+geom_bar(stat="identity")+
#geom_point(aes(y=resul))+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_blank(),panel.grid=element_blank(),axis.title=element_text(size=18),axis.text=element_text(size=16),legend.text=element_text(size=16))+
ggtitle("")+xlab("Phenologies (ordered by mean flight date)")+
ylab("Contribution to mean flight date shifts\n(day/year)")+labs(color="")+
scale_fill_manual(values=c("dodgerblue4","gold3"),labels=c('"Evolution"',"Plasticity"))+
guides(fill = guide_legend(override.aes = list(shape = NA),title=""))
dev.off();




bidon=dcast(bmaster,Species+ORDRE+FAMILLE+mfd+resul+
longmean2+latmean2~type,value.var="contrib")
model=lm(resul~scale(plast)+scale(adapt),data=bidon)
anova(model)
model=lm(resul~scale(adapt)+scale(plast),data=bidon)
anova(model)

bidon$plast2=abs(bidon$plast)
obj=expand.grid(seq(min(bidon$plast)-0.05,max(bidon$plast)+0.05,0.005),seq(min(bidon$adapt)-0.05,max(bidon$adapt)+0.05,0.005))
obj$value=obj[,1]+obj[,2]
pl1b=ggplot(data=obj,aes(y=Var1,x=Var2,fill=value))+geom_raster()+
geom_point(data=bidon,aes(y=plast,x=adapt),fill="black")+geom_vline(xintercept=0)+
geom_hline(yintercept=0)+geom_abline(intercept=0,slope=-1,linetype="dashed")+
scale_fill_gradient2()+labs(fill="Resulting\nMFD shift\n(day/year)")+theme_bw()+
theme(plot.title=element_text(size=14,face="bold"),panel.background=element_blank(),panel.grid=element_blank(),
panel.margin=unit(c(0,0,0,0), "cm"),panel.spacing = unit(0, "cm"))+
coord_fixed(ratio=1,expand = FALSE)+
ylab("Contribution of phenotypic plasticity\nto MFD shifts (day/year)")+
xlab("Contribution of year effect\nto MFD shifts (day/year)")+ggtitle("b")

grid.arrange(pl1,pl1b,pl2,layout_matrix=matrix(c(1,2,1,3),nrow=2))

png(paste0("contributions.png"),width=1500,height=1200,res=160)
grid.arrange(pl1,pl1b,pl2,layout_matrix=matrix(c(1,2,1,3),nrow=2))
dev.off();


pl1=ggplot(data=bidon,aes(x=plast,y=resul))+geom_point()+
stat_smooth(method="lm",col="black",fill="lightgrey")+
labs(fill="Resulting MFD shift\n(day/year)")+theme_bw()+
theme(plot.title=element_text(size=14,face="bold"),panel.background=element_blank(),panel.grid=element_blank())+
ylab("Long-term MFD shifts (day/year)")+
xlab("Contribution of phenotypic plasticity\nto MFD shifts (day/year)")+ggtitle("a")

pl2=ggplot(data=bidon,aes(x=adapt,y=resul))+geom_point()+
stat_smooth(method="lm",col="black",fill="lightgrey")+
labs(fill="Resulting MFD shift\n(day/year)")+theme_bw()+
theme(plot.title=element_text(size=14,face="bold"),panel.background=element_blank(),panel.grid=element_blank())+
ylab("Long-term MFD shifts (day/year)")+
xlab("Contribution of year effect\nto MFD shifts (day/year)")+ggtitle("b")


grid.arrange(pl1,pl2,nrow=1)

####par maille#
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/résultats simus")
dat=fread("trend_plasticite.txt")
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
tab=fread("plast_adapt_par_maille.txt",header=T)
tab=tab[tab$Species %in% dat$Esp,]
setwd(dir="C:/Users/Francois/Documents/land use change/European_costlines")
shp <- st_read(".", "Europe_coastline")
bidon=sf::st_transform(sf::st_as_sf(tab,coords = c('ctr_lon','ctr_lat'),crs=st_crs(shp)),crs=4326)
bidon=as.data.frame(st_coordinates(bidon))
tab$latitude=bidon[,2]
tab$longitude=bidon[,1]
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
tab=tab[,c("Speciesgen","Species","FAMILLE","ORDRE","maille","ctr_lon","ctr_lat","latitude","longitude",
"variable","value"),with=F]
nb=fread("nb_data_par_maille.txt",header=T)
tab=merge(tab,nb,by=c("Speciesgen","Species","FAMILLE","ORDRE","maille","ctr_lon","ctr_lat"))
tab=tab[grep("_NA",tab$Speciesgen),]
bidon_ann=subset(tab,variable=="Annee")
trendstemp=tab[grep("trend_temp",tab$variable),]
tab=tab[grep("trend_temp",tab$variable,invert=T),]
#tab=merge(tab,unique(resf[,c("Species","mfd")]),by="Species")



names(trendstemp)[names(trendstemp)=="value"]="trends"
trendstemp$variable=gsub("trend_","",trendstemp$variable)
bidon=tab[grep("temp",tab$variable),]
bidon2=merge(bidon,trendstemp,by=c("Speciesgen","Species","FAMILLE","ORDRE","maille","ctr_lon","ctr_lat",
"n","ny","variable","latitude","longitude"))
bidon2$type="plast"
bidon_ann$trends=1
bidon_ann$type="adapt"
bidonf=rbind(bidon2,bidon_ann[,names(bidon2),with=F])

bidonf$contrib=bidonf$value*bidonf$trend
bidonf$variable=gsub("temp_","",bidonf$variable)
bidonf$variable=gsub("_","-",bidonf$variable)
bidonf$date_fin=as.numeric(sapply(strsplit(as.character(bidonf$variable),"-"), function(l) l[[length(l)]]))


b=subset(bidonf,n>=30 & ny>=10) %>% dplyr::group_by(Species,ORDRE,FAMILLE,maille,ctr_lon,ctr_lat,n,ny,type,latitude,
longitude) %>%
dplyr::summarise(contrib=sum(contrib),est_sum=sum(value),trenda=mean(trends))
b=b %>% dplyr::group_by(maille,ctr_lon,ctr_lat,type) %>%
dplyr::mutate(nsp=n_distinct(Species))
b=b %>% dplyr::group_by(Species,ORDRE,FAMILLE,maille,ctr_lon,ctr_lat,n,ny,latitude,
longitude) %>%
dplyr::mutate(contrib_tot=sum(abs(contrib)))
b$contrib_perc=b$contrib/b$contrib_tot

liste=unique(b$Species)
for(i in 1:length(liste)){
bidon=subset(b,Species==liste[i])
model=lm(contrib~ctr_lat+ctr_lon,data=subset(bidon,type=="plast"))
res=as.data.frame(summary(model)$coeff)
res$varia=rownames(res)
res$Species=liste[i]
res$type="plast"
model=lm(contrib~ctr_lat+ctr_lon,data=subset(bidon,type=="adapt"))
res2=as.data.frame(summary(model)$coeff)
res2$varia=rownames(res2)
res2$Species=liste[i]
res2$type="adapt"
res=rbind(res,res2)
if(i==1){resf=res}else{resf=rbind(res,resf)}
}

names(resf)[2:4]=c("sde","t_val","p_val")
dim(subset(resf,varia=="ctr_lat" & p_val<0.05 & Estimate<0 & type=="plast"))
dim(subset(resf,varia=="ctr_lat" & p_val<0.05 & Estimate>0 & type=="plast"))
dim(subset(resf,varia=="ctr_lat" & p_val<0.05 & Estimate<0 & type=="adapt"))
dim(subset(resf,varia=="ctr_lat" & p_val<0.05 & Estimate>0 & type=="adapt"))


model=lmer(contrib~latitude*longitude+(1|Species),data=subset(b,type=="plast"))
pred=ggpredict(model,c("longitude","latitude"))
ggplot(data=pred,aes(x=x,y=predicted,color=group))+geom_line()

model2=lmer(contrib~latitude*longitude+(1|Species),data=subset(b,type=="adapt"))
pred=ggpredict(model2,c("latitude","longitude"))
ggplot(data=pred,aes(x=x,y=predicted,color=group))+geom_line()
b$fit=NA
b$fit[b$type=="adapt"]=predict(model2)
b$fit[b$type=="plast"]=predict(model)

b2=b %>% dplyr::group_by(maille,ctr_lon,ctr_lat,type) %>%
dplyr::summarise(contrib=mean(contrib),nsp=n_distinct(Species))



ggplot(data=b2,aes(x=ctr_lon,y=ctr_lat,col=contrib))+geom_point(aes(size=nsp))+
scale_color_gradientn(colors=rainbow(5))+
facet_wrap(~type)

ggplot(data=b,aes(x=ctr_lon,y=ctr_lat,col=contrib))+geom_point(aes(size=nsp))+scale_color_gradientn(rainbow(5))

model=lmer(contrib~(mfd+ctr_lat)*poly(date_fin,2)(1|Species),data=bidon2)

liste=unique(liste$Species)
for(i in 1:length(liste)){
bidon=subset(bidon2,Species==liste[i])
model=lmer(contrib~ctr_lat+ctr_lon,data=bidon2)



#### Viability ####

bidon=dcast(bmaster,Species+ORDRE+FAMILLE+mfd+trend.x+Std.Error+Mean_growth_rate+resul+
Precision+longmean2+latmean2+trend_err+trend_ang+abond~type,value.var="contrib")
bidon$Mean_growth_rate=as.numeric(bidon$Mean_growth_rate)
bidon$plast2=abs(bidon$plast)
bidon=subset(bidon,!is.na(Precision))
bidon$resid=residuals(lm(Mean_growth_rate~trend.x,data=bidon))
cor(bidon[,c("Mean_growth_rate","trend_ang","trend.x")])
bidon$mfd2="B"
bidon$mfd2[bidon$mfd<=180]="A"
bidon$mfd2[bidon$mfd>=220]="C"
model2=lmer(trend_ang~mfd*resul+latmean2+longmean2+(1|ORDRE/FAMILLE),data=bidon,
weights=sqrt(1/trend_err))

model=lmer(Mean_growth_rate~adapt+(plast)+latmean2+longmean2+(1|FAMILLE),data=subset(bidon,mfd>=200),
weights=sqrt(Precision))

b=ggpredict(model2,c("plast","mfd"))
ggplot(b,aes(x=x,y=predicted,color=group,group=group))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_line()

model=lmer(trend.x~adapt+mfd+plast+I(mfd^2)+longmean2+latmean2+(1|ORDRE/FAMILLE),data=bidon,
weights=sqrt(1/Std.Error))

b=ggeffect(model,c("plast"))
pl2=ggplot(b,aes(x=x,y=predicted))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_point(data=bidon,aes(x=plast,y=Mean_growth_rate,col=ORDRE),alpha=0.4)+geom_line(linetype="dashed")+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="none",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("c")+xlab("Contribution of plasticity\nto MFD shifts (day/year)")+
ylab("Mean growth rate")+labs(color="")+
scale_color_brewer(palette="Set1")

b=ggeffect(model,c("adapt"))
pl3=ggplot(b,aes(x=x,y=predicted))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_point(data=bidon,aes(x=adapt,y=Mean_growth_rate,col=ORDRE),alpha=0.4)+geom_line(linetype="dashed")+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="none",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("c")+xlab("Plasticity independent temporal trend in\nMFD shifts (day/year)")+
ylab("Mean growth rate")+labs(color="")+
scale_color_brewer(palette="Set1")

b=ggeffect(model,c("mfd"))
pl4=ggplot(b,aes(x=x,y=predicted))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_point(data=bidon,aes(x=mfd,y=Mean_growth_rate,col=ORDRE),alpha=0.4)+geom_line()+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("d")+xlab("Average mean flight date\n")+
ylab("Mean growth rate")+labs(color="")+
scale_color_brewer(palette="Set1")

b=ggeffect(model,c("latmean2"))
b$x=b$x+mean(resf$latitude)
pl5=ggplot(b,aes(x=x,y=predicted))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_point(data=bidon,aes(x=latmean2+mean(resf$latitude),y=Mean_growth_rate,col=ORDRE),alpha=0.4)+geom_line()+
theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("e")+xlab("Average latitude of records\n")+
ylab("Mean growth rate")+labs(color="")+
scale_color_brewer(palette="Set1")

grid.arrange(pl1,pl2,pl3,pl4,pl5,widths=c(1,1.4),heights=c(1,1,1),layout_matrix=matrix(c(1,2,3,1,4,5),nrow=3))

ggplot(data=subset(b,mfd>=200),aes(x=contrib_perc,y=trend.x,col=mfd))+geom_point()+facet_wrap(~type)


#PLOT PHYLOGENY
setwd(dir="C:/Users/Francois/Documents/papier 1 - données collection/data")
library(ape)
library(picante)
library(doBy)
library(ggtree)
library(phytools)
library(phylotools)
library(phangorn)
a=read.nexus("phylo_poly.nxs")
resf$FAMILLE=as.character(resf$FAMILLE)
resf$FAMILLE[which(resf$Espèce=="Tachina fera")]="Tachinidae"
resf$FAMILLE[which(resf$Espèce=="Bibio marci")]="Bibionidae"
resf$Species=gsub(" ","_",resf$Species,fixed=T)
library(caper)
a$node.label=1:length(a$node.label)
# lili=data.frame(ORDRE=unique(resf$ORDRE),lambda=NA,pval=NA,lambdavar=NA,pvalvar=NA)
# for(i in 1:nrow(lili)){
# resf2=subset(resf,ORDRE==lili$ORDRE[i])
# p3d=comparative.data(a,resf2[,c("Year_effect","diff_var","Espècegen","ORDRE")],Espècegen, vcv=TRUE)
# svl<-setNames(p3d$data[,c("Year_effect")],rownames(p3d$data))
# svl2<-setNames(p3d$data[,c("diff_var")],rownames(p3d$data))
# library(phylobase)
# library(phylosignal)
# p4d=phylo4d(as.phylo(p3d$phy),tip.data=svl)
# obj=phyloSignal(p4d, methods = c("Lambda"), reps=999, W = NULL)
# lili$lambda[i]=obj$stat[["Lambda"]]
# lili$pval[i]=obj$pvalue[["Lambda"]]
# p4d=phylo4d(as.phylo(p3d$phy),tip.data=svl2)
# obj=phyloSignal(p4d, methods = c("Lambda"), reps=999, W = NULL)
# lili$lambdavar[i]=obj$stat[["Lambda"]]
# lili$pvalvar[i]=obj$pvalue[["Lambda"]]
# }
a$edge.length[which(a$edge.length==0)]<-1e-7
a=force.ultrametric(a,method="nnls")
resf$contrib=NA
resf[,contrib:=Estimate*trend.y]
resf[type=="adapt",contrib:=Estimate]
b=subset(resf,!is.na(type)) %>% dplyr::group_by(Species,ORDRE,FAMILLE,type,trend.x,mfd,Std.Error,Mean_growth_rate,
Precision,longmean2,latmean2) %>%
dplyr::summarise(contrib=sum(contrib))
b=b %>% dplyr::group_by(Species) %>% dplyr::mutate(contrib_tot=sum(abs(contrib)),resul=sum(contrib))
b$contrib_perc=b$contrib/b$contrib_tot
b=b[order(b$mfd),]
b$Species=factor(b$Species,unique(b$Species))

bidon=dcast(b,Species+ORDRE+FAMILLE+mfd+trend.x+Std.Error+Mean_growth_rate+resul+
Precision+longmean2+latmean2~type,value.var="contrib_perc")
bidon$Mean_growth_rate=as.numeric(bidon$Mean_growth_rate)
bidon$plast2=abs(bidon$plast)
bidon=subset(bidon,ORDRE!="Hymenoptera")
p3d=comparative.data(a,bidon[,c("adapt","resul","mfd","plast","latmean2","longmean2","Species","ORDRE")],Species,
vcv=TRUE,na.omit=F)

model=pgls(adapt~1,data=p3d,lambda='ML')
p=ggtree(p3d$phy,layout="fan")

df=as.data.frame(p3d$data)
df$id=rownames(p3d$data)
df=df[,c("id",names(df)[1:(ncol(df)-1)])]
p1 <- p %<+% df +geom_tippoint(aes(color=adapt),size=2)+scale_color_viridis()

pl1=ggplot(data=subset(resf,!is.na(varia)),aes(x=varia,y=Estimate,col=ORDRE))+geom_hline(yintercept=0)+geom_boxplot()+theme_bw()+
theme(plot.title=element_text(size=14,face="bold"),legend.position="right",legend.title=element_blank(),
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("a")+xlab("Temperature indices (days of the year)")+ylab("Estimate phenotypic plasticity (day/°C)")+
scale_color_brewer(palette="Set1")

resf=resf %>% dplyr::group_by(ORDRE) %>% dplyr::mutate(nsp=length(unique(Species))) %>% ungroup()
b=resf %>% dplyr::group_by(varia,ORDRE) %>% dplyr::summarise(perc=length(unique(Species)),count=unique(nsp))
b=rbind(b,data.frame(varia=as.factor("210-300"),ORDRE="Coleoptera",perc=0,
count=unique(b$count[b$ORDRE=="Coleoptera"])))

pl2=ggplot(data=subset(b,!is.na(varia)),aes(x=varia,y=perc/count,col=NA,fill=ORDRE,group=ORDRE))+
geom_bar(stat="identity",position = "dodge2",width=0.5) + 
scale_y_continuous(labels=scales::percent)+theme_bw()+
theme(plot.title=element_text(size=14,face="bold"),legend.position="right",legend.title=element_blank(),
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("b")+xlab("Temperature indices (days of the year)")+ylab("Percentage of species responding")+
scale_color_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")

grid.arrange(pl1,pl2,ncol=1,heights=c(1.5,1))

png(paste0("fig_plast_season.png"),width=1400,height=1200,res=160)
grid.arrange(pl1,pl2,ncol=1,heights=c(1.5,1))
dev.off();







ggplot(data=subset(resf,!is.na(saison)),aes(x=saison,y=Estimate,col=ORDRE))+geom_boxplot()

resf$trend.y[is.na(resf$trend.y)]=1
b=subset(resf,!is.na(type)) %>% dplyr::group_by(Speciesgen,Species,ORDRE,type) %>%
dplyr::summarise(plast_contrib=sum(Estimate*trend.y,na.rm=T))

ggplot(data=b,aes(x=type,y=plast_contrib,col=ORDRE))+geom_boxplot()

b=subset(resf,!is.na(saison)) %>% dplyr::group_by(Speciesgen,Species,ORDRE,saison) %>%
dplyr::summarise(plast_contrib=sum(Estimate*trend.y,na.rm=T))
ggplot(data=b,aes(x=saison,y=plast_contrib,col=ORDRE))+geom_boxplot()


ggplot(data=resf,aes(x=trend.x,y=Estimate))+geom_point()+facet_wrap(~varia,scales="free")


library(raster)
library(sf)
library(rgdal)
setwd(dir="C:/Users/Francois/Documents/land use change/European_costlines")
shp <- readOGR(".", "Europe_coastline")
shp=spTransform(shp,CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 
+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
AG=fortify(shp)


vec=unique(resqualf$variable)
vec=vec[grep("temp",vec)]
vec=vec[grep("trend",vec,invert=T)]
bidon1=resqualf[resqualf$variable %in% vec,]
bidon2=resqualf[resqualf$variable %in% paste0("trend_",vec),]
names(bidon2)[names(bidon2)=="value"]="trend_temp"
bidon2$variable=gsub("trend_","",bidon2$variable)
bidon=merge(bidon1,bidon2[,c("Speciesgen","maille","trend_temp","variable")],
by=c("Speciesgen","maille","variable"))
bidon1=resqualf[resqualf$variable=="Annee",]
bidon1$trend_temp=NA
bidon=rbind(bidon,bidon1)
bidon$contrib=NA
bidon$contrib=bidon$value*bidon$trend_temp
bidon$contrib[bidon$varia=="Annee"]=bidon$value[bidon$varia=="Annee"]
bidon$type=NA
bidon$type[grep("temp",bidon$varia)]="plast"
bidon$type[grep("Annee",bidon$varia)]="Annee"
b=subset(bidon,n>50) %>% dplyr::group_by(Speciesgen,type,maille,ctr_lon,ctr_lat) %>%
dplyr::summarise(contrib=sum(contrib))
b2=b %>% dplyr::group_by(type,maille,ctr_lon,ctr_lat) %>% dplyr::summarise(nb=length(contrib),
contrib=mean(contrib))
AG2=AG %>% dplyr::filter(long>=min(b2$ctr_lon) & long<=max(b2$ctr_lon) &
lat>=min(b2$ctr_lat) & lat<=max(b2$ctr_lat))

pl1=ggplot()+geom_raster(data=subset(b2,type=="plast" & nb>=10),aes(x=ctr_lon,y=ctr_lat,fill=contrib))+
scale_fill_viridis()+coord_fixed(ratio = 1)+ggtitle("a - plasticity")+
geom_path(data=AG2,aes(x=long,y=lat,group=group),color="black")+labs(fill="days/year")

pl2=ggplot()+geom_raster(data=subset(b2,type=="Annee" & nb>=10),aes(x=ctr_lon,y=ctr_lat,fill=contrib))+
scale_fill_viridis()+coord_fixed(ratio = 1)+ggtitle("b - adaptation")+
geom_path(data=AG2,aes(x=long,y=lat,group=group),color="black")+labs(fill="days/year")
grid.arrange(pl1,pl2,ncol=2)
plot(value~trend_temp,data=subset(bidon,n>100))



resf$contrib=NA
resf$contrib[!is.na(resf$trend)]=resf$Estimate[!is.na(resf$trend)]*resf$trend[!is.na(resf$trend)]
resf$contrib[resf$varia=="Annee"]=resf$Estimate[resf$varia=="Annee"]
resf$type=NA
resf$type[grep("temp",resf$varia)]="plast"
resf$type[grep("Annee",resf$varia)]="Annee"
b=resf %>% dplyr::group_by(Speciesgen,type) %>% dplyr::summarise(contrib=sum(contrib))
ggplot(data=b,aes(y=contrib,x=type))+geom_boxplot()


resf=as.data.frame(resf)
resf$type="autre"
resf$type[resf$varia=="Annee"]="adapt"
resf$type[grep("temp",resf$varia)]="plast"
resf$contrib=NA
resf$contrib[resf$type=="plast"]=resf$Estimate[resf$type=="plast"]*resf$trend[resf$type=="plast"]
resf$contrib[resf$type=="adapt"]=resf$Estimate[resf$varia=="Annee"]
resf=resf %>% dplyr::group_by(Speciesgen,type) %>%
dplyr::mutate(contrib_moy=sum(contrib))

ggplot(data=resf,aes(x=type,y=contrib))+geom_boxplot()
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data_occ_UK")
occt=fread("trends.txt",sep="\t",header=T)
resf=merge(resf,liste,by="Speciesgen")
resf=merge(resf,occt,by="Species")

ggplot(data=subset(unique(resf[,c("plast_moy","adapt","Mean_growth_rate")]),adapt>-2),
aes(x=adapt,y=plast_moy,col=as.numeric(Mean_growth_rate)))+geom_point()+
stat_smooth(method="lm")+scale_color_gradient2()


form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+Annee+(1|Annee)+(1|maille)"))
model=fitme(form,data=bidon)

fwrite(liste,"temperature_select_var.txt",sep="\t")
