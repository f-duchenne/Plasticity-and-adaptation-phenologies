library(doBy)
library(geosphere)
library(dplyr)
library(fitdistrplus)
library(data.table)
library(lme4)
library(viridis)
library(MuMIn)
library(sf)
library(ggeffects)
library(RColorBrewer)
library(ggraph)
library(igraph)
library("ggExtra")
library(scales)
library(stringr)
library(metR)
library(cowplot)
library(car)
library(DHARMa)
library(glmmTMB)
library(dendextend)
library(ape)
library(ggtree)
library(ggnewscale)
library(geiger)
library(ggthemes)
library(ggplotify)
library(ggtern)
library(ggpubr)
library(ggh4x)

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
liste=fread("selec_temp_var_beta.txt",sep="\t")
liste=subset(liste,nb_pre1990>=50 & nb_post1990>=500)
vec=names(liste)[grep("temp",names(liste),fixed=T)]
nb=liste[,c("Speciesgen","Species","FAMILLE","ORDRE","nb_pre1990","nb_post1990")]

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
resf=NULL
lili=list.files()
lili=lili[grep("_qual.txt",lili,invert=T)]
lili=lili[grep("RData",lili,invert=T)]
bb=0
for(i in lili){
bb=bb+1
bidon=fread(i,sep="\t")
#if(bb!=1){resf
resf=rbind(resf,bidon)
}

resf=resf[,which(!duplicated(names(resf))),with=F]
resf=merge(resf,unique(liste[,c("Speciesgen","ORDRE","FAMILLE","Species")]),by="Speciesgen")

max(abs(resf$cor),na.rm=T)
min(resf$nb_annee)
min(resf$nbdata)

resf$var_std=resf$STD_maille/abs(resf$Estimate)

resf$varia2=as.character(resf$varia)
resf$varia2[resf$varia2=="Annee2"]="Year effect"
resf$varia2[grep("temp_",resf$varia2,fixed=T)]="Phenotypic plasticity"
resf=resf[,which(!duplicated(names(resf))),with=F]
resf$type=NA
resf$type[grep("temp",resf$varia)]="plast"
resf$type[grep("Annee",resf$varia)]="adapt"
names(resf)[names(resf)=="Std. Error"]="sde"
resf$lwr=resf$Estimate-1.96*resf$sde
resf$upr=resf$Estimate+1.96*resf$sde
resf$Est.signi=">0.05"
resf$Est.signi[resf$lwr>0]="<0.05"
resf$Est.signi[resf$upr<0]="<0.05"
resf$Est.signi=factor(resf$Est.signi,c(">0.05","<0.05"))

coucou=c("magenta2",brewer.pal(5,"Set2")[5],"red","dodgerblue4")

boxplot(var_std~varia2,data=subset(resf,abs(var_std)<10))

resf %>% group_by(ORDRE) %>% summarise(nb=length(unique(Speciesgen)))

######### FIGURE 3
#calculate contributions in day/year
resf[,contrib:=Estimate*trend]
resf[type=="adapt",contrib:=Estimate]

resf[,contrib_lwr:=apply(cbind(lwr*trend,upr*trend),1,min)]
resf[type=="adapt",contrib_lwr:=lwr]

resf[,contrib_upr:=apply(cbind(lwr*trend,upr*trend),1,max)]
resf[type=="adapt",contrib_upr:=upr]

contrib_tab=subset(resf,type %in% c("plast","adapt")) %>% dplyr::group_by(Speciesgen,ORDRE,FAMILLE,Species,mfd) %>% dplyr::summarise(contrib_plast=sum(contrib[type=="plast"],na.rm=T),contrib_plast_lwr=sum(contrib_lwr[type=="plast"],na.rm=T),
contrib_plast_upr=sum(contrib_upr[type=="plast"],na.rm=T),
contrib_adapt=sum(contrib[type=="adapt"],na.rm=T),contrib_adapt_lwr=sum(contrib_lwr[type=="adapt"],na.rm=T),contrib_adapt_upr=sum(contrib_upr[type=="adapt"],na.rm=T))

contrib_tab$Est.signi_plast=">0.05"
contrib_tab$Est.signi_plast[contrib_tab$contrib_plast_lwr>0]="<0.05"
contrib_tab$Est.signi_plast[contrib_tab$contrib_plast_upr<0]="<0.05"
contrib_tab$Est.signi_adapt=">0.05"
contrib_tab$Est.signi_adapt[contrib_tab$contrib_adapt_lwr>0]="<0.05"
contrib_tab$Est.signi_adapt[contrib_tab$contrib_adapt_upr<0]="<0.05"
contrib_tab$Est.signi_plast=factor(contrib_tab$Est.signi_plast,c(">0.05","<0.05"))

contrib_tab$contrib_adapt_lwr[contrib_tab$contrib_adapt_lwr<(min(contrib_tab$contrib_adapt)-0.05)]=min(contrib_tab$contrib_adapt)-0.05
contrib_tab$contrib_adapt_upr[contrib_tab$contrib_adapt_upr>(max(contrib_tab$contrib_adapt)+0.05)]=max(contrib_tab$contrib_adapt)+0.05
contrib_tab$contrib_plast_lwr[contrib_tab$contrib_plast_lwr<(min(contrib_tab$contrib_plast)-0.05)]=min(contrib_tab$contrib_plast)-0.05
contrib_tab$contrib_plast_upr[contrib_tab$contrib_plast_upr>(max(contrib_tab$contrib_plast)+0.05)]=max(contrib_tab$contrib_plast)+0.05

var(contrib_tab$contrib_plast)
var(contrib_tab$contrib_adapt)

#plast numbers
nrow(subset(contrib_tab,contrib_plast<0 & Est.signi_plast=="<0.05"))/nrow(contrib_tab)
nrow(subset(contrib_tab,contrib_plast>0 & Est.signi_plast=="<0.05"))/nrow(contrib_tab)
#adapt numbers
nrow(subset(contrib_tab,contrib_adapt<0 & Est.signi_adapt=="<0.05"))/nrow(contrib_tab)
nrow(subset(contrib_tab,contrib_adapt>0 & Est.signi_adapt=="<0.05"))/nrow(contrib_tab)

pl1=ggplot(data=contrib_tab,aes(x=contrib_plast,fill=Est.signi_plast))+geom_histogram(col="black")+
geom_vline(xintercept=0,col="red")+
scale_fill_manual(values=c("white","black"))+theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
legend.title = element_blank(),legend.position="none",plot.title=element_text(size=14,face="bold"),axis.title=element_blank())+
xlab("")+ylab("Number of species")+ggtitle("a",subtitle="Phenotypic plasticity")+coord_cartesian(expand=F)

pl2=ggplot(data=subset(resf,type=="adapt"),aes(x=Estimate,fill=Est.signi))+geom_histogram(col="black")+
geom_vline(xintercept=0,col="red")+
scale_fill_manual(values=c("white","black"))+theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
legend.title = element_blank(),legend.position="none",plot.title=element_text(size=14,face="bold"),axis.title=element_blank())+
xlab("")+ylab("Number of species")+ggtitle("b",subtitle="Evolution")+coord_cartesian(expand=F)


obj=expand.grid(seq(min(contrib_tab$contrib_adapt)-0.05,max(contrib_tab$contrib_adapt)+0.05,0.005),seq(min(contrib_tab$contrib_plast)-0.05,max(contrib_tab$contrib_plast)+0.05,0.005))
obj$value=obj[,1]+obj[,2]

contrib_tab$gradient="co-gradient"
contrib_tab$gradient[contrib_tab$Est.signi_adapt=="<0.05" & contrib_tab$Est.signi_plast=="<0.05" &
sign(contrib_tab$contrib_plast)!=sign(contrib_tab$contrib_adapt)]="counter-gradient"
contrib_tab$gradient[contrib_tab$Est.signi_adapt=="<0.05" & contrib_tab$Est.signi_plast==">0.05"]="zero or one mechanism significant only"
contrib_tab$gradient[contrib_tab$Est.signi_adapt==">0.05" & contrib_tab$Est.signi_plast=="<0.05"]="zero or one mechanism significant only"
contrib_tab$gradient[contrib_tab$Est.signi_adapt==">0.05" & contrib_tab$Est.signi_plast==">0.05"]="zero or one mechanism significant only"

pl=ggplot(data=contrib_tab,aes(y=contrib_plast,x=contrib_adapt))+
geom_raster(data=obj,aes(x=Var1,y=Var2,fill=value))+
scale_fill_gradient2(name="Resulting\npheno. shift\n(day/year)")+
geom_vline(xintercept=0,col="black")+geom_hline(yintercept=0,col="black")+geom_abline(intercept=0,slope=-1,linetype="dashed")+
geom_errorbarh(aes(xmin=contrib_adapt_lwr,xmax=contrib_adapt_upr),alpha=0.1)+
geom_errorbar(aes(ymin=contrib_plast_lwr,ymax=contrib_plast_upr),alpha=0.1)+
geom_point(col="black",alpha=0.7)+
theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"))+
ylab("Contribution of phenotypic plasticity to phenological shifts (day/year)")+
xlab("Contribution of evolution to phenological shifts (day/year)")+coord_fixed(ratio=1,expand=F)+ggtitle("c")+
stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,label.x=0.5,label.y=0.35)+
scale_fill_manual(values=c("white","black"))


haut=grid.arrange(pl1,pl2,ncol=2,bottom="Contribution to phenological shifts (day/year)",left="Number of species")

layout=rbind(c(1,1,1,1),c(3,3,3,4))
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
pdf("figure3.pdf",height=7,width=10)
grid.arrange(haut,pl,ncol=1,heights=c(1,2.6))
dev.off();

############## FIGURE 4
#subset(resf,Est.signi==">0.05" & varia2=="Year effect")[,c("Speciesgen","FAMILLE","ORDRE","nbdata")]
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
temp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")

dat=fread("selec_temp_var_beta.txt",sep="\t")
liste=data.frame(Speciesgen=unique(dat$Speciesgen),nb_pre1990=dat$nb_pre1990[!duplicated(dat$Speciesgen)],nb_post1990=dat$nb_post1990[!duplicated(dat$Speciesgen)])
liste$nb=liste$nb_pre1990+liste$nb_post1990

years_to_pred=c(1960,1990,2010)

#### EXAMPLE 1
sp="Abraxas sylvata_NA"
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
bidon=fread(paste0("data_per_species/",sp,".txt"))
bidon=subset(bidon, !is.na(elev) & !is.na(temp_0_90))
varia=subset(dat,Speciesgen==sp)$varia
moy=mean(bidon$Annee)
vec=bidon$Annee
bidon$Annee2=bidon$Annee-mean(bidon$Annee)
bidon$Annee=bidon$Annee-1960
ann_pre=years_to_pred-moy

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
load(paste0("model_",gsub(" ","_",sp,fixed=T),".RData"))

b=ggpredict(model,c("temp_75_165",paste0("Annee2[",paste(ann_pre,collapse=","),"]")))
b$group2=as.factor(years_to_pred)

newdat=bidon %>% dplyr::group_by(plyr::round_any(Annee2,5)) %>%
dplyr::summarise(temp_75_165=mean(temp_75_165),temp_285_315=mean(bidon$temp_285_315),Altitude=mean(bidon$Altitude),maille=NA,Annee=NA)
names(newdat)[1]="Annee2"
newdat=subset(newdat,Annee2 %in% plyr::round_any(c(years_to_pred-moy),5))
newdat$group2=as.factor(years_to_pred)
newdat$x=newdat$temp_75_165
newdat$predicted=predict(model,newdata=newdat,re.form=NA)
newdat2=newdat
newdat2$Annee2=newdat2$Annee2[1]
newdat2$predicted=predict(model,newdata=newdat2,re.form=NA)

newdat2$scenario="without evolution"
newdat$scenario="with evolution"
newdat=rbind(newdat,newdat2[-1,])

pl1=ggplot()+
geom_jitter(data=bidon,aes(x=temp_75_165,y=Jour.de.collecte,color=Annee+1960),height=0,width=0.05)+
scale_color_gradient(low=brewer.pal(8, "Blues")[c(2)],high=brewer.pal(9, "Blues")[9],name="Year of observation",breaks=c(1960,1980,2000))+
new_scale_colour()+
geom_ribbon(data=b,aes(x=x,y=predicted,color=group2,fill=group2,ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_line(data=b,aes(x=x,y=predicted,color=group2,fill=group2),size=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right")+xlab("Temperature anomalies (°C)")+coord_cartesian(expand=F)+
scale_y_continuous(name = "Predicted mean flight date (day of the year)")+
scale_color_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)],name="Predicted\nreaction norm")+scale_fill_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)],name="Predicted\nreaction norm")+
labs(col="",fill="")+ggtitle("a")+
geom_point(data=newdat,size=4,color="black",stroke = 2,aes(x=x,y=predicted,color=group2,fill=group2,shape=scenario))+
scale_shape_manual(values=c(21,24),name="Estimated position\non reaction norm")+labs(shape="")+
 guides(fill = guide_legend(order = 3,override.aes = list(linetype = 1)),color = guide_legend(order = 3,override.aes = list(linetype = 1,shape=NA)),
 shape=guide_legend(order=4))
##### EXAMPLE 2
sp="Andrena fulva_NA"
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
bidon=fread(paste0("data_per_species/",sp,".txt"))
bidon=subset(bidon, !is.na(elev) & !is.na(temp_0_90))
varia=subset(dat,Speciesgen==sp)$varia
moy=mean(bidon$Annee)
vec=bidon$Annee
bidon$Annee2=bidon$Annee-mean(bidon$Annee)
bidon$Annee=bidon$Annee-1960
ann_pre=sort(unique(bidon$Annee2[vec %in%years_to_pred]))

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
load(paste0("model_",gsub(" ","_",sp,fixed=T),".RData"))

b=ggpredict(model,c("temp_30_120",paste0("Annee2[",paste(ann_pre,collapse=","),"]")))
b$group2=as.factor(years_to_pred)

newdat=bidon %>% dplyr::group_by(plyr::round_any(Annee2,5)) %>%
dplyr::summarise(temp_30_120=mean(temp_30_120),temp_165_195=mean(bidon$temp_165_195),Altitude=mean(bidon$Altitude),maille=NA,Annee=NA)
names(newdat)[1]="Annee2"
newdat=subset(newdat,Annee2 %in% plyr::round_any(c(years_to_pred-moy),5))
newdat$group2=as.factor(years_to_pred)
newdat$x=newdat$temp_30_120
newdat$predicted=predict(model,newdata=newdat,re.form=NA)
newdat2=newdat
newdat2$Annee2=newdat2$Annee2[1]
newdat2$predicted=predict(model,newdata=newdat2,re.form=NA)

newdat2$scenario="without evolution"
newdat$scenario="with evolution"
newdat=rbind(newdat,newdat2[-1,])

pl2=ggplot()+
geom_jitter(data=bidon,aes(x=temp_30_120,y=Jour.de.collecte,color=Annee+1960),height=0,width=0.05)+
scale_color_gradient(low=brewer.pal(8, "Blues")[c(2)],high=brewer.pal(9, "Blues")[9],name="Year of observation",breaks=c(1960,1980,2000))+
new_scale_colour()+
geom_ribbon(data=b,aes(x=x,y=predicted,color=group2,fill=group2,ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_line(data=b,aes(x=x,y=predicted,color=group2,fill=group2),size=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right")+xlab("Temperature anomalies (°C)")+coord_cartesian(expand=F)+
scale_y_continuous(name = "Predicted mean flight date (day of the year)")+
scale_color_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)])+scale_fill_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)])+
labs(col="",fill="")+ggtitle("c")+
geom_point(data=newdat,size=4,color="black",stroke = 2,aes(x=x,y=predicted,color=group2,fill=group2,shape=scenario))+
scale_shape_manual(values=c(21,24),name="Estimated position\non reaction norm")+labs(shape="")+
 guides(fill = guide_legend(order = 3,override.aes = list(linetype = 1)),color = guide_legend(order = 3,override.aes = list(linetype = 1,shape=NA)),
 shape=guide_legend(order=4))
 
##### EXAMPLE 3
sp="Bombus vestalis_NA"
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
bidon=fread(paste0("data_per_species/",sp,".txt"))
bidon=subset(bidon, !is.na(elev) & !is.na(temp_0_90))
varia=subset(dat,Speciesgen==sp)$varia
moy=mean(bidon$Annee)
vec=bidon$Annee
bidon$Annee2=bidon$Annee-mean(bidon$Annee)
bidon$Annee=bidon$Annee-1960
ann_pre=years_to_pred-moy

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
load(paste0("model_",gsub(" ","_",sp,fixed=T),".RData"))

b=ggpredict(model,c("temp_30_120",paste0("Annee2[",paste(ann_pre,collapse=","),"]")))
b$group2=as.factor(years_to_pred)

newdat=bidon %>% dplyr::group_by(plyr::round_any(Annee2,5)) %>%
dplyr::summarise(temp_30_120=mean(temp_30_120),temp_315_345=mean(bidon$temp_315_345),Altitude=mean(bidon$Altitude),maille=NA,Annee=NA)
names(newdat)[1]="Annee2"
newdat=subset(newdat,Annee2 %in% plyr::round_any(c(years_to_pred-moy),5))
newdat$group2=as.factor(years_to_pred)
newdat$x=newdat$temp_30_120
newdat$predicted=predict(model,newdata=newdat,re.form=NA)
newdat2=newdat
newdat2$Annee2=newdat2$Annee2[1]
newdat2$predicted=predict(model,newdata=newdat2,re.form=NA)

newdat2$scenario="without evolution"
newdat$scenario="with evolution"
newdat=rbind(newdat,newdat2[-1,])

pl3=ggplot()+
geom_jitter(data=bidon,aes(x=temp_30_120,y=Jour.de.collecte,color=Annee+1960),height=0,width=0.02)+
scale_color_gradient(low=brewer.pal(8, "Blues")[c(2)],high=brewer.pal(9, "Blues")[9],name="Year of\nobservations",breaks=c(1960,1980,2000))+
new_scale_colour()+
geom_ribbon(data=b,aes(x=x,y=predicted,color=group2,fill=group2,ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA)+
geom_line(data=b,aes(x=x,y=predicted,color=group2,fill=group2),size=1)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right")+xlab("Temperature anomalies (°C)")+coord_cartesian(expand=F)+
scale_y_continuous(name = "Predicted mean flight date (day of the year)")+
scale_color_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)],name="Predicted\nreaction norm")+
scale_fill_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)],name="Predicted\nreaction norm")+
labs(col="",fill="")+ggtitle("b")+
geom_point(data=newdat,size=4,color="black",stroke = 2,aes(x=x,y=predicted,color=group2,fill=group2,shape=scenario))+
scale_shape_manual(values=c(21,24),name="Estimated position\non reaction norm")+labs(shape="")+
 guides(fill = guide_legend(order = 3,override.aes = list(linetype = 1)),color = guide_legend(order = 3,override.aes = list(linetype = 1,shape=NA)),
 shape=guide_legend(order=4))
 
leg <- ggpubr::as_ggplot(cowplot::get_legend(pl3))
pl1=pl1+theme(legend.position="none")
pl2=pl2+theme(legend.position="none")
pl3=pl3+theme(legend.position="none")

plot_grid(pl1,pl2,pl3,leg,align="hv",ncol=4,rel_widths=c(1,1,1,0.4))

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
pdf("figure4.pdf",height=4,width=12)
plot_grid(pl1,pl3,pl2,leg,align="hv",ncol=4,rel_widths=c(1,1,1,0.5))
dev.off();

###########################################################

######### FIGURE S8 & S9
setwd(dir="D:/land use change/European_costlines")
shp <- st_read(".", "Europe_coastline")
b=subset(resf,type=="plast")
bidon=as.data.frame(str_split_fixed(b$varia, "_", 3)[,2:3])
bidon=as.data.frame(sapply(bidon,as.numeric))
names(bidon)=c("debut","fin")
bidon$day=bidon$debut
bidon$day[bidon$debut>bidon$fin]=bidon$debut[bidon$debut>bidon$fin]-365
bidon$duration=NA
bidon$duration=bidon$fin-bidon$debut
bidon$duration[bidon$debut>bidon$fin]=bidon$fin[bidon$debut>bidon$fin]+365-bidon$debut[bidon$debut>bidon$fin]

b=cbind(b,bidon)
sfpts=sf::st_as_sf(b, coords = c("longmean","latmean"),crs=st_crs(shp))  %>% st_transform(st_crs(4326))
b=cbind(b,st_coordinates(sfpts))

model=lmer(Estimate~poly(day,3)*as.factor(duration)+mfd+X+Y+(1|ORDRE/FAMILLE),data=b,weights=1/sde)
Anova(model)
summary(model)

obj=ggpredict(model,c("day[all]","duration"))
obj[obj$group==90 & obj$x>275,c("predicted","conf.high","conf.low")]=NA
obj[obj$group==30 & obj$x<(-30),c("predicted","conf.high","conf.low")]=NA

obj[which.min(obj$predicted),]

pl1=ggplot()+geom_point(data=b,aes(x=day,y=Estimate,color=as.factor(duration)),alpha=0.4)+
xlab("Starting day of the temperature index")+ylab("Phenotypic plasticity")+
geom_ribbon(data=obj,aes(x=x,y=predicted,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+
geom_line(data=obj,aes(x=x,y=predicted,color=group),size=1.2)+labs(color="Duration",fill="Duration")+theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),legend.position="none",
plot.title=element_text(size=14,face="bold"))+ggtitle("a")+scale_color_manual(values=c("#3581B8","#FCB07E"))+scale_fill_manual(values=c("#3581B8","#FCB07E"))

obj=ggpredict(model,c("mfd[all]","duration"))

pl2=ggplot()+geom_point(data=b,aes(x=mfd,y=Estimate,color=as.factor(duration)),alpha=0.4)+
xlab("Mean flight date")+ylab("Phenotypic plasticity")+
geom_ribbon(data=obj,aes(x=x,y=predicted,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+
geom_line(data=obj,aes(x=x,y=predicted,color=group),size=1.2)+labs(color="Duration",fill="Duration")+theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),legend.position="none",
plot.title=element_text(size=14,face="bold"))+ggtitle("b")+scale_color_manual(values=c("#3581B8","#FCB07E"))+scale_fill_manual(values=c("#3581B8","#FCB07E"))

obj=ggpredict(model,c("Y[all]","duration"))

pl3=ggplot()+geom_point(data=b,aes(x=Y,y=Estimate,color=as.factor(duration)),alpha=0.4)+
ylab("Phenotypic plasticity")+
labs(color="Duration",fill="Duration")+theme_bw()+
geom_ribbon(data=obj,aes(x=x,y=predicted,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+
geom_line(data=obj,aes(x=x,y=predicted,color=group),size=1.2,linetype="dashed")+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),legend.position="none",
plot.title=element_text(size=14,face="bold"))+ggtitle("c")+scale_color_manual(values=c("#3581B8","#FCB07E"))+
scale_x_latitude()+scale_fill_manual(values=c("#3581B8","#FCB07E"))+scale_x_latitude(name="Average latitude",ticks =1)

obj=ggpredict(model,c("X[all]","duration"))

pl4=ggplot()+geom_point(data=b,aes(x=X,y=Estimate,color=as.factor(duration)),alpha=0.4)+
ylab("Phenotypic plasticity")+
labs(color="Duration",fill="Duration")+theme_bw()+
geom_ribbon(data=obj,aes(x=x,y=predicted,fill=group,ymin=conf.low,ymax=conf.high),alpha=0.2,color=NA)+
geom_line(data=obj,aes(x=x,y=predicted,color=group),size=1.2,linetype="dashed")+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
plot.title=element_text(size=14,face="bold"))+ggtitle("d")+scale_color_manual(values=c("#3581B8","#FCB07E"))+
scale_x_latitude()+scale_fill_manual(values=c("#3581B8","#FCB07E"))+scale_x_longitude(name="Average longitude",ticks =1)

grid.arrange(pl1,pl2,pl3,pl4,ncol=4,widths=c(1,1,1,1.2))
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
png("figure_S8.png",height=600,width=1600,res=120)
grid.arrange(pl1,pl2,pl3,pl4,ncol=4,widths=c(1,1,1,1.2))
dev.off();

obj=insight::get_variance(model)
b=data.frame(varia=c("Taxonomic (random effects)","Saptio-temporal (fixed effects)","Residual"),value=c(obj$var.random,obj$var.fixed,obj$var.residual))
contrib_tab$value=contrib_tab$value/sum(contrib_tab$value)
pl1=ggplot(data=b,aes(x=varia,y=value))+geom_bar(stat="identity",col="white")+theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
plot.title=element_text(size=14,face="bold"),axis.title.x=element_blank())+ylab("Percentage of variance in phenotypic plasticity explained")+
scale_y_continuous(labels=percent)+coord_cartesian(expand=F)

png("figure_S9.png",height=800,width=700,res=120)
pl1
dev.off();

################## FIGURE 5
setwd(dir="D:/land use change/European_costlines")
shp <- st_read(".", "Europe_coastline")
b=subset(resf,type=="adapt")

sfpts=sf::st_as_sf(b, coords = c("longmean","latmean"),crs=st_crs(shp))  %>% st_transform(st_crs(4326))
b=cbind(b,st_coordinates(sfpts))
model=lmer(Estimate~mfd+X+Y+(1|ORDRE/FAMILLE),data=b,weights=1/sde)
Anova(model)
summary(model)

pre=ggpredict(model,"mfd")

pl1=ggplot()+
geom_point(data=b,aes(x=mfd,y=Estimate))+geom_errorbar(data=b,aes(x=mfd,y=Estimate,ymax=contrib_upr,ymin=contrib_lwr),alpha=0.2)+
geom_ribbon(data=pre,aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high),alpha=0.2,col=NA,fill="lightgrey")+
geom_line(data=pre,aes(x=x,y=predicted),size=1,color="black")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),
plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="right")+xlab("Species mean flight date (day of the year)")+coord_cartesian(expand=F)+
scale_y_continuous(name = "Contribution of evolution to phenological shifts (day/year)")+
scale_color_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)])+scale_fill_manual(values=brewer.pal(8, "Blues")[c(3, 5, 8)])+
labs(col="",fill="")+ggtitle("a")

obj=insight::get_variance(model)
b=data.frame(varia=c("Taxonomic\n(random effects)","Saptio-temporal\n(fixed effects)","Residual"),value=c(obj$var.random,obj$var.fixed,obj$var.residual))
b$value=b$value/sum(b$value)
pl2=ggplot(data=b,aes(x=varia,y=value))+geom_bar(stat="identity",col="black")+theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
plot.title=element_text(size=14,face="bold"),axis.title.x=element_blank(),axis.text=element_text(angle=0,hjust=0.5))+ylab("Percentage of explained variance in evolution")+
scale_y_continuous(labels=percent)+coord_cartesian(expand=F)+ggtitle("b")


rr <- ranef(model)
dd <- as.data.frame(rr)
bb=transform(dd, lwr = condval - 1.96*condsd, upr = condval + 1.96*condsd)
bb=bb[order(bb$condval),]
bb=cbind(bb,as.data.frame(str_split_fixed(bb$grp, ":", 2)))
bb$Family=factor(bb$V1,levels=bb$V1)
pl3=ggplot(data=subset(bb,grpvar=="ORDRE"),aes(x=Family,y=condval,color=V1))+
geom_pointrange(aes(ymin=lwr,ymax=upr))+coord_flip()+geom_hline(yintercept=0,linetype="dashed")+ggtitle("c")+theme_bw()+
theme(strip.background = element_blank(),panel.border=element_blank(),panel.grid=element_blank(),axis.line=element_line(),
plot.title=element_text(size=14,face="bold"),legend.position="none")+ylab("Deviation from average in evolution")+scale_colour_colorblind()+xlab("Taxonomic order")


setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
pdf("figure5.pdf",height=5,width=12)
plot_grid(pl1,pl2,pl3,ncol=3,align="h")
dev.off();



pl3=ggplot(data=subset(bb,grpvar=="FAMILLE:ORDRE"),aes(x=Family,y=condval,color=V2))+
geom_pointrange(aes(ymin=lwr,ymax=upr))+coord_flip()+geom_hline(yintercept=0,linetype="dashed")+ggtitle("b")+theme_bw()+
theme(strip.background = element_blank(),panel.border=element_blank(),panel.grid=element_blank(),axis.line=element_line(),
plot.title=element_text(size=14,face="bold"))+ylab("Deviation from average in evolution")+scale_colour_colorblind()+labs(color="Order")


png("figure_S10.png",height=1000,width=700,res=120)
pl3
dev.off();

############FIGURE S6

plot(mean~Estimate,data=resf[resf$varia %in% c("Annee2"),])
abline(0,1)
cor(resf[(resf$varia %in% c("Annee2")),c("mean","Estimate")])

pl1=ggplot(data=resf[resf$varia %in% c("Annee2"),],aes(x=mean,y=Estimate))+
theme_bw()+geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
geom_point(alpha=0.55)+stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,label.x=-0.7)+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"),axis.line = element_line(colour = "black"))+
ylab("Year effect from LMER model (day/year)")+xlab("Year effect from INLA model (day/year)")+
coord_cartesian(expand=T)+ggtitle("a")

pl2=ggplot(data=resf[grep("temp_",resf$varia),],aes(x=mean,y=Estimate))+
theme_bw()+geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
geom_point(alpha=0.55)+stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,label.x=-0.7)+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"),axis.line = element_line(colour = "black"))+
ylab("Phenotypic plasticity estimates from LMER model (day/°C)")+xlab("Phenotypic plasticity estimates from INLA model (day/°C)")+
coord_cartesian(expand=T)+ggtitle("b")


setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
png("figure_S6.png",height=700,width=1200,res=120)
grid.arrange(pl1,pl2,ncol=2)
dev.off();

################## FIGURE S11
resf=resf %>% group_by(Speciesgen) %>% mutate(contrib_plast=sum(contrib[type=="plast"],na.rm=T),contrib_plast_lwr=sum(contrib_lwr[type=="plast"],na.rm=T),
contrib_plast_upr=sum(contrib_upr[type=="plast"],na.rm=T))
bb=subset(resf,type=="adapt")

pl1=ggplot(data=bb,aes(x=contrib_plast,y=longterm_mfd))+
theme_bw()+geom_vline(xintercept=0,col="black",linetype="dashed")+geom_hline(yintercept=0,col="black",linetype="dashed")+
geom_errorbarh(aes(xmin=contrib_plast_lwr,xmax=contrib_plast_upr),alpha=0.1)+
geom_errorbar(aes(ymin=longterm_mfd-1.96*longterm_mfd.se,ymax=longterm_mfd+1.96*longterm_mfd.se),alpha=0.1)+
geom_point(alpha=0.5)+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"),axis.line = element_line(colour = "black"))+
ylab("Long-term phenological shifts (day/year)")+xlab("Contribution of phenotypic plasticity(day/year)")+
coord_cartesian(expand=F)+ggtitle("a")+stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,label.x=-0.92)

pl2=ggplot(data=bb,aes(x=Estimate,y=longterm_mfd))+
theme_bw()+geom_vline(xintercept=0,col="black",linetype="dashed")+geom_hline(yintercept=0,col="black",linetype="dashed")+
geom_errorbarh(aes(xmin=Estimate-1.96*sde,xmax=Estimate+1.96*sde),alpha=0.1)+
geom_errorbar(aes(ymin=longterm_mfd-1.96*longterm_mfd.se,ymax=longterm_mfd+1.96*longterm_mfd.se),alpha=0.1)+
geom_point(alpha=0.55)+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"),axis.line = element_line(colour = "black"))+
ylab("Long-term phenological shifts (day/year)")+xlab("Contribution of evolution(day/year)")+
coord_cartesian(expand=F)+ggtitle("b")+stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,label.x=0.1,label.y=-1)


plot_grid(pl1,pl2,pl3,align="hv",ncol=3,rel_widths=c(1.2,1,1))

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
png("figure_S11.png",height=700,width=1200,res=120)
plot_grid(pl1,pl2,align="hv",ncol=2)
dev.off();