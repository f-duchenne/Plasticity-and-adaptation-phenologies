library(data.table)
library(dplyr)
library(raster)
library(Hmisc)
library(car)
library(doBy)
library(ggcorrplot)
library(lme4)
library(climateExtract)
library(MuMIn)

setwd(dir="C:/Users/anonym/Documents/plast_adaptation/data/data_publie_Ramakers")
tab=fread("Data_PeriodSpecificReactionNorm_and_Selection.txt",header=T)

#extract temperature for the coordinates of the site: 52.383333, 5.850000
# Points=data.frame(site_id="Ramaker",longitude=5.850000,latitude=52.383333)
# space= sf::st_as_sf(Points,coords=c('longitude','latitude'),crs="+proj=longlat +datum=WGS84 +no_defs")
# temp=NULL
# for(i in 1950:2019){
# climate_data <- extract_nc_value(i,i,local_file=T,clim_variable='mean temp',grid_size=0.1,file_path="tg_ens_mean_0.1deg_reg_v23.1e.nc",
# spatial_extent = space)
# aggrega=temporal_aggregate(x = climate_data,agg_function = "mean",variable_name = "average temp",time_step = "window",
# win_length = 3)
# bidon=as.data.frame(apply(as.data.frame(aggrega,xy=T)[,-c(1,2)],2,mean,na.rm=T))
# names(bidon)[1]="value"
# bidon$date_extract=as.Date(gsub("X","",row.names(bidon)),format="%Y.%m.%d")
# temp=rbind(temp,bidon)
# print(i)
# }

#OR use temperature data from the paper (http://climexp.knmi.nl/data/bteca2563.dat):
temp=fread("temperature_data.txt")
temp$date_extract=as.Date(paste(temp$y,temp$Month,temp$day,sep="."),format="%Y.%m.%d")


#CALCULATE TEMPERATURES INDICES
temp$jour=yday(temp$date_extract)
temp$y=year(temp$date_extract)
couples=data.frame(deb=rep(seq(0,365,15),2),fin= c(seq(0,365,15)+30,seq(0,365,15)+60))
couples$fin[which(couples$fin>365)]=couples$fin[which(couples$fin>365)]-365
couples=subset(couples,deb<360)

for(i in 1:nrow(couples)){
#subseter la période voulue
if(couples$fin[i]<=couples$deb[i]){
tab2=subset(temp,jour>couples$deb[i] | jour<=couples$fin[i])
}else{
temp2=subset(temp,jour>couples$deb[i] & jour<=couples$fin[i])}
if(couples$fin[i]<=couples$deb[i]){
temp2$y[which(temp2$jour>couples$deb[i])]=temp2$y[which(temp2$jour>couples$deb[i])]+1
}
obj=temp2 %>% dplyr::group_by(y) %>% dplyr::summarise(temperature=mean(value,na.rm=T))
bidon=melt(obj,id.var="y")
bidon=cbind(bidon,Points)
bidon$indice=paste("temp",couples$deb[i],couples$fin[i],sep="_")
if(i==1){tempf=bidon}else{tempf=rbind(tempf,bidon)}
print(i)}
names(tempf)[1]=c("Annee")
tempf=tempf %>% dplyr::group_by(site_id,indice) %>% dplyr::mutate(moy=mean(value,na.rm=T))
tempf=tempf %>% dplyr::group_by(Annee,indice) %>% filter (!duplicated(latitude,longitude))
tempf$value=tempf$value-tempf$moy
tempf2=reshape2::dcast(tempf,latitude+longitude+Annee~indice,value.var="value")
write.table(tempf2,"annual_mean_by_indice_ramakers.txt",sep="\t",row.names=F)

##################################################################################################################################
####LOAD DATA AND MERGE IT WITH TEMPERATURE INDICES
library(data.table)
library(dplyr)
library(raster)
library(Hmisc)
library(car)
library(doBy)
library(ggcorrplot)
library(lme4)
library(climateExtract)
library(MuMIn)
library(gridExtra)
library(cowplot)

setwd(dir="C:/Users/anonym/Documents/plast_adaptation/data/data_publie_Ramakers")
tab=fread("Data_PeriodSpecificReactionNorm_and_Selection.txt",header=T)
tab=subset(tab,Year>=1973)
tempf2=fread("annual_mean_by_indice_ramakers.txt",sep="\t")
tab=merge(tab,tempf2,by.x="Year",by.y="Annee")

####Select sensitivity of the optimal reaction norm to temperature (i.e., plasticity)
mati=cor(tab[,names(tab)[grep("temp_",names(tab),fixed=T)],with=F],use="complete.obs")
mati2=matrix(NA,ncol(mati),nrow(mati))
mati2[mati>0.8]=FALSE
mati2[mati<=0.8]=TRUE
mati2[upper.tri(mati2,diag=T)]=NA
colnames(mati2)=colnames(mati)
rownames(mati2)=rownames(mati)

liste=data.frame(model=names(tab)[grep("temp_",names(tab))],AICc=NA)
for(i in 1:nrow(liste)){
form=as.formula(paste0("LayDate~",liste$model[i],"+(1|Year)"))
model=lmer(form,data=tab)
liste$AICc[i]=AICc(model)
}


#FIT MODELS:
tab=as.data.frame(tab)
tab[,names(tab)[grep("temp_",names(tab),fixed=T)]]=as.data.frame(apply(tab[,names(tab)[grep("temp_",names(tab),fixed=T)]],2,scale))
mean(tab$Year)
tab$temp2=scale(tab$temp,scale=F)
tab$Year2=scale(tab$Year,scale=F)

#WITH TEMPERATURE INDEX FROM RAMAKERS ET AL.
model=lmer(LayDate~temp2*Year2+(1|Year),data=tab)
summary(model)

#WITH TEMPERATURE INDEX SELECTED HERE
form=as.formula(paste0("LayDate~(",liste$model[which.min(liste$AICc)],")*Year2+(1|Year)"))
model2=lmer(form,data=tab)
summary(model2)

adapt=data.frame(est=c("mean","lwr","upr"),varia=c("1988-2016"),
Ramakers=c(-2.34,-4.20,-0.48),sametemp=c(fixef(model)["Year2"]+
c(0,-1.96*coef(summary(model))["Year2","Std. Error"],1.96*coef(summary(model))["Year2","Std. Error"]))*
length(unique(tab$Year[tab$Year>=1988])),
difftemp=c(rep(fixef(model2)["Year2"],3)+c(0,-1.96*coef(summary(model2))["Year2","Std. Error"],
1.96*coef(summary(model2))["Year2","Std. Error"]))*length(unique(tab$Year[tab$Year>=1988])))
adapt$type="evolution of the reaction norm intercept"

plast=data.frame(est=c("mean","lwr","upr"),varia=c("1973-2016"),
Ramakers=c(-3.28,-3.92,-2.68),sametemp=c(fixef(model)["temp2"]+
c(0,-1.96*coef(summary(model))["temp2","Std. Error"],1.96*coef(summary(model))["temp2","Std. Error"]))*sd(tab$temp),
difftemp=fixef(model2)[paste(liste$model[which.min(liste$AICc)])]+
c(0,-1.96*coef(summary(model2))[paste(liste$model[which.min(liste$AICc)]),"Std. Error"],
1.96*coef(summary(model2))[paste(liste$model[which.min(liste$AICc)]),"Std. Error"])*sd(tab[,liste$model[which.min(liste$AICc)]]))
plast$type="plasticity"

evolplast=data.frame(est=c("mean","lwr","upr"),varia=c("1988-2016"),
Ramakers=c(-0.04,-0.09,0.02),sametemp=c(fixef(model)["temp2:Year2"]+
c(0,-1.96*coef(summary(model),)["temp2:Year2","Std. Error"],1.96*coef(summary(model))["temp2:Year2","Std. Error"]))*
sd(tab$temp)*length(unique(tab$Year[tab$Year>=1988])),
difftemp=c(fixef(model2)[paste0(liste$model[which.min(liste$AICc)],":Year2")])+
c(0,-1.96*coef(summary(model2))[paste0(liste$model[which.min(liste$AICc)],":Year2"),"Std. Error"],
1.96*coef(summary(model2))[paste0(liste$model[which.min(liste$AICc)],":Year2"),"Std. Error"])*
sd(tab[,liste$model[which.min(liste$AICc)]])*length(unique(tab$Year[tab$Year>=1988])))
evolplast$type="evolution of plasticity"

b=rbind(adapt,plast,evolplast)
b2=dcast(melt(b,id.vars=c("type","est","varia")),type+variable+varia~est)
b2$variable=as.character(b2$variable)
b2$variable[b2$variable=="Ramakers"]="Estimate from Ramakers et al."
b2$variable[b2$variable=="sametemp"]="Our method applied with\nRamakers et al. temperature index"
b2$variable[b2$variable=="difftemp"]="Our method applied while\nre-estimating the best temperature index"
b2$variable=factor(b2$variable,levels=c("Estimate from Ramakers et al.",
"Our method applied with\nRamakers et al. temperature index",
"Our method applied while\nre-estimating the best temperature index"))
b2$varia[b2$type=="plasticity"]="1973-2016"

pl1=ggplot(data=subset(b2,type=="evolution of the reaction norm intercept"),aes(x=variable,y=mean,color=varia))+
geom_hline(yintercept=0,linetype="dashed")+geom_point(size=4)+
geom_errorbar(aes(ymin=lwr,ymax=upr),width=0,alpha=0.7)+theme_bw()+ggtitle("a")+
theme(plot.title=element_text(size=14,face="bold"),panel.grid=element_blank(),
axis.title.y=element_blank(),axis.title.x=element_text(size=10),legend.position="bottom",
legend.title=element_text(size=10))+ylab("Evolution of the reaction norm intercept (days)")+
coord_flip()+labs(color="Period on which changes are estimated:")+
scale_color_manual(values=c("black","grey"),drop=F,breaks=c("1988-2016","1973-2016"))+guides(colour = guide_legend(title.position="top"))

pl2=ggplot(data=subset(b2,type=="plasticity"),aes(x=variable,y=mean,color=varia))+
geom_hline(yintercept=0,linetype="dashed")+geom_point(size=4)+
geom_errorbar(aes(ymin=lwr,ymax=upr),width=0,alpha=0.7)+theme_bw()+ggtitle("b")+
theme(plot.title=element_text(size=14,face="bold"),panel.grid=element_blank(),
axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=10),legend.position="bottom",
legend.title=element_text(size=10))+ylab("Phenotypic plasticity to temperature (days/°C)")+
coord_flip()+labs(color="Period used for estimation:")+
scale_color_manual(values=c("grey"))+guides(colour = guide_legend(title.position="top"))

pl3=ggplot(data=subset(b2,type=="evolution of plasticity"),aes(x=variable,y=mean,color=varia))+
geom_hline(yintercept=0,linetype="dashed")+geom_point(size=4)+
geom_errorbar(aes(ymin=lwr,ymax=upr),width=0,alpha=0.7)+theme_bw()+ggtitle("c")+
theme(plot.title=element_text(size=14,face="bold"),panel.grid=element_blank(),
axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=10),legend.position="bottom",
legend.title=element_text(size=10))+ylab("Evolution of phenotypic plasticity (days/°C)")+
coord_flip()+labs(color="Period on which changes are estimated:")+
scale_color_manual(values=c("black","grey"))+guides(colour = guide_legend(title.position="top"))

pl4=ggplot(data=b2,aes(x=variable,y=mean,color=varia))+
geom_hline(yintercept=0,linetype="dashed")+geom_point(size=4)+
geom_errorbar(aes(ymin=lwr,ymax=upr),width=0,alpha=0.7)+theme_bw()+ggtitle("a")+
theme(plot.title=element_text(size=14,face="bold"),panel.grid=element_blank(),
axis.title.y=element_blank(),axis.title.x=element_text(size=10),legend.position="bottom",
legend.title=element_text(size=10))+ylab("Evolution of the reaction norm intercept (days)")+
coord_flip()+labs(color="Period on which changes are estimated:")+
scale_color_manual(values=c("grey","black"),drop=F,breaks=c("1973-2016","1988-2016"))+guides(colour = guide_legend(title.position="top"))

leg <- ggpubr::as_ggplot(cowplot::get_legend(pl4))
pl1=pl1+theme(legend.position="none")
pl2=pl2+theme(legend.position="none")
pl3=pl3+theme(legend.position="none")
blank=ggplot() + theme_void()


plot_grid(pl1,pl2,pl3,blank,leg,blank,ncol=3,rel_widths=c(1.7,1,1),rel_heights=c(6,1),align="h")

setwd(dir="C:/Users/anonym/Documents/plast_adaptation")
pdf("figureS5.pdf",height=4,width=12)
plot_grid(pl1,pl2,pl3,blank,leg,blank,ncol=3,rel_widths=c(1.7,1,1),rel_heights=c(6,1),align="h")
dev.off();