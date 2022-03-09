pkgs <- c("ggplot2", "mgcv", "MASS","car","gplots","doBy","dplyr","lme4","mgcv","numDeriv","nleqslv","lubridate","reshape2","Hmisc",
"fitdistrplus","geosphere","pastecs","SDMtools","jtools","gridExtra","stats","ggmap")
lapply(pkgs, require, character.only = TRUE,quietly=T)
setwd(dir="C:/Users/Francois/Documents/papier 1 - données collection/data")
library(mgcv)
library(numDeriv)
library(nleqslv)
library(lubridate)
library(Hmisc)
library(stats)
library(pastecs)
library(SDMTools)
library(geosphere)
library(fitdistrplus)
library(reshape2)
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

# setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/résultats simus")
# dat=fread("trend_plasticite.txt")
# dat=subset(dat,param=="as.numeric(Annee)")
# hist(dat$Estimate)
# names(dat)[names(dat)=="Esp"]="Species"
# names(dat)[names(dat)=="Estimate"]="trend"
# setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data_occ_UK")
# occt=fread("trends.txt",sep="\t",header=T)
# resf=merge(dat,occt[,c("Species","Mean_growth_rate","Precision")],by="Species",all.x=T,all.y=F)
# setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/SUMMARY_TABLES")
# lili=list.files()
# vec=c()
# for(i in unique(dat$Species)){
# vec=c(vec,grep(i,lili))
# }
# lili=lili[vec]
# objf=NULL
# for(i in lili){
# obj=fread(i,header=T)
# obj=obj[obj$Region %in% c("UK","GB"),]
# obj$Mean2=car::logit(obj$Mean)
# obj$Lower_CI2=car::logit(obj$Lower_CI)
# obj$Upper_CI2=car::logit(obj$Upper_CI)
# obj$wei=obj$Upper_CI2-obj$Lower_CI2
# obj$wei[obj$wei<1e-7]=0.001
# model=lm(Mean2~Year,data=obj,weights=1/wei)
# obj$trend_ang=model$coeff[2]
# obj$trend_err=summary(model)$coeff[2,2]
# objf=rbind(obj,objf)
# }

# obf2=unique(objf[,c("Species","trend_ang","trend_err")])
# obf2=objf %>% dplyr::group_by(Species,trend_ang,trend_err) %>% dplyr::summarise(abond=mean(Mean))
# resf=merge(resf,obf2,by="Species",all.x=T,all.y=F)
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
plast=fread("plast_adapt.txt",sep="\t",header=T)
plast=subset(plat,!is.na(core))
plast$data_number=500
plast$data_number[plast$nbdata>=750]=1000
plast$data_number[plast$nbdata>=1250]=1500
plast$data_number[plast$nbdata>=2000]=2500
plast$data_number=paste0("n = ",plast$data_number)
plast$correlation=plast$core
plast$nbit=plast$nb_annee
plast$data_number=factor(plast$data_number,c("n = 500","n = 1000","n = 1500","n = 2500"))

setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/résultats simus")
library(plyr)
res2b=fread("resultats_simulations_sans_g_evol_sans_plast_evol500.txt",header=T,sep="\t")
res2c=fread("resultats_simulations_sans_g_evol_sans_plast_evol1000.txt",header=T,sep="\t")
res2c2=fread("resultats_simulations_sans_g_evol_sans_plast_evol1500.txt",header=T,sep="\t")
res2d=fread("resultats_simulations_sans_g_evol_sans_plast_evol2500.txt",header=T,sep="\t")
res2=rbind(res2d,res2b,res2c,res2c2)
res2$diffea=abs(res2$adapt-res2$t_adapt)
res2$diffep=abs(res2$plast-res2$t_plast)
res2$diffei=abs(res2$evol_plast-0)
res2$data_number=paste("n = ",round_any(res2$data_number,100),sep="")
res2$data_number=factor(res2$data_number,c("n = 500","n = 1000","n = 1500","n = 2500"))
res2$correlation=plyr::round_any(res2$correlation,0.02)


b=summaryBy(diffea~correlation+nbit+data_number,data=subset(res2,correlation>(-0.1) & correlation<1),
FUN=quantile,probs=0.95,keep.names=T)
b$diffea[b$diffea>=0.4]=0.4
pl1=ggplot(data=b,aes(x=correlation,y=nbit))+
geom_raster(interpolate=T,aes(fill=diffea,col=diffea))+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),axis.title.x=element_blank(),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colours=rev(c(rainbow(10)[1:8],"darkblue")),na.value="white",limits=c(0,max(b$diffea)))+
labs(fill="Error in\nday/year")+ylab("")+ggtitle("a")+
geom_point(data=plast,alpha=0.7,size=0.8,col="red")

b=summaryBy(diffep~correlation+nbit+data_number,data=subset(res2,correlation>(-0.1) & correlation<1),
FUN=quantile,probs=0.95,keep.names=T)
b$diffep[b$diffep>=10]=10
pl2=ggplot(data=b,aes(x=correlation,y=nbit))+
geom_raster(interpolate=T,aes(fill=diffep,col=diffep))+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.x=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colours=rev(c(rainbow(10)[1:8],"darkblue")),na.value="white",limits=c(0,max(b$diffep)))+
labs(fill="Error in\nday/°C")+ylab("Number of years")+ggtitle("b")+
geom_point(data=plast,alpha=0.7,size=0.8,col="red")


b=summaryBy(diffei~correlation+nbit+data_number,data=subset(res2,correlation>(-0.1) & correlation<1),
FUN=quantile,probs=0.95,keep.names=T)
b$diffei[b$diffei>=1]=1
pl3=ggplot(data=b,aes(x=correlation,y=nbit))+
geom_raster(interpolate=T,aes(fill=diffei,col=diffei))+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),,
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colours=rev(c(rainbow(10)[1:8],"darkblue")),na.value="white",limits=c(0,max(b$diffei)))+
labs(fill="Error in\nday/°C/year")+xlab("Time-temperature correlation")+ylab("")+ggtitle("c")+
geom_point(data=plast,alpha=0.7,size=0.8,col="red")


cowplot::plot_grid(pl1,pl2,pl3,nrow=3,align = "v")


pdf("figure_2.pdf",width=11,height=8)
cowplot::plot_grid(pl1,pl2,pl3,nrow=3,align = "v")
dev.off();
