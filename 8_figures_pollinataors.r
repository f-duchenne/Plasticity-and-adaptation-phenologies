setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data/resultats")
library(doBy)
library(geosphere)
library(dplyr)
library(fitdistrplus)
library(data.table)
library(lem4)
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

resf=resf[,which(!duplicated(names(resf))),with=F]
resf=merge(resf,unique(liste[,c("Speciesgen","ORDRE","FAMILLE","Species")]),by="Speciesgen")

max(abs(resf$cor),na.rm=T)
min(resf$nb_annee)
min(resf$nbdata)

resf$var_std=resf$STD_maille/abs(resf$Estimate)

resf$varia2=as.character(resf$varia)
resf$varia2[resf$varia2=="Annee2"]="Year effect"
resf$varia2[grep("temp_",resf$varia2,fixed=T)]="Phenotypic plasticity"

boxplot(var_std~varia2,data=subset(resf,abs(var_std)<10))

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


###########################################################
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

######### FIGURE S3 & S4
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
scale_x_latitude()+scale_fill_manual(values=c("#3581B8","#FCB07E"))+scale_x_longitude(name="Average latitude",ticks =1)

grid.arrange(pl1,pl2,pl3,pl4,ncol=4,widths=c(1,1,1,1.2))
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
png("figure_S3.png",height=600,width=1600,res=120)
grid.arrange(pl1,pl2,pl3,pl4,ncol=4,widths=c(1,1,1,1.2))
dev.off();

obj=insight::get_variance(model)
b=data.frame(varia=c("Taxonomic (random effects)","Saptio-temporal (fixed effects)","Residual"),value=c(obj$var.random,obj$var.fixed,obj$var.residual))
b$value=b$value/sum(b$value)
pl1=ggplot(data=b,aes(x=varia,y=value))+geom_bar(stat="identity",col="white")+theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
plot.title=element_text(size=14,face="bold"),axis.title.x=element_blank())+ylab("Percentage of variance in phenotypic plasticity explained")+
scale_y_continuous(labels=percent)+coord_cartesian(expand=F)

png("figure_S4.png",height=800,width=700,res=120)
pl1
dev.off();

######### FIGURE 5
#calculate contributions in day/year
resf[,contrib:=Estimate*trend]
resf[type=="adapt",contrib:=Estimate]

resf[,contrib_lwr:=apply(cbind(lwr*trend,upr*trend),1,min)]
resf[type=="adapt",contrib_lwr:=lwr]

resf[,contrib_upr:=apply(cbind(lwr*trend,upr*trend),1,max)]
resf[type=="adapt",contrib_upr:=upr]

b=subset(resf,type %in% c("plast","adapt")) %>% dplyr::group_by(Speciesgen,ORDRE,FAMILLE,Species,mfd) %>% dplyr::summarise(contrib_plast=sum(contrib[type=="plast"],na.rm=T),contrib_plast_lwr=sum(contrib_lwr[type=="plast"],na.rm=T),
contrib_plast_upr=sum(contrib_upr[type=="plast"],na.rm=T),
contrib_adapt=sum(contrib[type=="adapt"],na.rm=T),contrib_adapt_lwr=sum(contrib_lwr[type=="adapt"],na.rm=T),contrib_adapt_upr=sum(contrib_upr[type=="adapt"],na.rm=T))

b$Est.signi_plast=">0.05"
b$Est.signi_plast[b$contrib_plast_lwr>0]="<0.05"
b$Est.signi_plast[b$contrib_plast_upr<0]="<0.05"
b$Est.signi_adapt=">0.05"
b$Est.signi_adapt[b$contrib_adapt_lwr>0]="<0.05"
b$Est.signi_adapt[b$contrib_adapt_upr<0]="<0.05"
b$Est.signi_plast=factor(b$Est.signi_plast,c(">0.05","<0.05"))

b$contrib_adapt_lwr[b$contrib_adapt_lwr<(min(b$contrib_adapt)-0.05)]=min(b$contrib_adapt)-0.05
b$contrib_adapt_upr[b$contrib_adapt_upr>(max(b$contrib_adapt)+0.05)]=max(b$contrib_adapt)+0.05
b$contrib_plast_lwr[b$contrib_plast_lwr<(min(b$contrib_plast)-0.05)]=min(b$contrib_plast)-0.05
b$contrib_plast_upr[b$contrib_plast_upr>(max(b$contrib_plast)+0.05)]=max(b$contrib_plast)+0.05

#plast numbers
nrow(subset(b,contrib_plast<0 & Est.signi_plast=="<0.05"))/nrow(b)
nrow(subset(b,contrib_plast>0 & Est.signi_plast=="<0.05"))/nrow(b)
#adapt numbers
nrow(subset(b,contrib_adapt<0 & Est.signi_adapt=="<0.05"))/nrow(b)
nrow(subset(b,contrib_adapt>0 & Est.signi_adapt=="<0.05"))/nrow(b)

pl1=ggplot(data=b,aes(x=contrib_plast,fill=Est.signi_plast))+geom_histogram(col="black")+
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
xlab("")+ylab("Number of species")+ggtitle("b",subtitle="Year effect")+coord_cartesian(expand=F)


obj=expand.grid(seq(min(b$contrib_adapt)-0.05,max(b$contrib_adapt)+0.05,0.005),seq(min(b$contrib_plast)-0.05,max(b$contrib_plast)+0.05,0.005))
obj$value=obj[,1]+obj[,2]

pl=ggplot()+
geom_raster(data=obj,aes(x=Var1,y=Var2,fill=value))+
scale_fill_gradient2(name="Resulting\npheno. shift\n(day/year)")+
geom_vline(xintercept=0,col="black")+geom_hline(yintercept=0,col="black")+geom_abline(intercept=0,slope=-1,linetype="dashed")+
geom_errorbarh(data=b,aes(y=contrib_plast,x=contrib_adapt,xmin=contrib_adapt_lwr,xmax=contrib_adapt_upr),alpha=0.1)+
geom_errorbar(data=b,aes(y=contrib_plast,x=contrib_adapt,ymin=contrib_plast_lwr,ymax=contrib_plast_upr),alpha=0.1)+
geom_point(data=b,aes(y=contrib_plast,x=contrib_adapt),col="black",alpha=0.7)+
theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"))+
ylab("Contribution of phenotypic plasticity to phenological shifts (day/year)")+
xlab("Contribution of year effect to phenological shifts (day/year)")+coord_fixed(ratio=1,expand=F)+ggtitle("c")


haut=grid.arrange(pl1,pl2,ncol=2,bottom="Contribution to phenological shifts (day/year)",left="Number of species")

layout=rbind(c(1,1,1,1),c(3,3,3,4))
setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
pdf("figure5.pdf",height=7,width=10)
grid.arrange(haut,pl,ncol=1,heights=c(1,2.6))
dev.off();


################## FIGURE 6

b$contrib_to_total=abs(b$contrib_adapt)/(abs(b$contrib_plast)+abs(b$contrib_adapt))
b$contrib_to_total_s=ifelse(b$Est.signi_adapt=="<0.05",abs(b$contrib_adapt),0)/(ifelse(b$Est.signi_plast=="<0.05",abs(b$contrib_plast),0)+
ifelse(b$Est.signi_adapt=="<0.05",abs(b$contrib_adapt),0))
boxplot(b[,c("contrib_to_total","contrib_to_total_s")])
b$direction="similar"
b$direction[sign(b$contrib_adapt)!=sign(b$contrib_plast)]="opposite"

b2=melt(b,id.vars=c("Speciesgen","mfd","direction"),measure.vars=c("contrib_to_total","contrib_to_total_s"))
b2=subset(b2,!is.na(value))

give.n <- function(x){
  return(data.frame(label = paste0("n=",length(x)),y=1.02)) 
  # experiment with the multiplier to find the perfect position
}

b2 %>% dplyr::group_by(variable) %>% dplyr::summarise(mean(value),sd(value))

pl1=ggplot(data=b2,aes(x=variable,y=value,color=direction))+#geom_violin(fill="dodgerblue3")+
geom_boxplot(width=0.5)+scale_x_discrete(labels=c("All\ncontributions","Only significant\ncontirbutions"))+
theme_bw()+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"),axis.line = element_line(colour = "black"),axis.title.x=element_blank(),axis.text.x=element_text(angle=0,hjust=0.5))+
ylab("Relative contribution of\nthe year effect to phenological shifts")+scale_y_continuous(labels=percent,lim=c(0,1.05))+
stat_summary(fun.data = give.n, geom = "text",position = position_dodge(width = 0.75),show_guide  = FALSE,size=3)+
coord_cartesian(expand=F)+scale_color_manual(values=c("darkslategray4","gold3"))+labs(col="Direction of\ncontributions")+ggtitle("a")


resf=resf %>% group_by(Speciesgen) %>% mutate(contrib_plast=sum(contrib[type=="plast"],na.rm=T),contrib_plast_lwr=sum(contrib_lwr[type=="plast"],na.rm=T),
contrib_plast_upr=sum(contrib_upr[type=="plast"],na.rm=T))
bb=subset(resf,type=="adapt")


var(bb$contrib_plast)
var(bb$Estimate)

pl2=ggplot(data=bb,aes(x=contrib_plast,y=longterm_mfd))+
theme_bw()+geom_vline(xintercept=0,col="black",linetype="dashed")+geom_hline(yintercept=0,col="black",linetype="dashed")+
geom_errorbarh(aes(xmin=contrib_plast_lwr,xmax=contrib_plast_upr),alpha=0.1)+
geom_errorbar(aes(ymin=longterm_mfd-1.96*longterm_mfd.se,ymax=longterm_mfd+1.96*longterm_mfd.se),alpha=0.1)+
geom_point(alpha=0.5)+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"),axis.line = element_line(colour = "black"))+
ylab("Long-term phenological shifts (day/year)")+xlab("Contribution of phenotypic plasticity(day/year)")+
coord_cartesian(expand=F)+ggtitle("b")+stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,label.x=-0.92)

pl3=ggplot(data=bb,aes(x=Estimate,y=longterm_mfd))+
theme_bw()+geom_vline(xintercept=0,col="black",linetype="dashed")+geom_hline(yintercept=0,col="black",linetype="dashed")+
geom_errorbarh(aes(xmin=Estimate-1.96*sde,xmax=Estimate+1.96*sde),alpha=0.1)+
geom_errorbar(aes(ymin=longterm_mfd-1.96*longterm_mfd.se,ymax=longterm_mfd+1.96*longterm_mfd.se),alpha=0.1)+
geom_point(alpha=0.55)+
theme(panel.grid=element_blank(),strip.background = element_blank(),panel.border=element_blank(),
plot.title=element_text(size=14,face="bold"),axis.line = element_line(colour = "black"))+
ylab("Long-term phenological shifts (day/year)")+xlab("Contribution of year effect(day/year)")+
coord_cartesian(expand=F)+ggtitle("c")+stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,label.x=0.1,label.y=-1)


plot_grid(pl1,pl2,pl3,align="hv",ncol=3,rel_widths=c(1.2,1,1))

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation")
pdf("figure6.pdf",height=5,width=12)
plot_grid(pl1,pl2,pl3,align="hv",ncol=3,rel_widths=c(1.2,1,1))
dev.off();



