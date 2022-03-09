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
library(viridis)
library(MuMIn)
library(moments)
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
donnf2=fread("dataf_avec_elev.txt",sep="\t")
tabtemp=fread("annual_mean_by_indice_and_by_point_mat.txt",sep="\t")
tabtemp=tabtemp[,-which(names(tabtemp)=="temp_360_85"),with=F]
tabtemp$latitude=round(tabtemp$latitude,digits=5)
tabtemp$longitude=round(tabtemp$longitude,digits=5)
donnf2$latitude=round(donnf2$latitude,digits=5)
donnf2$longitude=round(donnf2$longitude,digits=5)
donnf2=merge(donnf2,tabtemp,by=c("latitude","longitude","Annee"))
liste=unique(donnf2[,c("Speciesgen","Species","FAMILLE","ORDRE")])
liste[,names(tabtemp)[4:ncol(tabtemp)]]=NA
liste$rsq=NA
donnf2=subset(donnf2,!is.na(elev))
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/résultats simus")
liste=fread("selec_temp_var_beta.txt",sep="\t")
liste=subset(liste,nb_pre1990>=50 & nb_post1990>=500)
vec=names(tabtemp)[4:ncol(tabtemp)]
nb=liste[,c("Speciesgen","Species","FAMILLE","ORDRE","nb_pre1990","nb_post1990")]


for(i in 1:nrow(liste)){
bidon=subset(donnf2,Speciesgen==liste$Speciesgen[i] & !is.na(temp_0_90))
a=0
bidon$elev2=bidon$elev
bidon$elev=(bidon$elev-mean(bidon$elev,na.rm=T))/1000
bidon$Annee=bidon$Annee-1960
varia=vec[which(!is.na(liste[i,vec,with=F]))]
if(length(unique(bidon$maille))>1){
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+Annee+(1+",paste(varia,collapse="+"),"+Annee","|maille)+(1|Annee)"))
model=fitme(form,data=bidon)
}else{
form=formula(paste0("Jour.de.collecte~",paste(varia,collapse="+"),
"+elev+Annee+(1|Annee)"))
model=fitme(form,data=bidon)
}	
bidon$resid=residuals(model)
pl1=ggplot(data=bidon,aes(x=X,y=Y,color=sqrt(abs(resid))*sign(resid)))+
geom_point(size=2)+scale_color_viridis()
res=as.data.frame(summary(model)$beta_table)
res$varia=rownames(res)
res$core=NA
res$core[res$varia %in% varia]=t(cor(bidon[,c("Annee",varia),with=F])[-1,1])
res$Var=NA
res$Var[res$varia %in% summary(model)$lambda_table$Term]=summary(model)$lambda_table$Var
res$Speciesgen=liste$Speciesgen[i]
res$latmean=mean(bidon$Y)
res$longmean=mean(bidon$X)
res$elevmean=mean(bidon$elev2,na.rm=T)
res$mfd=mean(bidon$Jour.de.collecte)
res$nb_annee=length(unique(bidon$Annee))
res$nbdata=nrow(bidon)
res$trend=NA
bidon$Annee2=plyr::round_any(bidon$Annee,2)
b=bidon %>% dplyr::group_by(Speciesgen,Annee2) %>% dplyr::summarise(skew=skewness(resid),nbd=length(resid),
kurt=kurtosis(resid))
modelskew=lm(skew~Annee2,data=b,weights=sqrt(nbd))
res$trend_skew=modelskew$coeff[2]
res$trend_skew_err=summary(modelskew)$coeff[2,2]
modelkurt=lm(kurt~Annee2,data=b,weights=sqrt(nbd))
res$trend_kurt=modelkurt$coeff[2]
res$trend_kurt_err=summary(modelkurt)$coeff[2,2]

resqual=bidon %>% dplyr::group_by(Speciesgen,Species,FAMILLE,ORDRE,maille,ctr_lon,ctr_lat) %>%
dplyr::summarise(n=n())
if(length(unique(bidon$maille))>1){bibi=as.data.frame(ranef(model)[[1]])
}else{bibi=as.data.frame(matrix(0,1,length(varia)+2))}
bibi=bibi+do.call("rbind",replicate(nrow(bibi),as.data.frame(res$Estimate[res$varia %in% colnames(ranef(model)[[1]])]),
))
if(length(unique(bidon$maille))>1){
for(j in varia){
mod=fitme(formula(paste0(j,"~Annee+(Annee","|maille)")),
data=unique(bidon[,c("Annee","latitude","longitude","maille",j),with=F]))
res$trend[res$varia==j]=mod$fixef[2]
bibi1=as.data.frame(ranef(mod)[[1]][,2])+mod$fixef[2]
names(bibi1)=paste0("trend_",j)
bibi=cbind(bibi,bibi1)
}
}else{
for(j in varia){
mod=lm(formula(paste0(j,"~Annee")),
data=unique(bidon[,c("Annee","latitude","longitude",j),with=F]))
res$trend[res$varia==j]=mod$coeff[2]
}
}

bibi$maille=as.numeric(rownames(bibi))
resqual=merge(resqual,bibi,by="maille")
print(ggplot(data=resqual,aes(x=ctr_lon,y=ctr_lat,col=Annee,size=n))+geom_point()+scale_color_gradient2())
resqual=melt(resqual,id.vars=c("Speciesgen","Species","FAMILLE","ORDRE","maille","ctr_lon","ctr_lat","n"))
if(i==1){
resf=res
resqualf=resqual
}else{
resf=rbind(resf,res)
resqualf=rbind(resqualf,resqual)
}
print(i)
}
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
fwrite(resqualf,"plast_adapt_par_maille.txt",sep="\t")
fwrite(resf,"plast_adapt.txt",sep="\t")

######################################
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
elev=fread("elev_by_species.txt",header=T)
plast=fread("plast_adapt.txt",sep="\t",header=T)
liste=fread("nb_data.txt",sep="\t")
plast=merge(plast,liste,by="Speciesgen")
resf=merge(plast,elev,by="Speciesgen")
#resf=merge(resf,plast,by="Species")
#resf=resf[grep("_NA",resf$Speciesgen),]
resf$type=NA
resf$type[grep("temp",resf$varia)]="plast"
resf$type[grep("Annee",resf$varia)]="adapt"
resf$saison=NA
resf$saison[resf$varia %in% c("temp_300_25","temp_330_55","temp_0_90")]="hiver"
resf$saison[resf$varia %in% c("temp_60_150","temp_120_210","temp_180_270","temp_30_120",
"temp_270_360","temp_90_180","temp_210_300","temp_240_330")]="été"
resf$saison[resf$varia %in% c("temp_60_150","temp_30_120","temp_90_180")]="printemps"
resf$saison[resf$varia %in% c("temp_240_330","temp_270_360","temp_210_300")]="automne"
resf$varia=factor(resf$varia,c("temp_300_25","temp_330_55","temp_0_90","temp_30_120","temp_60_150","temp_90_180",
"temp_120_210","temp_150_240","temp_180_270","temp_210_300","temp_240_330","temp_270_360"))
resf$varia=gsub("temp_","",resf$varia)
resf$varia=gsub("_","-",resf$varia)
resf$varia=factor(resf$varia,c("300-25","330-55","0-90","30-120","60-150","90-180",
"120-210","150-240","180-270","210-300","240-330","270-360"))
resf$Est.lwr=resf$Estimate-1.96*resf$"Cond. SE"
resf$Est.upr=resf$Estimate+1.96*resf$"Cond. SE"
resf$Est.signi=">0.05"
resf$Est.signi[resf$Est.lwr>0]="<0.05"
resf$Est.signi[resf$Est.upr<0]="<0.05"
resf$date_fin=as.numeric(sapply(strsplit(as.character(resf$varia),"-"), function(l) l[[length(l)]]))
names(resf)[names(resf)=="Cond. SE"]="Cond.SE"
names(resf)[names(resf)=="Std. Error"]="Std.Error"

b=subset(resf,!is.na(varia)) %>% group_by(Speciesgen,Est.signi) %>% dplyr::count()
b %>% group_by(Est.signi,n) %>% dplyr::summarise(nb=n(),perc=n()/length(unique(resf$Speciesgen)))


setwd(dir="C:/Users/Francois/Documents/land use change/European_costlines")
shp <- st_read(".", "Europe_coastline")
bidon=sf::st_transform(sf::st_as_sf(resf,coords = c('longmean','latmean'),crs=st_crs(shp)),crs=4326)
bidon2=as.data.frame(st_coordinates(bidon))
resf$latitude=bidon2[,2]
resf$longitude=bidon2[,1]
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
coucou=c("magenta2",brewer.pal(5,"Set2")[5],"red","dodgerblue4")
####### PLAST ######
library(ggraph)
summary(subset(resf,varia=="90-180")$Estimate)
sd(subset(resf,varia=="90-180")$Estimate)/sqrt(nrow(subset(resf,varia=="90-180")))
resf$longmean2=scale(resf$longitude,scale=F)
resf$latmean2=scale(resf$latitude,scale=F)
bidon=subset(resf,!is.na(varia))
bidon$contrib=bidon$Estimate
mati=as.data.frame(dcast(bidon,ORDRE+Speciesgen+mfd~varia,value.var="Estimate"))
mati2=as.matrix(mati[,-c(1:3)])
rownames(mati2)=c(as.character(mati[,2]))
colnames(mati2)=names(mati[,-c(1:3)])
mati2[is.na(mati2)]=0
graph=graph_from_incidence_matrix(mati2, directed = FALSE,weighted =T)
plot(graph,vertex.label=NA,layout=layout_as_bipartite)
dipt=unique(mati$Speciesgen[mati$ORDRE=="Diptera"])
lep=unique(mati$Speciesgen[mati$ORDRE=="Lepidoptera"])
cole=unique(mati$Speciesgen[mati$ORDRE=="Coleoptera"])
hym=unique(mati$Speciesgen[mati$ORDRE=="Hymenoptera"])
tab_temp=unique(bidon[,c("varia","date_fin")])

V(graph)$mfd=NA
V(graph)$mfd=ifelse(V(graph)$name %in% tab_temp$varia,tab_temp$date_fin[match(tab_temp$varia,V(graph)$name)],
mati$mfd[match(mati$Speciesgen,V(graph)$name)])
V(graph)$ORDRE <- NA
V(graph)$ORDRE <- ifelse(V(graph)$name %in% dipt, 'Diptera',ifelse(V(graph)$name %in% lep, 'Lepidoptera',
ifelse(V(graph)$name %in% cole, 'Coleoptera',ifelse(V(graph)$name %in% hym, 'Hymenoptera'," Temperature indices"))))
#V(graph)$ORDRE=factor(V(graph)$ORDRE,c("Temperature indices","Coleoptera","Diptera","Hymenoptera","Lepidoptera"))
zeropos=abs(min(E(graph)$weight))/(max(E(graph)$weight)-min(E(graph)$weight))
graph1=subgraph(graph, which(V(graph)$ORDRE %in% c("Coleoptera","Diptera"," Temperature indices")))
pl1=ggraph(graph1, 'hive', axis=ORDRE,sort.by=mfd) + 
geom_edge_hive(aes(color=weight,alpha=abs(weight),width=abs(weight)))+
scale_edge_colour_gradientn(colours=c("firebrick4","pink","white","cadetblue1","dodgerblue4"),
values=c(0,zeropos-0.01,zeropos,zeropos+0.01,1),name="Effect (day/°C)",
limits=c(min(E(graph)$weight),max(E(graph)$weight)))+
scale_edge_width_continuous(range = c(0.1, 2),guide=F)+
scale_edge_alpha_continuous(range = c(0,1),guide=F)+
geom_axis_hive(color="black",alpha=1)+theme(panel.background=element_rect(fill="white"),
plot.title=element_text(size=14,face="bold"),legend.position="none")+
labs(color="Effect")+ggtitle("a")

graph2=subgraph(graph, which(V(graph)$ORDRE %in% c("Lepidoptera","Hymenoptera"," Temperature indices")))
pl2=ggraph(graph2, 'hive', axis=ORDRE,sort.by=mfd) + 
geom_edge_hive(aes(color=weight,alpha=abs(weight),width=abs(weight)))+
scale_edge_colour_gradientn(colours=c("firebrick4","pink","white","cadetblue1","dodgerblue4"),
values=c(0,zeropos-0.01,zeropos,zeropos+0.01,1),name="Effect (day/°C)",
limits=c(min(E(graph)$weight),max(E(graph)$weight)))+
scale_edge_width_continuous(range = c(0.1, 2),guide=F)+
scale_edge_alpha_continuous(range = c(0,1),guide=F)+
geom_axis_hive(color="black",alpha=1)+theme(panel.background=element_rect(fill="white"),
plot.title=element_text(size=14,face="bold"))+
labs(color="Effect")+ggtitle("b")

model=lmer(Estimate~longmean2+elev+latmean2+mfd*poly(date_fin,2)+(1|ORDRE)+(1|FAMILLE)+(1|Species)+
(1|Speciesgen),
data=bidon,weights=sqrt(1/Cond.SE))


b=ggeffect(model,c("date_fin","mfd"))
b$group=round(as.numeric(as.character(b$group)))
pl3=ggplot(b,aes(x=x,y=predicted,group=group,color=as.factor(group),fill=as.factor(group)))+
geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.5,col=NA)+
scale_color_manual(values=c(brewer.pal(9,"YlOrRd")[c(4,6,9)]))+
scale_fill_manual(values=c(brewer.pal(9,"YlOrRd")[c(4,6,9)]))+
geom_line()+theme_bw()+theme(plot.title=element_text(size=14,face="bold"),legend.position="right",
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("c")+
xlab("3-month temperature indices (in day of the year),\nfrom winter (left) to autumnal (right) temperatures")+
ylab("Predicted phenotypic plasticity (day/°C)\n")+labs(color="MFD\n(day of the year)",fill="MFD\n(day of the year)")+
scale_x_continuous(breaks=tab_temp$date_fin,labels=tab_temp$varia)

grid.arrange(pl1,pl2,pl3,layout_matrix=matrix(c(1,3,2,3),nrow=2),widths=c(1,1.5))

setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/data")
pdf(paste0("fig_plast_season.pdf"),width=8,height=7)
grid.arrange(pl1,pl2,pl3,layout_matrix=matrix(c(1,3,2,3),nrow=2),widths=c(1,1.5))
dev.off();


leg <- cowplot::get_legend(pl1)
pl1=pl1+theme(legend.position="none")
# Convert to a ggplot and print

resf=resf %>% dplyr::mutate(nsp=length(unique(Species))) %>% ungroup()
b=resf %>% dplyr::group_by(varia,ORDRE) %>% dplyr::summarise(perc=length(unique(Species)),count=unique(nsp))

pl2=ggplot(data=subset(b,!is.na(varia)),aes(x=varia,y=perc/count,col=NA,fill=ORDRE))+
geom_bar(stat="identity",width=0.5) + 
scale_y_continuous(breaks=c(0,0.25,0.5),labels=scales::percent(c(0,0.25,0.5)))+theme_bw()+
theme(plot.title=element_text(size=14,face="bold"),legend.position="right",legend.title=element_blank(),
axis.text.x=element_text(angle = 45,hjust=0.9))+
ggtitle("b")+xlab("3-month temperature indices (day of the year)")+ylab("Percentage of species responding")+
scale_color_manual(values=coucou)+
scale_fill_manual(values=coucou)
leg2 <- cowplot::get_legend(pl2)
pl2=pl2+theme(legend.position="none")

blan=patchwork::plot_spacer()+theme_void()
grid.arrange(pl1,pl2,leg,leg2,heights=c(1.5,1),widths=c(3.5,1),layout_matrix=matrix(c(1,2,3,4),nrow=2))


####### ADAPT ######

resf$Est.signi=factor(resf$Est.signi,c(">0.05","<0.05"))
b=subset(resf,type=="adapt" & Est.signi=="<0.05") %>% group_by(Speciesgen) %>% dplyr::count()
p=ggplot(data=subset(resf,type=="adapt"),aes(x=Estimate,fill=Est.signi))+geom_histogram(col="black")+
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
