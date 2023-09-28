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
library(ggforce)
library(grid)


setwd(dir="C:/Users/anonym/Documents/plast_adaptation/data/variations")
resf=NULL
lili=list.files()
for(i in lili){
bidon=fread(i,sep="\t")
resf=rbind(resf,bidon)
}


b=resf %>% dplyr::group_by(Speciesgen,varia) %>% dplyr::summarise(INLA=sd(mean)/abs(mean(mean)),LMER=sd(Estimate)/abs(mean(Estimate)))
b2=melt(b,id.vars=c("Speciesgen","varia"))

b2$varia=as.character(b2$varia)
b2$varia[b2$varia=="beta0"]="Intercept"
b2$varia[b2$varia=="Annee2"]="Year effect"
b2$varia[grep("temp_",b2$varia,fixed=T)]="Phenotypic plasticity"

pl1=ggplot(data=b2,aes(x=varia,y=value,color=variable))+geom_boxplot()+theme_bw()+
theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),axis.text.x=element_text(angle=45,hjust=1),
plot.title=element_text(size=14,face="bold"),axis.title.x=element_blank(),legend.position="left")+ylab("Coefficient of variation over 3 runs")+
labs(color="Model:")+
facet_zoom(y=TRUE,ylim =c(0,5),split = TRUE)


pl1

setwd(dir="C:/Users/anonym/Documents/plast_adaptation")
png("figure_S7.png",height=800,width=900,res=120)
pl1
dev.off();

