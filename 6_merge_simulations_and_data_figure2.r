library(data.table)
library(rootSolve)
library(ggplot2)
library(doBy)
library(gridExtra)
library(MCMCglmm)
library(plot3D)
library(reshape2)
library(plyr)
library(Rmisc)
library(svMisc)
library(viridis)
library(ggplot2)
library(car)
library(lme4)
library(colorBlindness)
library(spaMM)

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/data")
liste=fread("selec_temp_var_beta.txt",sep="\t")
liste=subset(liste,!is.na(cor))
liste$nbdata=liste$nb_pre1990+liste$nb_post1990
liste$data_number=500
liste$data_number[liste$nbdata>=750]=1000
liste$data_number[liste$nbdata>=1250]=1500
liste$data_number[liste$nbdata>=2000]=2500
liste$data_number=paste0("n = ",liste$data_number)
liste=liste %>% dplyr::group_by(Speciesgen) %>% dplyr::mutate(correlation=max(abs(cor)))
liste$nbit=liste$nb_annee
liste$data_number=factor(liste$data_number,c("n = 500","n = 1000","n = 1500","n = 2500"))


nbit_vec=c(25,30,35,40,45,50,55)
nbdata_vec=c(500,1000,1500,2500)
mechanisms=c("without_adapt","without_plast","without_evol_plast","total","without_evolution")
setwd(dir="C:/Users/duchenne/Documents/plast_adaptation/data/resultats_simues")
tabtab=expand.grid(nbdata_vec,nbit_vec,mechanisms)
res2=NULL
for(i in 1:nrow(tabtab)){
nb_data=tabtab[i,1] 
nbit=tabtab[i,2]  
mecha=tabtab[i,3]
sim=fread(paste("resultats_simulations_sans_evol_de_g_",nb_data,"_",mecha,"_",nbit,".txt",sep=""))
sim$mechanism=as.character(sim$mechanism)
sim$mechanism=as.character(mecha)
res2=rbind(res2,sim)
}
library(plyr)
res2$diffea=abs(res2$adapt-res2$t_adapt)
res2$diffep=abs(res2$plast-res2$t_plast)
res2$diffei=abs(res2$evol_plast-0)
res2$data_number=paste("n = ",round_any(res2$data_number,100),sep="")
res2$data_number=factor(res2$data_number,c("n = 500","n = 1000","n = 1500","n = 2500"))
res2$correlation=plyr::round_any(res2$correlation,0.04)


b=summaryBy(diffea~correlation+nbit+data_number,data=subset(res2,correlation>=min(liste$correlation) & correlation<=1),
FUN=mean,na.rm=T,keep.names=T)
b$diffea[b$diffea>=0.4]=0.4
pl1=ggplot()+
geom_raster(data=b,aes(x=correlation,y=nbit,fill=diffea,col=diffea),interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),axis.title.x=element_blank(),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colors=Blue2DarkRed18Steps,na.value="black",limits=c(0,max(b$diffea)))+
labs(fill="Error in\nday/year")+ylab("")+ggtitle("a - Evolution of reaction norm's elevation")+
geom_point(data=liste,aes(x=correlation,y=nbit),size=1,col="black")+coord_cartesian(expand=F)

b=summaryBy(diffep~correlation+nbit+data_number,data=subset(res2,correlation>=min(liste$correlation) & correlation<=1),
FUN=mean,na.rm=T,keep.names=T)
b$diffep[b$diffep>=10]=10
pl2=ggplot()+
geom_raster(data=b,aes(x=correlation,y=nbit,fill=diffep,col=diffep),interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.x=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colors=Blue2DarkRed18Steps,na.value="black",limits=c(0,max(b$diffep)))+
labs(fill="Error in\nday/°C")+ylab("Number of years")+ggtitle("b - Phenotypic plasticity")+
geom_point(data=liste,aes(x=correlation,y=nbit),size=1,col="black")+coord_cartesian(expand=F)


b=summaryBy(diffei~correlation+nbit+data_number,data=subset(res2,correlation>=min(liste$correlation) & correlation<=1),
FUN=mean,na.rm=T,keep.names=T)
b$diffei[b$diffei>=1]=1
pl3=ggplot()+
geom_raster(data=b,aes(x=correlation,y=nbit,fill=diffei,col=diffei),interpolate=T,na.rm = TRUE)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),,
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colors=Blue2DarkRed18Steps,na.value="black",limits=c(0,max(b$diffei)))+
labs(fill="Error in\nday/°C/year")+xlab("Time-temperature correlation")+ylab("")+ggtitle("c - Evolution of phenotypic plasticity")+
geom_point(data=liste,aes(x=correlation,y=nbit),size=1,col="black")+coord_cartesian(expand=F)

cowplot::plot_grid(pl1,pl2,pl3,nrow=3,align = "v")

setwd(dir="C:/Users/Duchenne/Documents/plast_adaptation/")
pdf("figure_2.pdf",width=11,height=8)
cowplot::plot_grid(pl1,pl2,pl3,nrow=3,align = "v",rel_heights=c(1,1,1.1))
dev.off();
