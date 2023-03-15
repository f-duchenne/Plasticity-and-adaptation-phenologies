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

#FIGURE PLOT SURFACE DES SUPPLEMENTARY:
setwd(dir="C:/Users/duchenne/Documents/plast_adaptation/data/resultats_simues")
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

############# FIGURE 3
res2b=subset(res2,correlation<0.8 & correlation>0.2 & nbit>20)
pl1=ggplot(data=res2b,aes(x=t_adapt,y=adapt,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("a")+scale_color_viridis(n.breaks =4)+ylab("Estimation of evolution on\nthe reaction norm elevation (days/year)")+
xlab("True value (from simulation) of evolution on the reaction norm elevation (days/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl2=ggplot(data=res2b,aes(x=t_plast,y=plast,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("b")+scale_color_viridis(n.breaks =4)+ylab("Estimation of\nphenotypic plasticity (days/°C)")+
xlab("True value (from simulation) of phenotypic plasticity (days/°C)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

#grid.arrange(pl1,pl2)

png("figure_3.png",width=1100,height=700,res=120)
grid.arrange(pl1,pl2)
dev.off();

############# FIGURE S1
setwd(dir="C:/Users/duchenne/Documents/plast_adaptation/data/resultats_simues")
library(plyr)
res2b=fread("resultats_simulations_sans_g_evol_avec_plast_evol500.txt",header=T,sep="\t",fill = TRUE)
res2c=fread("resultats_simulations_sans_g_evol_avec_plast_evol1000.txt",header=T,sep="\t",fill = TRUE)
res2c2=fread("resultats_simulations_sans_g_evol_avec_plast_evol1500.txt",header=T,sep="\t",fill = TRUE)
res2d=fread("resultats_simulations_sans_g_evol_avec_plast_evol2500.txt",header=T,sep="\t",fill = TRUE)
res2=rbind(res2d,res2b,res2c,res2c2)
res2$diffea=abs(res2$adapt-res2$t_adapt)
res2$diffep=abs(res2$plast-res2$t_plast)
res2$diffei=abs(res2$evol_plast-res2$t_evolplast)
res2$data_number=paste("n = ",round_any(res2$data_number,100),sep="")
res2$data_number=factor(res2$data_number,c("n = 500","n = 1000","n = 1500","n = 2500"))
#png("figure supp.png",width=1100, height=800,res=130)
res2$correlation=round_any(res2$correlation,0.02)

res2b=subset(res2,correlation<0.8 & correlation>0.2 & nbit>20)
pl1=ggplot(data=res2b,aes(x=t_adapt,y=adapt,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("a")+scale_color_viridis(n.breaks =4)+ylab("Estimation of evolution on\nthe reaction norm elevation (days/year)")+
xlab("True value (from simulation) of evolution on the reaction norm elevation (days/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl2=ggplot(data=res2b,aes(x=t_plast,y=plast,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("b")+scale_color_viridis(n.breaks =4)+ylab("Estimation of\nphenotypic plasticity (days/°C)")+
xlab("True value (from simulation) of phenotypic plasticity (days/°C)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl3=ggplot(data=res2b,aes(x=t_evolplast,y=evol_plast,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("c")+scale_color_viridis(n.breaks =4)+ylab("Estimation of evolution on\nphenotypic plasticity (days/°C/year)")+
xlab("True value (from simulation) of evolution on phenotypic plasticity (days/°C/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

png("figure_S1.png",width=1100,height=1100,res=120)
grid.arrange(pl1,pl2,pl3)
dev.off();


############# FIGURE S2
pl1=ggplot(data=res2b,aes(x=t_adapt,y=adapt_lmer,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("a")+scale_color_viridis(n.breaks =4)+ylab("Estimation of evolution on\nthe reaction norm elevation (days/year)")+
xlab("True value (from simulation) of evolution on the reaction norm elevation (days/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl2=ggplot(data=res2b,aes(x=t_plast,y=plast_lmer,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("b")+scale_color_viridis(n.breaks =4)+ylab("Estimation of\nphenotypic plasticity (days/°C)")+
xlab("True value (from simulation) of phenotypic plasticity (days/°C)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")
#grid.arrange(pl1,pl2)

png("figure_S2.png",width=1100,height=700,res=120)
grid.arrange(pl1,pl2)
dev.off();


############# FIGURE LMER
res2_lm=res2
res2b=fread("resultats_simulations_sans_g_evol_avec_plast_evol_lmer_500.txt",header=T,sep="\t")
res2c=fread("resultats_simulations_sans_g_evol_avec_plast_evol_lmer_1000.txt",header=T,sep="\t")
res2c2=fread("resultats_simulations_sans_g_evol_avec_plast_evol_lmer_1500.txt",header=T,sep="\t")
res2d=fread("resultats_simulations_sans_g_evol_avec_plast_evol_lmer_2500.txt",header=T,sep="\t")
res2_lmer=rbind(res2d,res2b,res2c,res2c2)
res2_lmer$diffea=abs(res2_lmer$adapt-res2_lmer$t_adapt)
res2_lmer$diffep=abs(res2_lmer$plast-res2_lmer$t_plast)
res2_lmer$diffei=abs(res2_lmer$evol_plast-res2_lmer$t_evolplast)
res2_lmer$data_number=paste("n = ",round_any(res2_lmer$data_number,100),sep="")
res2_lmer$data_number=factor(res2_lmer$data_number,c("n = 500","n = 1000","n = 1500","n = 2500"))
#png("figure supp.png",width=1100, height=800,res=130)
res2_lmer$correlation=round_any(res2_lmer$correlation,0.02)

b=res2_lmer %>% dplyr::group_by(nbit,data_number,correlation) %>% dplyr::summarise(diffea_lmer=mean(diffea),diffep_lmer=mean(diffep),diffei_lmer=mean(diffei))
b2=res2_lm %>% dplyr::group_by(nbit,data_number,correlation) %>% dplyr::summarise(diffea_lm=mean(diffea),diffep_lm=mean(diffep),diffei_lm=mean(diffei))

bf=merge(b,b2,by=c("nbit","data_number","correlation"))

ggplot(data=bf,aes(x=diffea_lm,y=diffea_lmer,col=nbit))+
geom_point(size=1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("b")+scale_color_viridis(n.breaks =4)+ylab("Estimation of\nphenotypic plasticity (days/°C)")+
xlab("True value (from simulation) of phenotypic plasticity (days/°C)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

