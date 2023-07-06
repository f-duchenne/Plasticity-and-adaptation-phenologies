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

res2$diffea=abs(res2$adapt-res2$t_adapt)
res2$diffep=abs(res2$plast-res2$t_plast)
res2$diffei=abs(res2$evol_plast-0)
res2$data_number=paste("n = ",round_any(res2$data_number,100),sep="")
res2$data_number=factor(res2$data_number,c("n = 500","n = 1000","n = 1500","n = 2500"))

############# FIGURE S1
res2b=subset(res2,correlation<0.8 & correlation>0.2 & nbit>20 & mechanism=="without_evol_plast")
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
setwd(dir="C:/Users/duchenne/Documents/plast_adaptation")
png("figure_S1.png",width=1100,height=700,res=120)
grid.arrange(pl1,pl2)
dev.off();

############# FIGURE S2
setwd(dir="C:/Users/duchenne/Documents/plast_adaptation/data/resultats_simues")
library(plyr)
res2b=subset(res2,correlation<0.8 & correlation>0.2 & nbit>20 & mechanism=="total")
res2b$correlation=round_any(res2b$correlation,0.02)

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

setwd(dir="C:/Users/duchenne/Documents/plast_adaptation")
png("figure_S2.png",width=1100,height=1100,res=120)
grid.arrange(pl1,pl2,pl3)
dev.off();


############# FIGURE S3
pl1=ggplot(data=res2b,aes(x=t_adapt,y=adapt_wi,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("a")+scale_color_viridis(n.breaks =4)+ylab("Estimation of evolution on\nthe reaction norm elevation (days/year)")+
xlab("True value (from simulation) of evolution on the reaction norm elevation (days/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl2=ggplot(data=res2b,aes(x=t_plast,y=plast_wi,col=nbit))+
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

setwd(dir="C:/Users/duchenne/Documents/plast_adaptation")
png("figure_S3.png",width=1100,height=700,res=120)
grid.arrange(pl1,pl2)
dev.off();


############# FIGURE S4
alpha=0.05
res2$cate=as.character(NA)
res2$cate[res2$t_adapt_pval<alpha & res2$pval_adapt<alpha]="true positive"
res2$cate[res2$t_adapt_pval>alpha & res2$pval_adapt>alpha]="true negative"
res2$cate[res2$t_adapt_pval>alpha & res2$pval_adapt<alpha]="false positive"
res2$cate[res2$t_adapt_pval>alpha & res2$pval_adapt<alpha & res2$adapt<0.01]="false positive, but very small estimate"
res2$cate[res2$t_adapt_pval<alpha & res2$pval_adapt>alpha & sign(res2$t_adapt)!=sign(res2$adapt)]="opposite"
res2$cate[res2$t_adapt_pval<alpha & res2$pval_adapt>alpha]="false negative"

b=subset(res2,abs(correlation)<0.7) %>% dplyr::group_by(mechanism,nbit,data_number,cate) %>% dplyr::count()
b$cate=factor(b$cate,levels=rev(c("false positive","false positive, but very small estimate","true positive","true negative","false negative")))

pl1=ggplot(data=subset(b,mechanism=="total"),aes(x=nbit,fill=cate,y=n))+geom_bar(stat="identity",position="fill")+facet_wrap(~data_number)+theme_bw()+
scale_y_continuous(labels=scales::percent)+
geom_hline(yintercept=0.05)+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+coord_cartesian(expand=F)+
scale_fill_manual(values = rev(c('red', 'pink',"chartreuse3","lightgreen", 'dodgerblue3')),
breaks=rev(c("false positive","false positive, but very small estimate","true positive","true negative","false negative")),name="")+
ggtitle("a - all mechanisms",subtitle=expression(paste("[",Gaa %~% U(1,5)," , ", Gbb %~% U(0.5,2)," , ",b[t==0] %~% U(-6,6),"]")))+
ylab("Percentage of simulations")+xlab("Number of years")

pl2=ggplot(data=subset(b,mechanism=="without_plast"),aes(x=nbit,fill=cate,y=n))+geom_bar(stat="identity",position="fill")+facet_wrap(~data_number)+theme_bw()+
scale_y_continuous(labels=scales::percent)+
geom_hline(yintercept=0.05)+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+coord_cartesian(expand=F)+
scale_fill_manual(values = rev(c('red', 'pink',"chartreuse3","lightgreen", 'dodgerblue3')),
breaks=rev(c("false positive","false positive, but very small estimate","true positive","true negative","false negative")),name="")+
ggtitle("b - without plasticity",subtitle=expression(paste("[",Gaa %~% U(1,5)," , ", Gbb == 0," , ",b[t==0] == 0,"]")))+
ylab("Percentage of simulations")+xlab("Number of years")

pl3=ggplot(data=subset(b,mechanism=="without_adapt"),aes(x=nbit,fill=cate,y=n))+geom_bar(stat="identity",position="fill")+facet_wrap(~data_number)+theme_bw()+
scale_y_continuous(labels=scales::percent)+
geom_hline(yintercept=0.05)+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+coord_cartesian(expand=F)+
scale_fill_manual(values = rev(c('red', 'pink',"chartreuse3","lightgreen", 'dodgerblue3')),
breaks=rev(c("false positive","false positive, but very small estimate","true positive","true negative","false negative")),name="")+
ggtitle("c - without evolution of elevation",subtitle=expression(paste("[",Gaa == 0," , ", Gbb %~% U(0.5,2)," , ",b[t==0] %~% U(-6,6),"]")))+
ylab("Percentage of simulations")+xlab("Number of years")

pl4=ggplot(data=subset(b,mechanism=="without_evol_plast"),aes(x=nbit,fill=cate,y=n))+geom_bar(stat="identity",position="fill")+facet_wrap(~data_number)+theme_bw()+
scale_y_continuous(labels=scales::percent)+
geom_hline(yintercept=0.05)+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+coord_cartesian(expand=F)+
scale_fill_manual(values = rev(c('red', 'pink',"chartreuse3","lightgreen", 'dodgerblue3')),
breaks=rev(c("false positive","false positive, but very small estimate","true positive","true negative","false negative")),name="")+
ggtitle("d - without evolution of plasticity",subtitle=expression(paste("[",Gaa %~% U(1,5)," , ", Gbb == 0," , ",b[t==0] %~% U(-6,6),"]")))+
ylab("Percentage of simulations")+xlab("Number of years")

pl5=ggplot(data=subset(b,mechanism=="without_evolution"),aes(x=nbit,fill=cate,y=n))+geom_bar(stat="identity",position="fill")+facet_wrap(~data_number)+theme_bw()+
scale_y_continuous(labels=scales::percent)+
geom_hline(yintercept=0.05)+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+coord_cartesian(expand=F)+
scale_fill_manual(values = rev(c('red', 'pink',"chartreuse3","lightgreen", 'dodgerblue3')),
breaks=rev(c("false positive","false positive, but very small estimate","true positive","true negative","false negative")),name="")+
ggtitle("e - without evolution",subtitle=expression(paste("[",Gaa == 0," , ", Gbb == 0," , ",b[t==0] %~% U(-6,6),"]")))+
ylab("Percentage of simulations")+xlab("Number of years")


leg <- ggpubr::as_ggplot(cowplot::get_legend(pl1))
pl1=pl1+theme(legend.position="none")
pl2=pl2+theme(legend.position="none")
pl3=pl3+theme(legend.position="none")
pl4=pl4+theme(legend.position="none")
pl5=pl5+theme(legend.position="none")


setwd(dir="C:/Users/duchenne/Documents/plast_adaptation")
png("figure_S4.png",width=1000,height=1200,res=120)
grid.arrange(pl1,pl2,pl3,pl4,pl5,leg,layout_matrix=rbind(c(1,2),c(3,4),c(5,56)),heights=c(1,1,1))
dev.off();

