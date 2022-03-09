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
library(spaMM)


t_init=12 #initial annual mean temperature
mfdi=180 #initial MFD
omega=5 #width of the fitness function
omegab=10 #width of the fitness function on plasticity
bruit_pheno=7

vec_z=c(0.03,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.4,0.5) #standard deviation of the gaussian used to draw temperature
#définir le nombre de chaine
nb_chaine=2000
nbit_vec=c(25,30,35,40,45,50,55)
nbdata_vec=c(500,1000,1500,2500)

ri=2 #the degree of the function defined for the increase of the record number across year. 1 = linear increase, 2 = quadratic increase...
#fonction pour définir le nombre de données par an à partir de la fonction :
func=function(x){sum(ceiling((c(1:nbit)/x)^ri))-nb_data}


#fitness function:
fitness=function(x){exp(-(x-mfd_opt)^2/(2*omega))}

#div: fait varier le nombre de données
#nb: fait varier le nombre d'années
#z fait varier la variance de la normale dans laquelle on tire la température -> faire varier la corrélation temps-température
#k répète le processus stochastique autant de fois qu'il y a de chaine
#i fait avance la chaine d'année en année
div=4
nb=3
z=6
k=1

wtd.mean=function(x,y){sum(x*y)/sum(y)}

for (div in 1:length(nbdata_vec)){
res2=as.data.frame(matrix(NA,length(nbit_vec)*nb_chaine*length(vec_z),33))
names(res2)=c("nbit","sigma","chaine","correlation","data_number","ti",
"t_adapt","t_plast","t_evolplast",
"adapt","plast","evol_plast","error_adapt","error_plast","error_inter",
"pval_plast","pval_adapt","pval_inter",
"adapt_lmer","plast_lmer","evol_plast_lmer","error_adapt_lmer","error_plast_lmer","error_inter_lmer",
"pval_plast_lmer","pval_adapt_lmer","pval_inter_lmer",
"pval_t_adapt","pval_t_evol_plast","annee_eval","temp_eval","G_adapt","G_plast")

nb_data=nbdata_vec[div] #number of records

for (nb in 1:length(nbit_vec)){
nbit=nbit_vec[nb] #number of years

#spread records by years :
# vec_nb=ceiling((c(1:nbit)/(uniroot(func,c(0,1000))[[1]]))^ri)
vec_nb=rep(floor(nb_data/nbit),nbit)
vec_nb[length(vec_nb)]=nb_data-sum(vec_nb[1:(length(vec_nb)-1)])

for (z in 1:length(vec_z)){
sigmat=vec_z[z] #standard deviation of the gaussian used to draw temperature



for (k in 1:nb_chaine){
res=as.data.frame(matrix(NA,sum(vec_nb),12))
names(res)=c("records","mfd_opt","temp","index","temp_variance","chaine",
"adaptation","mfd_adapt","evol_plast","plast","delta_temp","mfd_true")

b=0
adapt=runif(1,-6,0) #changes in the optimal flight date in function of the temperature (day/°C)
plast=runif(1,-6,6) #initial plasticity slope (day/°C)
mfd=mfdi #initial MFD (julian days)
G=matrix(c(runif(1,1,5),0,0,runif(1,1.5,2)),2,2) #covariance matrix of the reaction norm parameters
ti=0.015
for (i in 1:nbit){
a=b+1
#à chaque itération la température est tirée dans une gaussienne dont la moyenne augmente de 0.01 par an
tep=rnorm(1,t_init+i*ti,sigmat)
delta=tep-t_init
# la nouvelle date moyenne de vol optimale est définie en fonction de la température, par la valeur adapt:
mfd_opt=mfdi+delta*(adapt)
#pheno moyen:
xbar=c(mfd,plast)
zbar=c(1,delta)%*%xbar #MFD moyenne (MFd du pas d'avant + plasticité)
vec_c=matrix(c(1,delta),2,1) #vec_c = c(1,1) car on considère la variance de la plasticité constante à travers l'environnement
sigma=t(vec_c)%*%G%*%vec_c+bruit_pheno #variance intercept + variance de plasticité + bruit
#simuler les données pour l'année i:
veco=rnorm(vec_nb[i],zbar,sigma)
#couts de la plast impact xbar et G, si evol de la plast:
# L=matrix(c(0,0,0,1),2,2)
# xbar=xbar-1/(omegab+G[2,2])*G%*%L%*%xbar
# G=G-1/(omegab+G[2,2])*G%*%L%*%G
#selection:
betar=c(-1/(omega+sigma))*c(zbar-mfd_opt)*c(1,delta)
deltamu=G%*%betar
#2 cas soit on a des données, soit il y a zéro données cette année la:
if (length(veco)>0){
b=a+length(veco)-1
res[a:b,1]=veco
res[a:b,2]=mfd_opt
res[a:b,3]=tep
res[a:b,4]=i
res[a:b,5]=vec_z[z]
res[a:b,6]=k
res[a:b,7]=deltamu[1,]
res[a:b,8]=mfd
res[a:b,9]=deltamu[2,]
res[a:b,10]=plast
res[a:b,11]=delta
res[a:b,12]=zbar[1,1]
}else{
res[a,1]=NA
res[a,2]=mfd_opt
res[a,3]=tep
res[a,4]=i
res[a,5]=vec_z[z]
res[a,6]=k
res[a,7]=deltamu[1,]
res[a,8]=mfd
res[a,9]=deltamu[2,]
res[a,10]=plast
res[a,11]=delta
res[a,12]=zbar
b=b+1
}
mfd=mfd+deltamu[1,]
plast=plast+deltamu[2,] #si evol de la plast
#evolution de G
#G=G-G%*%vec_c%*%t(vec_c)%*%G/(omega+c(sigma)) #utiliser c(1,0) au lieu de vec_c si on veut changer que Gaa
}

bibi=res %>% dplyr::group_by(index,mfd_adapt,mfd_true,adaptation,plast,evol_plast,delta_temp,
temp) %>% dplyr::summarise(n=length(records),moy=mean(records))

annee_eval=mean(res$index)
temp_eval=mean(res$temp)
res$index=res$index-mean(bibi$index)
res$temp=res$temp-mean(bibi$temp)
model=lm(records~temp*index,data=res)
suma=summary(model)

model2=lm(records~temp+index,data=res)
suma2=summary(model2)
# mean(bibi$adaptation)
# bibi$plast[bibi$index==round(annee_eval)]
# mean(bibi$evol_plast)#evol_plast

# plot(adaptation~index,data=bibi)
# abline(summary(model)$coeff[3,1],0,col="blue")
# abline(summary(model)$coeff[3,1]+max(res$temp)*summary(model)$coeff[4,1],0,col="blue",lty="dashed")
# abline(summary(model)$coeff[3,1]+min(res$temp)*summary(model)$coeff[4,1],0,col="blue",lty="dashed")
# abline(mean(bibi$adaptation),0,col="red")
# abline(summary(model2)$coeff[3,1],0,col="orange",lty="dashed")

ano=Anova(model)
res$temp2=0
res$temp2[res$index>0]=res$temp[res$index>0]
ano2=Anova(model2)

# aval=suma2$coeff["index","Estimate"]+(bibi$index-mean(bibi$index))*(bibi$temp-mean(bibi$temp))*
# suma2$coeff["temp:index","Estimate"]
# mean(aval)
# mean(bibi$adaptation)

#le nombre d'année
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),1]=nbit
#l'écart type de la gaussienne utilisée pour la température
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),2]=vec_z[z]
# le numéro de la chaine
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),3]=k
# la corrélation temps-température
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),4]=cor(res$index,res$temp)
#le nombre de données totale
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),5]=sum(vec_nb)
#temperature increase
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),6]=ti

#les valeurs d'adaptation et de plasticité
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),7]=mean(bibi$adaptation)
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),8]=bibi$plast[bibi$index==round(annee_eval)]
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),9]=mean(bibi$evol_plast)#evol_plast

#estimation de l'adaptation
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),10]=suma$coeff["index","Estimate"]
#estimation de la plasticité
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),11]=suma$coeff["temp","Estimate"]
#estimation de l'interaction
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),12]=if("temp:index" %in% rownames(suma$coeff)){
suma$coeff["temp:index","Estimate"]}else{NA}
#les erreurs faites sur les estimation d'adaptation et de plasticité
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),13]=suma$coeff["index","Std. Error"]
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),14]=suma$coeff["temp","Std. Error"]
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),15]=if("temp:index" %in% rownames(suma$coeff)){
suma$coeff["temp:index","Std. Error"]}else{NA}
#les pvalues d'adaptation et de plasticité
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),16:18]=t(ano[1:3,4])

#estimation de l'adaptation lmer
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),19]=suma2$coeff["index","Estimate"]
#estimation de la plasticité lmer
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),20]=suma2$coeff["temp","Estimate"]
#estimation de l'interaction
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),21]=if("temp:index" %in% rownames(suma2$coeff)){
suma2$coeff["temp2","Estimate"]}else{NA}
#les erreurs faites sur les estimation d'adaptation et de plasticité lmer
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),22]=suma2$coeff["index","Std. Error"]
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),23]=suma2$coeff["temp","Std. Error"]
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),24]=if("temp:index2" %in% rownames(suma2$coeff)){
suma2$coeff["temp:index2","Std. Error"]}else{NA}
#les pvalues d'adaptation et de plasticité
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),25:27]=t(ano2[1:3,3])

#pval t_adapt
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),28]=summary(lm(mfd_adapt~index,data=bibi))$coeff[2,4]
#pval t_plast
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),29]=summary(lm(plast~index,data=bibi))$coeff[2,4]

#annee_eval
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),30]=annee_eval
#temp_eval
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),31]=temp_eval
#G
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),32]=G[1,1]
res2[k+(z-1)*nb_chaine+nb_chaine*length(vec_z)*(nb-1),33]=G[2,2]
}

# print(ggplot()+geom_point(data=res,aes(x=index,y=records),alpha=0.4)+geom_point(data=res,aes(x=index,y=mfd_true),col="firebrick",size=1.4)+
# geom_point(data=res,aes(x=index,y=mfd_adapt),col="dodgerblue",size=1.4)+geom_pointrange(data=newres,aes(x=index,y=fit,ymin=lwr,ymax=upr),col="green")+
# geom_abline(intercept=model$coef[1],slope=model$coef[3]))
# print(paste("data_number = ",nbdata_vec[div]," - number of years = ",nbit," - standard deviation = ",vec_z[z],sep="")) 
}
}
setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/résultats simus")
write.table(res2,paste("resultats_simulations_sans_g_evol_avec_plast_evol",nb_data,".txt",sep=""),row.names=F,sep="\t")
}

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
library(spaMM)

#FIGURE PLOT SURFACE DES SUPPLEMENTARY:
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


res2b=subset(res2,correlation<0.8 & correlation>0.2 & nbit>20)
pl1=ggplot(data=res2b,aes(x=t_adapt,y=adapt,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("a")+scale_color_viridis()+ylab("Estimation of evolution on\nthe reaction norm elevation (days/year)")+
xlab("True value (from simulation) of evolution on the reaction norm elevation (days/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl2=ggplot(data=res2b,aes(x=t_plast,y=plast,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("b")+scale_color_viridis()+ylab("Estimation of\nphenotypic plasticity (days/°C)")+
xlab("True value (from simulation) of phenotypic plasticity (days/°C)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

#grid.arrange(pl1,pl2)

png("figure_1.png",width=1100,height=700,res=120)
grid.arrange(pl1,pl2)
dev.off();



#png("figure supp.png",width=1100, height=800,res=130)
res2$correlation=round_any(res2$correlation,0.02)
b=summaryBy(diffea~correlation+nbit+data_number,data=subset(res2,correlation>0.2 & correlation<1),
FUN=quantile,probs=0.95,keep.names=T)
b$diffea[b$diffea>=0.4]=0.4
pl1=ggplot(data=b,aes(x=correlation,y=nbit,fill=diffea,col=diffea))+
geom_raster(interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),axis.title.x=element_blank(),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colours=rev(c(rainbow(10)[1:8],"darkblue")),na.value="white",limits=c(0,max(b$diffea)))+
labs(fill="Error in\nday/year")+ylab("")+ggtitle("a")

b=summaryBy(diffep~correlation+nbit+data_number,data=subset(res2,correlation>0.2 & correlation<1),
FUN=quantile,probs=0.95,keep.names=T)
b$diffep[b$diffep>=10]=10
pl2=ggplot(data=b,aes(x=correlation,y=nbit,fill=diffep,col=diffep))+
geom_raster(interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.x=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colours=rev(c(rainbow(10)[1:8],"darkblue")),na.value="white",limits=c(0,max(b$diffep)))+
labs(fill="Error in\nday/°C")+ylab("Number of years")+ggtitle("b")


b=summaryBy(diffei~correlation+nbit+data_number,data=subset(res2,correlation>0.2 & correlation<1),
FUN=quantile,probs=0.95,keep.names=T)
b$diffei[b$diffei>=1]=1
pl3=ggplot(data=b,aes(x=correlation,y=nbit,fill=diffei,col=diffei))+
geom_raster(interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),,
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colours=rev(c(rainbow(10)[1:8],"darkblue")),na.value="white",limits=c(0,max(b$diffei)))+
labs(fill="Error in\nday/°C/year")+xlab("Time-temperature correlation")+ylab("")+ggtitle("c")
cowplot::plot_grid(pl1,pl2,pl3,nrow=3,align = "v")

pdf("figure_S1.pdf",width=11,height=8)
cowplot::plot_grid(pl1,pl2,pl3,nrow=3,align = "v")
dev.off();


setwd(dir="C:/Users/Francois/Documents/papier 2 - plasticité- adaptation/résultats simus")
library(plyr)
res2b=fread("resultats_simulations_sans_g_evol_avec_plast_evol500.txt",header=T,sep="\t")
res2c=fread("resultats_simulations_sans_g_evol_avec_plast_evol1000.txt",header=T,sep="\t")
res2c2=fread("resultats_simulations_sans_g_evol_avec_plast_evol1500.txt",header=T,sep="\t")
res2d=fread("resultats_simulations_sans_g_evol_avec_plast_evol2500.txt",header=T,sep="\t")
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
ggtitle("a")+scale_color_viridis()+ylab("Estimation of evolution on\nthe reaction norm elevation (days/year)")+
xlab("True value (from simulation) of evolution on the reaction norm elevation (days/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl2=ggplot(data=res2b,aes(x=t_plast,y=plast,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("b")+scale_color_viridis()+ylab("Estimation of\nphenotypic plasticity (days/°C)")+
xlab("True value (from simulation) of phenotypic plasticity (days/°C)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl3=ggplot(data=res2b,aes(x=t_evolplast,y=evol_plast,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("c")+scale_color_viridis()+ylab("Estimation of evolution on\nphenotypic plasticity (days/°C/year)")+
xlab("True value (from simulation) of evolution on phenotypic plasticity (days/°C/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")


#grid.arrange(pl1,pl2,pl3)

png("figure_S2.png",width=1100,height=1100,res=120)
grid.arrange(pl1,pl2,pl3)
dev.off();

pl1=ggplot(data=res2b,aes(x=t_adapt,y=adapt_lmer,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("a")+scale_color_viridis()+ylab("Estimation of evolution on\nthe reaction norm elevation (days/year)")+
xlab("True value (from simulation) of evolution on the reaction norm elevation (days/year)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")

pl2=ggplot(data=res2b,aes(x=t_plast,y=plast_lmer,col=nbit))+
geom_point(size=0.1,alpha=0.2)+theme_bw()+stat_smooth(method="lm",col="darkgrey",size=1.3)+
geom_abline(intercept=0,slope=1,col="black",linetype="dashed")+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
ggtitle("b")+scale_color_viridis()+ylab("Estimation of\nphenotypic plasticity (days/°C)")+
xlab("True value (from simulation) of phenotypic plasticity (days/°C)")+
geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+
labs(col="Nb. of years")
#grid.arrange(pl1,pl2)

png("figure_S3.png",width=1100,height=700,res=120)
grid.arrange(pl1,pl2)
dev.off();


b=summaryBy(diffei~correlation+nbit+data_number,data=subset(res2,correlation>0.2),
FUN=quantile,probs=0.95,keep.names=T)
b$diffei[b$diffei>=1]=1
pl3=ggplot(data=b,aes(x=correlation,y=nbit,fill=diffei,col=diffei))+
geom_raster(interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),
legend.position="right",strip.text=element_text(size=12),,
plot.title=element_text(size=14,face="bold"))+facet_wrap(~data_number,ncol=4)+
scale_fill_gradientn(colours=rev(c(rainbow(10)[1:8],"darkblue")),limits=c(0,max(b$diffei)))+
labs(fill="Error in\nday/°C/year")+xlab("Time-temperature correlation")+ylab("")+ggtitle("c")

#en fonction de température-temps corrélation
pl1=ggplot(data=subset(res2, nbit==55),aes(x=correlation,y=(diffea)))+geom_point(alpha=0.3)+geom_vline(xintercept=0.76,color="red")+
theme_bw()+theme(axis.title.x=element_blank())+ylab("error in adaptation estimation (%)")+facet_wrap(~data_number)+stat_smooth()
pl2=ggplot(data=subset(res2, nbit==55),aes(x=correlation,y=(diffep)))+geom_point(alpha=0.3)+geom_vline(xintercept=0.76,color="red")+
theme_bw()+theme(axis.title.x=element_blank())+ylab("error in plastic response estimation (%)")+facet_wrap(~data_number)+stat_smooth()
grid.arrange(pl1,pl2,nrow=2,bottom="time-temperature correlation coefficient")
#en fonction de t_adapt
pl1=ggplot(data=subset(res2,nbit==55 & correlation<0.8),aes(x=t_evolplast,y=(diffea)))+geom_point(alpha=0.3)+
theme_bw()+theme(axis.title.x=element_blank())+ylab("error in adaptation estimation (%)")+
facet_wrap(~data_number,scales="free_y")+stat_smooth()
pl2=ggplot(data=subset(res2,nbit==55 & correlation<0.8),aes(x=t_evolplast,y=(diffep)))+geom_point(alpha=0.3)+
theme_bw()+theme(axis.title.x=element_blank())+ylab("error in plastic response estimation (%)")+
facet_wrap(~data_number,scales="free_y")+stat_smooth()
grid.arrange(pl1,pl2,nrow=2,bottom="time-temperature correlation coefficient")
#density
pl1=ggplot(data=subset(res2,correlation<0.9),aes(x=(diffea),col=as.factor(nbit)))+geom_density(alpha=0.3)+
theme_bw()+theme(axis.title.x=element_blank())+xlab("error in adaptation estimation")+facet_wrap(~data_number)+
geom_vline(xintercept=0)
pl2=ggplot(data=subset(res2,correlation<0.9),aes(x=(diffep),col=as.factor(nbit)))+geom_density(alpha=0.3)+
theme_bw()+theme(axis.title.x=element_blank())+
xlab("error in plastic response estimation")+facet_wrap(~data_number)+xlim(c(-20,20))+geom_vline(xintercept=0)
grid.arrange(pl1,pl2,nrow=2)
X11()
pl1=ggplot(data=subset(res2,correlation<0.9),aes(x=(diffea2),col=as.factor(nbit)))+geom_density(alpha=0.3)+
theme_bw()+theme(axis.title.x=element_blank())+xlab("error in adaptation estimation")+facet_wrap(~data_number)
pl2=ggplot(data=subset(res2,correlation<0.9),aes(x=(diffep2),col=as.factor(nbit)))+geom_density(alpha=0.3)+
theme_bw()+theme(axis.title.x=element_blank())+xlab("error in plastic response estimation")+facet_wrap(~data_number)
grid.arrange(pl1,pl2,nrow=2)



ggplot(data=subset(res2,correlation<0.9),aes(x=diffea2,y=diffea,col=nbit))+geom_point()+stat_smooth(method="lm")+
theme_bw()+theme(axis.title.x=element_blank())+facet_wrap(~data_number)+geom_abline(intercept=0,slope=1)

res2$lwr_adapt=res2$adapt-qnorm(1-1/(2*res2$data_number))*res2$error_adapt
res2$upr_adapt=res2$adapt+qnorm(1-1/(2*res2$data_number))*res2$error_adapt
res2$lwr_plast=res2$plast-qnorm(1-1/(2*res2$data_number))*res2$error_plast
res2$upr_plast=res2$plast+qnorm(1-1/(2*res2$data_number))*res2$error_plast
res2$correlationb=plyr::round_any(res2$correlation,0.05)
res2$data_number=plyr::round_any(res2$data_number,100)
res2$tib=round(res2$ti,digits=3)
#faux positifs:
res2$ok_adapt=1 #erreur de type I (faux positifs)
res2$ok_adapt[which(res2$pval_adapt<0.05 & res2$pval_t_adapt<0.05 & res2$t_adapt>0 & res2$adapt>0)]=0 #vrais positifs
res2$ok_adapt[which(res2$pval_adapt<0.05 & res2$pval_t_adapt<0.05 & res2$t_adapt<0 & res2$adapt<0)]=0 #vrais positifs
res2$ok_adapt[which(res2$pval_adapt>0.05 & res2$pval_t_adapt>0.05)]=0 #vrais négatifs
res2$ok_adapt[which(res2$pval_adapt>0.05 & res2$pval_t_adapt<0.05)]=0 #faux négatifs
b=summaryBy(ok_adapt~nbit+data_number+correlationb,data=res2,FUN=c(mean,length),keep.names=T)
ggplot(data=subset(b,ok_adapt.length>20 & correlationb>0),aes(x=correlationb,y=ok_adapt.mean,col=as.factor(nbit)))+geom_line()+facet_wrap(~data_number)+theme_bw()+
ylab("Type I error (false positive)")+xlab("Time-Temperature correlation")
ggplot(data=subset(b,ok_adapt.length>20 & correlationb>0),aes(x=data_number,y=ok_adapt.mean,col=correlationb,group=correlationb))+geom_line()+facet_wrap(~nbit,ncol=5)+theme_bw()+
ylab("Type I error (false positive)")+xlab("Number of records")+theme(legend.position="bottom")+labs(color="Time-Temperature correlation")+
scale_colour_gradientn(colours = viridis(20))
#puissance:
res2$ok_adapt=NA #erreur de type I (faux positifs)
res2$ok_adapt[which(res2$pval_adapt<0.05 & res2$pval_t_adapt<0.05 & res2$t_adapt>0 & res2$adapt>0)]=1 #vrais positifs
res2$ok_adapt[which(res2$pval_adapt<0.05 & res2$pval_t_adapt<0.05 & res2$t_adapt<0 & res2$adapt<0)]=1 #vrais positifs
res2$ok_adapt[which(res2$pval_adapt>0.05 & res2$pval_t_adapt>0.05)]=NA #vrais négatifs
res2$ok_adapt[which(res2$pval_adapt>0.05 & res2$pval_t_adapt<0.05)]=0 #faux négatifs
b=summaryBy(ok_adapt~data_number+correlationb+nbit,data=subset(res2,!is.na(ok_adapt)),FUN=c(mean,length))
ggplot(data=subset(b,ok_adapt.length>3 & correlationb>0),aes(x=correlationb,y=ok_adapt.mean,col=as.factor(tib)))+geom_line()+facet_wrap(~data_number)+theme_bw()+
ylab("Type II error (false negative)")+xlab("Time-Temperature correlation")
ggplot(data=subset(b,ok_adapt.length>10 & correlationb>0),aes(x=data_number,y=ok_adapt.mean,col=correlationb,group=correlationb))+geom_line()+facet_wrap(~nbit,ncol=5)+theme_bw()+
ylab("Power (true positives detected)")+xlab("Number of records")+theme(legend.position="bottom")+labs(color="Time-Temperature correlation")+
scale_colour_gradientn(colours = viridis(20))
#estimate:
ggplot(data=res2,aes(x=correlationb,y=abs(diffea),col=as.factor(nbit),fill=as.factor(nbit)))+stat_smooth(alpha=0.4)+facet_wrap(~data_number)+theme_bw()+
ylab("Estimation error")+xlab("Time-Temperature correlation")





#importance
res2$imp_adapt=abs(res2$importance_adapt)/apply(abs(res2[,c("importance_adapt","importance_plast")]),1,sum)
res2$t_imp_adapt=abs(res2$t_importance_adapt)/apply(abs(res2[,c("t_importance_plast","t_importance_adapt")]),1,sum)
res2$diffimp=res2$t_imp_adapt-res2$imp_adapt
b=summaryBy(diffimp~nbit+data_number+correlationb,data=res2,FUN=mean)
ggplot(data=subset(res2,nbit==50),aes(x=as.factor(correlationb),y=abs(diffimp)))+geom_boxplot()+facet_wrap(~data_number)+theme_bw()

res2$data_number=round_any(res2$data_number,100)
res2$datan=paste("n=",res2$data_number,sep="")
res2$datan=factor(res2$datan,c("n=100","n=500","n=1500"))
pl1=ggplot(data=subset(res2,correlation>0.45 & correlation<0.55 & nbit!=20 & nbit!=50),aes(x=plasticite,group=as.factor(nbit),fill=as.factor(nbit)))+geom_density(alpha=0.5)+
theme_bw()+geom_vline(xintercept=-0.09)+theme(legend.position="none")+xlab("Adaptation estimation")+labs(fill="Number of years")+facet_wrap(~datan,scales="free")
pl2=ggplot(data=subset(res2,correlation>0.45 & correlation<0.55 & nbit!=20 & nbit!=50),aes(x=adapt,group=as.factor(nbit),fill=as.factor(nbit)))+geom_density(alpha=0.5)+facet_wrap(~datan,scales="free")+
theme_bw()+geom_vline(xintercept=-6.5)+labs(fill="Number of years")+xlab("Plasticity estimation")+theme(legend.position="bottom")
png("figure supp2.png",width=1100, height=800,res=120)
grid.arrange(pl1,pl2,nrow=2,heights=c(1.8,2))
dev.off ()
gA <- ggplotGrob(pl1)
gB <- ggplotGrob(pl2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))
group.CI(adapt~datan+nbit,data=res2)


model=randomForest(abs(diffep)~nbit+data_number+tib+correlation+t_adapt+t_plast,data=res2,importance=T)

res2$correlation2=round_any(res2$correlation,0.05)
veci=unique(res2$correlation2[which(res2$correlation2<0.75)])
par(mfrow=c(1,2))
boxplot(abs(diffea)~correlation2,data=res2,ylimits=c(0,100),border=c(rep("black",length(veci)),"red",rep("black",4)))
boxplot(abs(res2$diffep)~res2$correlation2,ylimits=c(0,100),border=c(rep("black",length(veci)),"red",rep("black",4)))

#FIGURE PLOT SURFACE DES SUPPLEMENTARY:
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
res2$data_number=round_any(res2$data_number,100)
#png("figure supp.png",width=1100, height=800,res=130)
res2$correlation=round_any(res2$correlation,0.02)
b=summaryBy(diffea~correlation+nbit+data_number,data=subset(res2,correlation>0.2),
FUN=quantile,probs=0.95,keep.names=T)
pl1=ggplot(data=b,aes(x=correlation,y=nbit,fill=diffea,col=diffea))+
geom_raster(interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.y=element_blank(),
legend.position="right",strip.text=element_text(size=12),,axis.title=element_blank(),
plot.subtitle=element_text(size=14))+facet_wrap(~data_number,ncol=4)+scale_fill_viridis(option = "plasma")

pl2=ggplot(data=b,aes(x=correlation,y=nbit,fill=diffep,col=diffep))+
geom_raster(interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.y=element_blank(),
legend.position="right",strip.text=element_text(size=12),,axis.title=element_blank(),
plot.subtitle=element_text(size=14))+facet_wrap(~data_number,ncol=4)+scale_fill_viridis(option = "plasma")

pl3=ggplot(data=b,aes(x=correlation,y=nbit,fill=diffei,col=diffei))+
geom_raster(interpolate=T)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.y=element_blank(),
legend.position="right",strip.text=element_text(size=12),,axis.title=element_blank(),
plot.subtitle=element_text(size=14))+facet_wrap(~data_number,ncol=4)+scale_fill_viridis(option = "plasma")

grid.arrange(pl1,pl2,pl3,nrow=3)


layout(matrix(c(1:12),nrow=3,byrow=T))
par(mar=c(2, 3, 2, 2))
res3=subset(res2,correlation>0.2 & data_number==500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffea~correlation+nbit+data_number,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffea")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>0.5]=0.5
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="Number of years")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.15)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.05)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.1)
res3=subset(res2,correlation>0.2 & data_number==1000)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffea~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffea")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>0.5]=0.5
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="Time-Temperature correlation",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.15)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.05)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.1)
res3=subset(res2,correlation>0.2 & data_number==1500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffea~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffea")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>0.5]=0.5
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.15)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.05)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.1)
res3=subset(res2,correlation>0.2 & data_number==2500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffea~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffea")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>0.5]=0.5
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.15)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.05)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.1)

#####PLASTICITE#########
res3=subset(res2,correlation>0.2 & data_number==500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffep~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffep")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="Number of years")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==1000)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffep~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffep")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="Time-Temperature correlation",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==1500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffep~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffep")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==2500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffep~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffep")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)



#####EVOL PLASTICITE#########
res3=subset(res2,correlation>0.2 & data_number==500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffei~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffei")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>1]=1
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="Number of years")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==1000)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffei~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffei")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>1]=1
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="Time-Temperature correlation",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==1500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffei~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffei")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>1]=1
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==2500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffei~correlation+nbit,data=res3,FUN=quantile,probs=0.95,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffei")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>1]=1
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)



dev.off ();




#FIGURE PLOT SURFACE DES SUPPLEMENTARY MOYENNES:
setwd(dir="C:/Users/Francois/Documents/papier 1 - données collection/résultats simus")
res2=read.table("resultats_simulations_100.txt",header=T,sep="\t")
res2b=read.table("resultats_simulations_500.txt",header=T,sep="\t")
res2c=read.table("resultats_simulations_1000.txt",header=T,sep="\t")
res2=rbind(res2,res2b,res2c)
res2$diffea=abs(res2$adapta-res2$t_adapt)
res2$diffep=abs(res2$plast-res2$t_plast)
res2$diffei=abs(res2$evol_plast-res2$t_evolplast)
res2$data_number=round_any(res2$data_number,100)
#png("figure supp.png",width=1100, height=800,res=130)
layout(matrix(c(1,4,7,2,5,8,3,6,9),nrow=3))
par(mar=c(2, 3, 2, 2))
res3=subset(res2,correlation>0.2 & data_number==100)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffea~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffea")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>0.5]=0.5
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="Number of years")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.15)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.05)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.1)
res3=subset(res2,correlation>0.2 & data_number==500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffea~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffea")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>0.5]=0.5
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="Time-Temperature correlation",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.15)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.05)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.1)
res3=subset(res2,correlation>0.2 & data_number==1000)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffea~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffea")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>0.5]=0.5
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.15)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.05)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=0.1)
#####PLASTICITE#########
res3=subset(res2,correlation>0.2 & data_number==100)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffep~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffep")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="Number of years")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffep~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffep")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="Time-Temperature correlation",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==1000)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffep~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffep")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
#####EVOL PLASTICITE#########
res3=subset(res2,correlation>0.2 & data_number==100)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffei~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffei")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>1]=1
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="Number of years")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==500)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffei~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffei")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>1]=1
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="Time-Temperature correlation",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
res3=subset(res2,correlation>0.2 & data_number==1000)
res3$correlation=round_any(res3$correlation,0.02)
b=summaryBy(diffei~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffei")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>1]=1
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,xlab="",ylab="")
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=2)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=4)
# contour2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)), col = "black", labcex=1,lwd=2, alpha=1,resfac=0.5, add = TRUE,levels=1)
dev.off ();



#FIGURE PLOT SURFACE DES SUPPLEMENTARY EN BLANC:
res2=read.table("resultats_simulations_100.txt",header=T,sep="\t")
res2b=read.table("resultats_simulations_500.txt",header=T,sep="\t")
res2c=read.table("resultats_simulations_1500.txt",header=T,sep="\t")
res2=rbind(res2,res2b,res2c)
res2$diffea=abs(res2$adapt-res2$t_adapt)-1.96*res2$error_adapt
res2$diffep=abs(res2$plast-res2$t_plast)-1.96*res2$error_plast
res2$diffeab=0
res2$diffeab[which(res2$diffea<=0)]=1
res2$diffepb=0
res2$diffepb[which(res2$diffep<=0)]=1
res3=subset(res2,correlation>0.2)
res3$correlation=round_any(res3$correlation,0.02)
layout(matrix(c(1,2),nrow=2),heights=c(10,10))
par(mar=c(4,4,1,1.5))
b=summaryBy(diffeab~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffeab")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,ylab="Number of years",xlab="")
b=summaryBy(diffepb~correlation+nbit,data=res3,FUN=mean,keep.names=T)
volc=dcast(b,correlation~nbit,value.var="diffepb")
rownames(volc)=volc$correlation
volc=volc[,-1]
volc=as.matrix(volc)
volc[volc>10]=10
image2D(volc,y=as.numeric(colnames(volc)),x=as.numeric(rownames(volc)),resfac=10,ylab="Number of years",xlab="Time-Temperature correlation")



par(mfrow=c(1,2))
#boxplot(res2$diffe~res2$sigma)
plot(res2$diffea~res2$correlation,pch=".",ylim=c(-10,10),cex=2,xlab="corrélation temperature-temps",ylab="erreur d'estimation de l'adaptation en %")
summaryBy(diffea~sigma,data=res2,FUN=mean)
abline(0,0,col="red")
abline(v=0.77,col="green")
summary((subset(res2,sigma>0.2)$diffea/subset(res2,sigma>0.2)$adaptation)*100)
#boxplot(res2$diffe~res2$sigma)
plot(res2$diffep~res2$correlation,pch=".",ylim=c(-10,10),cex=2,xlab="corrélation temperature-temps",ylab="erreur d'estimation de la plasticité en %")
summaryBy(diffep~sigma,data=res2,FUN=mean)
summary((subset(res2,sigma>0.2)$diffep/subset(res2,sigma>0.2)$plasticite)*100)
abline(0,0,col="red")
abline(v=0.77,col="green")


par(mfrow=c(1,2))
res2$diffe=res2$adapt-res2$adaptation
plot(res2$diffe~res2$correlation)
abline(0,0,col="red")
res2$diffe=res2$plast2-res2$plasticite
plot(res2$diffe~res2$correlation)
abline(0,0,col="red")
