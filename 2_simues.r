library(dplyr)
library(car)

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])
print(i)

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
mechanisms=c("without_adapt","without_plast","without_evol_plast","total","without_evolution")

tabtab=expand.grid(nbdata_vec,nbit_vec,mechanisms)

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


res2=as.data.frame(matrix(NA,nb_chaine*length(vec_z),37))
names(res2)=c("nbit","sigma","chaine","correlation","data_number","ti",
"t_adapt","t_plast","t_evolplast",
"t_adapt_pval","t_plast_pval","t_evolplast_pval",
"adapt","plast","evol_plast","error_adapt","error_plast","error_inter",
"pval_plast","pval_adapt","pval_inter",
"adapt_wi","plast_wi","evol_plast_wi","error_adapt_wi","error_plast_wi","error_inter_wi",
"pval_plast_wi","pval_adapt_wi","pval_inter_wi",
"pval_t_adapt","pval_t_evol_plast","annee_eval","temp_eval","G_adapt","G_plast","mechanism")

nb_data=tabtab[i,1] #number of records

nbit=tabtab[i,2]  #number of years

#spread records by years :
# vec_nb=ceiling((c(1:nbit)/(uniroot(func,c(0,1000))[[1]]))^ri)
vec_nb=rep(floor(nb_data/nbit),nbit)
vec_nb[length(vec_nb)]=nb_data-sum(vec_nb[1:(length(vec_nb)-1)])

mecha=tabtab[i,3]

for (z in 1:length(vec_z)){
sigmat=vec_z[z] #standard deviation of the gaussian used to draw temperature



for (k in 1:nb_chaine){
res=as.data.frame(matrix(NA,sum(vec_nb),12))
names(res)=c("records","mfd_opt","temp","index","temp_variance","chaine",
"adaptation","mfd_adapt","evol_plast","plast","delta_temp","mfd_true")

b=0
adapt=runif(1,-6,0) #changes in the optimal flight date in function of the temperature (day/°C)
plast=runif(1,-6,6) #initial plasticity slope (day/°C)
if(mecha=="without_plast"){plast=0}
mfd=mfdi #initial MFD (julian days)
if(mecha=="without_adapt"){G=matrix(c(0,0,0,runif(1,0.5,2)),2,2)} #covariance matrix of the reaction norm parameters
if(mecha=="without_evolution"){G=matrix(c(0,0,0,0),2,2)} #covariance matrix of the reaction norm parameters
if(mecha %in% c("without_plast","without_evol_plast")){G=matrix(c(runif(1,1,5),0,0,0),2,2)}
if(mecha=="total"){G=matrix(c(runif(1,1,5),0,0,runif(1,0.5,2)),2,2)}

ti=0.015 #temperature increase
for (i in 1:nbit){
a=b+1
#à chaque itération la température est tirée dans une gaussienne dont la moyenne augmente de ti par an
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

index_row=k+(z-1)*nb_chaine
#le nombre d'année
res2[index_row,"nbit"]=nbit
#l'écart type de la gaussienne utilisée pour la température
res2[index_row,"sigma"]=vec_z[z]
# le numéro de la chaine
res2[index_row,"chaine"]=k
# la corrélation temps-température
res2[index_row,"correlation"]=cor(res$index,res$temp)
#le nombre de données totale
res2[index_row,"data_number"]=sum(vec_nb)
#temperature increase
res2[index_row,"ti"]=ti

#les valeurs d'adaptation et de plasticité
res2[index_row,"t_adapt"]=mean(bibi$adaptation)
res2[index_row,"t_plast"]=bibi$plast[bibi$index==round(annee_eval)]
res2[index_row,"t_evolplast"]=mean(bibi$evol_plast)#evol_plast
#les erreurs d'adaptation et de plasticité
res2[index_row,"t_adapt_pval"]=if(mecha %in% c("without_adapt","without_evolution")){1}else{summary(lm(mfd_adapt~index,data=bibi))$coeff["index","Pr(>|t|)"]}
res2[index_row,"t_plast_pval"]=NA
res2[index_row,"t_evolplast_pval"]=if(mecha %in% c("total","without_adapt")){summary(lm(plast~index,data=bibi))$coeff["index","Pr(>|t|)"]}else{1}

#estimation de l'adaptation
res2[index_row,"adapt"]=suma$coeff["index","Estimate"]
#estimation de la plasticité
res2[index_row,"plast"]=suma$coeff["temp","Estimate"]
#estimation de l'interaction
res2[index_row,"evol_plast"]=if("temp:index" %in% rownames(suma$coeff)){
suma$coeff["temp:index","Estimate"]}else{NA}
#les erreurs faites sur les estimation d'adaptation et de plasticité
res2[index_row,"error_adapt"]=suma$coeff["index","Std. Error"]
res2[index_row,"error_plast"]=suma$coeff["temp","Std. Error"]
res2[index_row,"error_inter"]=if("temp:index" %in% rownames(suma$coeff)){
suma$coeff["temp:index","Std. Error"]}else{NA}
#les pvalues d'adaptation et de plasticité
res2[index_row,c("pval_plast","pval_adapt","pval_inter")]=t(ano[1:3,4])

#estimation de l'adaptation lmer
res2[index_row,"adapt_wi"]=suma2$coeff["index","Estimate"]
#estimation de la plasticité lmer
res2[index_row,"plast_wi"]=suma2$coeff["temp","Estimate"]
#estimation de l'interaction
res2[index_row,"evol_plast_wi"]=if("temp:index" %in% rownames(suma2$coeff)){
suma2$coeff["temp2","Estimate"]}else{NA}
#les erreurs faites sur les estimation d'adaptation et de plasticité lmer
res2[index_row,"error_adapt_wi"]=suma2$coeff["index","Std. Error"]
res2[index_row,"error_plast_wi"]=suma2$coeff["temp","Std. Error"]
res2[index_row,"error_inter_wi"]=if("temp:index2" %in% rownames(suma2$coeff)){
suma2$coeff["temp:index2","Std. Error"]}else{NA}
#les pvalues d'adaptation et de plasticité
res2[index_row,c("pval_plast_wi","pval_adapt_wi","pval_inter_wi")]=t(ano2[1:3,4])

#pval t_adapt
res2[index_row,"pval_t_adapt"]=summary(lm(mfd_adapt~index,data=bibi))$coeff[2,4]
#pval t_plast
res2[index_row,"pval_t_evol_plast"]=summary(lm(plast~index,data=bibi))$coeff[2,4]

#annee_eval
res2[index_row,"annee_eval"]=annee_eval
#temp_eval
res2[index_row,"temp_eval"]=temp_eval
#G
res2[index_row,"G_adapt"]=G[1,1]
res2[index_row,"G_plast"]=G[2,2]
res2[index_row,"mechanism"]=mecha
}
}
# print(ggplot()+geom_point(data=res,aes(x=index,y=records),alpha=0.4)+geom_point(data=res,aes(x=index,y=mfd_true),col="firebrick",size=1.4)+
# geom_point(data=res,aes(x=index,y=mfd_adapt),col="dodgerblue",size=1.4)+geom_pointrange(data=newres,aes(x=index,y=fit,ymin=lwr,ymax=upr),col="green")+
# geom_abline(intercept=model$coef[1],slope=model$coef[3]))
# print(paste("data_number = ",nbdata_vec[div]," - number of years = ",nbit," - standard deviation = ",vec_z[z],sep="")) 


setwd(dir="/home/duchenne/plast/resultats_simues/")
write.table(res2,paste("resultats_simulations_sans_evol_de_g_",nb_data,"_",mecha,"_",nbit,".txt",sep=""),row.names=F,sep="\t")

