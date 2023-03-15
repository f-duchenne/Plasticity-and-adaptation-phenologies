library(data.table)
library(dplyr)
library(R2jags)
# setup parallel backend to use all available CPUs
#ncpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
#registerDoParallel(ncpus)

# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])
print(i)

setwd(dir="/home/duchenne/plast/")
dat=fread("selec_temp_var_beta.txt",sep="\t")
liste=data.frame(Speciesgen=unique(dat$Speciesgen),nb_pre1990=dat$nb_pre1990[!duplicated(dat$Speciesgen)],nb_post1990=dat$nb_post1990[!duplicated(dat$Speciesgen)])
liste$nb=liste$nb_pre1990+liste$nb_post1990

resf=NULL

bidon=fread(paste0("data_per_species/",liste$Speciesgen[i],".txt"))
bidon=subset(bidon, !is.na(elev) & !is.na(temp_0_90))

bidon$Altitude2=bidon$Altitude
bidon$Altitude=(bidon$Altitude-mean(bidon$Altitude,na.rm=T))/1000
bidon$Annee2=bidon$Annee-1990
bidon$Annee=bidon$Annee-min(bidon$Annee)+1
bidon$maille_slope=bidon$maille
varia=subset(dat,Speciesgen==liste$Speciesgen[i])$varia
bidon$maillen=as.numeric(as.factor(bidon$maille))

sites=unique(bidon[,c("maille","ctr_lon","ctr_lat")])
sites=sites[order(sites$maille),]
Distance=dist(sites[,c("ctr_lon","ctr_lat")], method="euclidean", diag=TRUE, upper=TRUE)

bidon=as.data.frame(bidon)
dat=c(as.list(bidon),list(Ndata=nrow(bidon),Ntemp=max(bidon$Annee),mailleu=unique(bidon$maille),b=c(1,1),Distance=as.matrix(Distance),Nsites=nrow(sites),
varia1=bidon[,varia[1]],varia2=bidon[,varia[2]]))


######################################################################################################
######################################################################################################
#### JAGS model file written by runjags version 2.2.0-2 on 2021-11-23 15:22:45 
######################################################################################################
######################################################################################################

### Model template as follows - ensure this is syntactically correct before running the model!
if(length(varia)==2){
model_string="
model{

# Process model
for(i in 1:Ndata){
	mu[i] <- beta0 +  (beta1+beta1_delta[maillen[i]])* varia1[i] +  (beta2+beta2_delta[maillen[i]])* varia2[i] + (beta3+beta3_delta[maillen[i]])*Annee2[i]+beta4*Altitude[i]+
	theta_A[Annee[i]]+theta_m[maillen[i]]
	Jour.de.collecte[i] ~ dnorm(mu[i],tau.resid)
	}

#Jour.de.collecte[1:Ndata] ~ dmnorm.vcov(mu[1:Ndata], spatialvariance * D.covar[1:Nsites,1:Nsites])


beta0 ~ dunif(0,365)
beta1 ~ dnorm(0,0.01)
beta2 ~ dnorm(0,0.01)
beta3 ~ dnorm(0,1)
beta4 ~ dnorm(0,0.01)

#spatial autocorrelation:
edec0 ~ dgamma(3, 0.1)
edec1 ~ dgamma(3, 0.1)
edec2 ~ dgamma(3, 0.1)
edec3 ~ dgamma(3, 0.1)
for(j in 1:Nsites){
for(i in 1:Nsites){
D.covar0[i,j] <- exp(-edec0*Distance[i,j])
D.covar1[i,j] <- exp(-edec1*Distance[i,j])
D.covar2[i,j] <- exp(-edec2*Distance[i,j])
D.covar3[i,j] <- exp(-edec3*Distance[i,j])
}}


tau.resid=1/(sd.resid * sd.resid)
sd.resid ~ dt(0, 1, 1)T(0,)

theta_m[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.maille * sd.maille) * D.covar0[1:Nsites,1:Nsites])
beta1_delta[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.beta1 * sd.beta1) * D.covar1[1:Nsites,1:Nsites])
beta2_delta[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.beta2 * sd.beta2) * D.covar2[1:Nsites,1:Nsites])
beta3_delta[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.beta3 * sd.beta3) * D.covar3[1:Nsites,1:Nsites])

# for(j in mailleu){
# theta_m[j] ~ dnorm(0, tau.maille)
# beta1_delta[j]~ dnorm(0, tau.beta1)
# beta2_delta[j]~ dnorm(0, tau.beta2)
# beta3_delta[j]~ dnorm(0, tau.beta3)
# }
tau.maille<- 1/(sd.maille * sd.maille)
sd.maille ~ dt(0, 1, 1)T(0,)
tau.beta1<- 1/(sd.beta1 * sd.beta1)
sd.beta1 ~ dt(0, 1, 1)T(0,)
tau.beta2<- 1/(sd.beta2 * sd.beta2)
sd.beta2 ~ dt(0, 1, 1)T(0,)
tau.beta3<- 1/(sd.beta3 * sd.beta3)
sd.beta3 ~ dt(0, 1, 1)T(0,)

theta_A[1] ~ dnorm(0, tau.temp)
for(j in 2:Ntemp){
theta_A[j] ~ dnorm(theta_A[j-1], tau.temp)
}
tau.temp<- 1/(sd.temp * sd.temp)
sd.temp ~ dt(0, 1, 1)T(0,)

beta1_deltat=beta1+beta1_delta
beta2_deltat=beta2+beta2_delta
beta3_deltat=beta3+beta3_delta
theta_mt=theta_m+beta0

}
"
}

if(length(varia)==1){
model_string="
model{

# Process model
for(i in 1:Ndata){
	mu[i] <- beta0 +  (beta1+beta1_delta[maillen[i]])* varia1[i] + (beta3+beta3_delta[maillen[i]])*Annee2[i]+beta4*Altitude[i]+
	theta_A[Annee[i]]+theta_m[maillen[i]]
	Jour.de.collecte[i] ~ dnorm(mu[i],tau.resid)
	}

#Jour.de.collecte[1:Ndata] ~ dmnorm.vcov(mu[1:Ndata], spatialvariance * D.covar[1:Nsites,1:Nsites])


beta0 ~ dunif(0,365)
beta1 ~ dnorm(0,0.01)
beta3 ~ dnorm(0,1)
beta4 ~ dnorm(0,0.01)

#spatial autocorrelation:
edec0 ~ dgamma(3, 0.1)
edec1 ~ dgamma(3, 0.1)
edec3 ~ dgamma(3, 0.1)
for(j in 1:Nsites){
for(i in 1:Nsites){
D.covar0[i,j] <- exp(-edec0*Distance[i,j])
D.covar1[i,j] <- exp(-edec1*Distance[i,j])
D.covar3[i,j] <- exp(-edec3*Distance[i,j])
}}


tau.resid=1/(sd.resid * sd.resid)
sd.resid ~ dt(0, 1, 1)T(0,)

theta_m[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.maille * sd.maille) * D.covar0[1:Nsites,1:Nsites])
beta1_delta[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.beta1 * sd.beta1) * D.covar1[1:Nsites,1:Nsites])
beta3_delta[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.beta3 * sd.beta3) * D.covar3[1:Nsites,1:Nsites])

# for(j in mailleu){
# theta_m[j] ~ dnorm(0, tau.maille)
# beta1_delta[j]~ dnorm(0, tau.beta1)
# beta2_delta[j]~ dnorm(0, tau.beta2)
# beta3_delta[j]~ dnorm(0, tau.beta3)
# }
tau.maille<- 1/(sd.maille * sd.maille)
sd.maille ~ dt(0, 1, 1)T(0,)
tau.beta1<- 1/(sd.beta1 * sd.beta1)
sd.beta1 ~ dt(0, 1, 1)T(0,)
tau.beta3<- 1/(sd.beta3 * sd.beta3)
sd.beta3 ~ dt(0, 1, 1)T(0,)

theta_A[1] ~ dnorm(0, tau.temp)
for(j in 2:Ntemp){
theta_A[j] ~ dnorm(theta_A[j-1], tau.temp)
}
tau.temp<- 1/(sd.temp * sd.temp)
sd.temp ~ dt(0, 1, 1)T(0,)

beta1_deltat=beta1+beta1_delta
beta3_deltat=beta3+beta3_delta
theta_mt=theta_m+beta0

}
"
}

if(length(varia)==0){
model_string="
model{

# Process model
for(i in 1:Ndata){
	mu[i] <- beta0 + (beta3+beta3_delta[maillen[i]])*Annee2[i]+beta4*Altitude[i]+
	theta_A[Annee[i]]+theta_m[maillen[i]]
	Jour.de.collecte[i] ~ dnorm(mu[i],tau.resid)
	}

#Jour.de.collecte[1:Ndata] ~ dmnorm.vcov(mu[1:Ndata], spatialvariance * D.covar[1:Nsites,1:Nsites])


beta0 ~ dunif(0,365)
beta3 ~ dnorm(0,1)
beta4 ~ dnorm(0,0.01)

#spatial autocorrelation:
edec0 ~ dgamma(3, 0.1)
edec3 ~ dgamma(3, 0.1)
for(j in 1:Nsites){
for(i in 1:Nsites){
D.covar0[i,j] <- exp(-edec0*Distance[i,j])
D.covar3[i,j] <- exp(-edec3*Distance[i,j])
}}


tau.resid=1/(sd.resid * sd.resid)
sd.resid ~ dt(0, 1, 1)T(0,)

theta_m[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.maille * sd.maille) * D.covar0[1:Nsites,1:Nsites])
beta3_delta[1:Nsites] ~ dmnorm.vcov(rep(0,Nsites), (sd.beta3 * sd.beta3) * D.covar3[1:Nsites,1:Nsites])

# for(j in mailleu){
# theta_m[j] ~ dnorm(0, tau.maille)
# beta1_delta[j]~ dnorm(0, tau.beta1)
# beta2_delta[j]~ dnorm(0, tau.beta2)
# beta3_delta[j]~ dnorm(0, tau.beta3)
# }
tau.maille<- 1/(sd.maille * sd.maille)
sd.maille ~ dt(0, 1, 1)T(0,)
tau.beta3<- 1/(sd.beta3 * sd.beta3)
sd.beta3 ~ dt(0, 1, 1)T(0,)

theta_A[1] ~ dnorm(0, tau.temp)
for(j in 2:Ntemp){
theta_A[j] ~ dnorm(theta_A[j-1], tau.temp)
}
tau.temp<- 1/(sd.temp * sd.temp)
sd.temp ~ dt(0, 1, 1)T(0,)

beta3_deltat=beta3+beta3_delta
theta_mt=theta_m+beta0

}
"
}


setwd(dir="/home/duchenne/plast/resultats")
writeLines(model_string,con=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),".txt"))

ParsStage <- c("beta0","beta1","beta2","beta3","beta4","sd.maille","sd.temp","sd.beta1","sd.beta2","sd.beta3","theta_A","theta_mt","edec0","edec1","edec2","edec3",
"beta1_deltat","beta2_deltat","beta3_deltat")

Inits <- function(){list()}

results1 <- jags(model.file=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),".txt"), parameters.to.save=ParsStage, n.chains=1, 
data=dat,n.burnin = 15000,n.iter = 30000, n.thin = 3,inits =Inits)


save(results1,file=paste0("model_",gsub(" ","_",liste$Speciesgen[i],fixed=T),"_",liste$nchain[i],".RData"))