library(data.table)
setwd(dir="/home/duchenne/plast/species_temp_sel")
liste=list.files()
liste=liste[grep("selec_temp",liste,fixed=T)]
bidonf=NULL

for(i in liste){
bidon=fread(i)
bidonf=rbind(bidonf,bidon)
}

setwd(dir="/home/duchenne/plast/")
fwrite(bidonf,"selec_temp_var_beta.txt",sep="\t")