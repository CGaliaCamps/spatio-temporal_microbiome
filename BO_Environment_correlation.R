
##Primer aneu a Session, set working directory i seleccioneu el directori on es trobaran tots els fitxers d'entrada i sortida

##A continuaci? carreguem els paquets que farem servir

if (!require("vcfR")) {
  install.packages("vcfR", dependencies = TRUE)
  library(vcfR)
}
if (!require("adegenet")) {
  install.packages("adegenet", dependencies = TRUE)
  library(adegenet)
}
if (!require("hierfstat")) {
  install.packages("hierfstat", dependencies = TRUE)
  library(hierfstat)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("vegan")) {
  install.packages("vegan", dependencies = TRUE)
  library(vegan)
} 
if (!require("poppr")) {
  install.packages("poppr", dependencies = TRUE)
  library(poppr)
}
if (!require("pegas")) {
  install.packages("pegas", dependencies = TRUE)
  library(pegas)
}
if (!require("psych")) {
  install.packages("psych", dependencies = TRUE)
  library(psych)
}


setwd("G:/styela/Molecular/Microbioma/mirobioma_temporal")

###Read the pollutants file
env <- read.table("metadata_pollutants_aigua.txt", sep = "\t", header = TRUE)
row.names(env)<-env$ID
colnames(env)[5] <- "ºC"


#check for correlated values
#TE
pdf(file="correlacions_pollutants.pdf", width=7, height=7)
pairs.panels(env[7:15], scale=T, density=T, stars=F, ellipse=F, hist.col="grey80", 
             lm=F, smooth=T, digits=1, rug=T, cex=1.1, cex.cor=2, smoother=T, method="pearson") 
dev.off()

#Salinity and Temperature
pdf(file="correlacions_environment.pdf", width=3, height=3)
pairs.panels(env[5:6], scale=T, density=T, stars=F, ellipse=F, hist.col="grey80", 
             lm=F, smooth=T, digits=1, rug=T, cex=1.1, cex.cor=2, smoother=T, method="pearson") 
dev.off()


