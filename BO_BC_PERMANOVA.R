
library(devtools)
library(vegan)
library(scales)
library(RColorBrewer)
library(VennDiagram)
library(gplots)
library(car)
library(ggVennDiagram)
library(ggplot2)
library(eulerr)
library(ggpubr)
library(dplyr)
library(reshape)
library(reshape2)
library(data.table)
library(tidyverse)
library(lme4)
library(rsq)
library(emmeans)
library(stringr)
library(pheatmap)
library(gdata)
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
if (!require("rstudioapi")){
  install.packages("rstudioapi", dependencies = TRUE)
  library(rstudioapi)}


#### Get your table ####
setwd(dirname(getActiveDocumentContext()$path))
codes <- read.delim("metadata_temporal_bo.txt", header = TRUE, sep = "\t", dec = ".")
micro <- read.table("asv_taxo_temporal_bo.txt", sep="\t", header=T, dec = ".")
Class <- micro$ID
micro <- micro[,c(1:210)]
micro <- mutate_all(micro, function(x) as.numeric(as.character(x)))
microrf<-as.data.frame(prop.table(as.matrix(micro),2))

#### Get your subtables####
Aigua<-microrf[,which((codes$TISSUE== "WATER"))]
codesA<-codes[which((codes$TISSUE== "WATER")),]
Aigua <- Aigua[apply(Aigua[,-1], 1, function(x) !all(x==0)),]
Aigua<-prop.table(as.matrix(Aigua),2)
Aigua<-as.data.frame(Aigua)

Tunica<-micro[,which((codes$TISSUE== "TUNIC"))]
codesT<-codes[which((codes$TISSUE== "TUNIC")),]
Tunica <- Tunica[apply(Tunica[,-1], 1, function(x) !all(x==0)),]
Tunica<-prop.table(as.matrix(Tunica),2)
Tunica<-as.data.frame(Tunica)

Branquia<-microrf[,which((codes$TISSUE== "GILL"))]
codesB<-codes[which((codes$TISSUE== "GILL")),]
Branquia <- Branquia[apply(Branquia[,-1], 1, function(x) !all(x==0)),]
Branquia<-prop.table(as.matrix(Branquia),2)
Branquia<-as.data.frame(Branquia)

TunicaB<-micro[,which((codes$TISSUE== "TUNIC" & codes$POP== "BLA" ))]
codesTB<-codes[which((codes$TISSUE== "TUNIC" & codes$POP== "BLA" )),]
TunicaB <- TunicaB[apply(TunicaB[,-1], 1, function(x) !all(x==0)),]
TunicaB<-prop.table(as.matrix(TunicaB),2)
TunicaB<-as.data.frame(TunicaB)

BranquiaB<-micro[,which((codes$TISSUE== "GILL" & codes$POP== "BLA" ))]
codesBB<-codes[which((codes$TISSUE== "GILL" & codes$POP== "BLA" )),]
BranquiaB <- BranquiaB[apply(BranquiaB[,-1], 1, function(x) !all(x==0)),]
BranquiaB<-prop.table(as.matrix(BranquiaB),2)
BranquiaB<-as.data.frame(BranquiaB)

TunicaV<-micro[,which((codes$TISSUE== "TUNIC" & codes$POP== "VIL" ))]
codesTV<-codes[which((codes$TISSUE== "TUNIC" & codes$POP== "VIL" )),]
TunicaV <- TunicaV[apply(TunicaV[,-1], 1, function(x) !all(x==0)),]
TunicaV<-prop.table(as.matrix(TunicaV),2)
TunicaV<-as.data.frame(TunicaV)

BranquiaV<-micro[,which((codes$TISSUE== "GILL" & codes$POP== "VIL" ))]
codesBV<-codes[which((codes$TISSUE== "GILL" & codes$POP== "VIL" )),]
BranquiaV <- BranquiaV[apply(BranquiaV[,-1], 1, function(x) !all(x==0)),]
BranquiaV<-prop.table(as.matrix(BranquiaV),2)
BranquiaV<-as.data.frame(BranquiaV)

#### Calculate your BC distances ####

BC_all<-vegdist(t(microrf))
BCA<-vegdist(t(Aigua))
BCT<-vegdist(t(Tunica))
BCB<-vegdist(t(Branquia))
BCTB<-vegdist(t(TunicaB))
BCBB<-vegdist(t(BranquiaB))
BCTV<-vegdist(t(TunicaV))
BCBV<-vegdist(t(BranquiaV))

#### Perform your PERMANOVAs ####
perm_all<- adonis2(BC_all~codes$TISSUE*codes$POP*codes$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BC_all, codes$POP, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BC_all,codes$TISSUE, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BC_all,codes$POP, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BC_all,codes$INDVAL, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

perm_aigua<- adonis2(BCA~codesA$POP*codesA$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BCA, codesA$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BCA,codesA$POP, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCA,codesA$MONTH, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCA,codesA$YEAR, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCA,codesA$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

perm_tun<-adonis2(BCT~codesT$POP*codesT$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BCT, codesT$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BCT,codesT$POP, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCT,codesT$MONTH, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCT,codesT$YEAR, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCT,codesT$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

perm_gill<-adonis2(BCB~codesB$POP*codesB$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BCB, codesB$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BCB,codesB$POP, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCB,codesB$MONTH, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCB,codesB$YEAR, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCB,codesB$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

perm_tun_bla<-adonis2(BCTB~codesTB$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BCTB, codesTB$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BCTB,codesTB$MONTH, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCTB,codesTB$YEAR, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCTB,codesTB$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

perm_bra_bla<- adonis2(BCBB~codesBB$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BCBB, codesBB$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BCBB,codesBB$MONTH, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCBB,codesBB$YEAR, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCBB,codesBB$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

perm_bra_vil<-adonis2(BCBV~codesBV$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BCBV, codesBV$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BCBV,codesBV$MONTH, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCBV,codesBV$YEAR, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCBV,codesBV$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

perm_tun_vil<-adonis2(BCTV~codesTV$INDVAL, sqrt.dist=F, by="terms")
permdisp<- betadisper(BCTV, codesTV$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdisp)
boxplot(permdisp)

permutest(betadisper(BCTV,codesTV$MONTH, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCTV,codesTV$YEAR, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BCTV,codesTV$SAMPLING, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

####Do pairwise####

codes$inter4 <- paste(codes$TISSUE, codes$POP, codes$INDVAL, sep ="-")
codesA$inter3 <-  paste(codesA$POP, codesA$INDVAL,  sep ="-") 
codesT$inter3 <-  paste(codesT$POP, sep ="-") 
codesB$inter3 <-  paste(codesB$POP, codesB$INDVAL, sep ="-") 

pair_BC_all <- pairwise.adonis2(BC_all~inter4, codes, p.adjust.m = "BY")
pair_BCA <- pairwise.adonis2(BCA~inter3, codesA, p.adjust.m = "BY")
pair_BCT <- pairwise.adonis2(BCT~inter3, codesT, p.adjust.m = "BY")
pair_BCB <- pairwise.adonis2(BCB~inter3, codesB, p.adjust.m = "BY")

#### Export relevant tables ####

write.table(perm_all, "permanova_all.txt", sep ="\t", dec=".")
write.table(perm_aigua, "permanova_aigua.txt", sep ="\t", dec=".")
write.table(perm_gill, "permanova_branquia.txt", sep ="\t", dec=".")
write.table(perm_tun, "permanova_tunica.txt", sep ="\t", dec=".")
write.table(perm_bra_bla, "permanova_bla_bra.txt", sep ="\t", dec=".")
write.table(perm_tun_bla, "permanova_bla_tun.txt", sep ="\t", dec=".")
write.table(perm_bra_vil, "permanova_vil_bra.txt", sep ="\t", dec=".")
write.table(perm_tun_vil, "permanova_vil_tun.txt", sep ="\t", dec=".")

write.table(pair_BC_all, "pairwise_all.txt", sep ="\t", dec=".")
write.table(pair_BCA, "pairwise_aigua.txt", sep ="\t", dec=".")
write.table(pair_BCB, "pairwise_branquia.txt", sep ="\t", dec=".")
write.table(pair_BCT, "pairwise_tunica.txt", sep ="\t", dec=".")