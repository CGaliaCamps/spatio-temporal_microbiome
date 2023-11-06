## Comencem amb les llibreries i el directori
if (!require("rstudioapi")){
  install.packages("rstudioapi", dependencies = TRUE)
  library(rstudioapi)}
if (!require("indicspecies")){
  install.packages("indicspecies", dependencies = TRUE)
  library(indicspecies)}
if (!require("readxl")){
  install.packages("readxl", dependencies = TRUE)
  library(readxl)} 

##Read abundance table
setwd("E:/styela/Molecular/Microbioma/mirobioma_temporal/Coses Liam")
sample <-read.table("asv_taxo_temporal_bo.txt",stringsAsFactors = F, header=T, sep="\t")
rownames(sample) <- sample$ASV
# Remove non-numeric columns
sample <- sample[,1:210]

# Relative Frequencies
sample<-as.matrix(sample)
sample<-prop.table(sample,2)
sample <- as.data.frame(sample)

# Metadata time
codes <- read.table("metadata_temporal_bo.txt",header=T,stringsAsFactors = F)
rownames(codes)<-codes$ID
codes <- codes[order(codes$ID),]

#Select each tissue
sample_T<-sample[,which((codes$TISSUE=="TUNIC"))]
sample_G<-sample[,which((codes$TISSUE=="GILL"))]
sample_W<-sample[,which((codes$TISSUE=="WATER"))]

# Remove all missing data
sample_T<-sample_T[which(rowSums(sample_T)>0),]
sample_G<-sample_G[which(rowSums(sample_G)>0),]
sample_W<-sample_W[which(rowSums(sample_W)>0),]

#transpose
sample<-t(sample)
sample_T<-t(sample_T)
sample_G<-t(sample_G)
sample_W<-t(sample_W)

codes_T<-codes[which(codes$TISSUE=="TUNIC"),]
codes_G<-codes[which(codes$TISSUE=="GILL"),]
codes_W<-codes[which(codes$TISSUE=="WATER"),]

group<-factor(codes$TISSUE,levels=unique(codes$TISSUE))
groupT<-factor(codes_T$INDVAL,levels=unique(codes_T$INDVAL))
groupG<-factor(codes_G$INDVAL,levels=unique(codes_G$INDVAL))
groupW<-factor(codes_W$INDVAL,levels=unique(codes_W$INDVAL))

 indval_tissue <- multipatt(sample,group,control=how(nperm=999),max.order=1, duleg=TRUE)
 indvalT<-multipatt(sample_T,groupT,control=how(nperm=999),max.order=1, duleg=TRUE)
 indvalG<-multipatt(sample_G,groupG,control=how(nperm=999),max.order=1, duleg=TRUE)
 indvalW<-multipatt(sample_W,groupW,control=how(nperm=999),max.order=1, duleg=TRUE)


 #Export tables
  sink("indval_tissue.txt")
  cat(summary(indval_tissue,indvalcomp = T))
  sink()

  sink("indvalT.txt")
  cat(summary(indvalT,indvalcomp = T))
  sink()

  sink("indvalG.txt")
  cat(summary(indvalG,indvalcomp = T))
  sink()

 sink("indvalW.txt")
 cat(summary(indvalW,indvalcomp = T))
 sink()


#add taxononomy to the matrix
taxo <- read.csv("taxonomia_final.csv",header=T,sep=",", stringsAsFactors = F)
#taxo <- gsub("indet","Unknown", taxo)
rownames(taxo) <- taxo$ASV.ID
taxo<- taxo[,c(3:10)]
sample$ID <- taxo$ASV.ID

#Load Indval data and keep only thos ASV with A & B above 0.7. (70%  reads in 70% samples)

IndvalG <- read.table("IndvalGill.txt", header=T, sep="\t")
IndvalT <- read.table("IndvalTunic.txt", header=T, sep="\t")
IndvalW <- read.table("IndvalWater.txt", header=T, sep="\t")

IndvalGW <- read.table("IndvalGW.txt", header=T, sep="\t")
IndvalGC <- read.table("IndvalGC.txt", header=T, sep="\t")
IndvalTW <- read.table("IndvalTW.txt", header=T, sep="\t")
IndvalTC <- read.table("IndvalTC.txt", header=T, sep="\t")
IndvalWW <- read.table("IndvalWW.txt", header=T, sep="\t")
IndvalWC <- read.table("IndvalWC.txt", header=T, sep="\t")

indvalG8 <-IndvalG[which(IndvalG$A>=0.7 & IndvalG$B>=0.7),] 
indvalT8 <-IndvalT[which(IndvalT$A>=0.7 & IndvalT$B>=0.7),] 
indvalW8 <-IndvalW[which(IndvalW$A>=0.7 & IndvalW$B>=0.7),] 
indvalTW8 <-IndvalTW[which(IndvalTW$A>=0.7 & IndvalTW$B>=0.7),] 
indvalTC8 <-IndvalTC[which(IndvalTC$A>=0.7 & IndvalTC$B>=0.7),] 
indvalGW8 <-IndvalGW[which(IndvalGW$A>=0.7 & IndvalGW$B>=0.7),] 
indvalGC8 <-IndvalGC[which(IndvalGC$A>=0.7 & IndvalGC$B>=0.7),] 
indvalWW8 <-IndvalWW[which(IndvalWW$A>=0.7 & IndvalWW$B>=0.7),] 
indvalWC8 <-IndvalWC[which(IndvalWC$A>=0.7 & IndvalWC$B>=0.7),] 


write.table(indvalTW8, "IndvalTW7.txt")
write.table(indvalTC8, "IndvalTC7.txt")
write.table(indvalGW8, "IndvalGW7.txt")
write.table(indvalGC8, "IndvalGC7.txt")
write.table(indvalWW8, "IndvalWW7.txt")
write.table(indvalWC8, "IndvalWC7.txt")


#Plot your venn diagrams
library(ggplot2)
library(VennDiagram)
Gill <- venn.diagram(x = list(IndvalG$ID, IndvalGW$ID, IndvalGC$ID), 
                        filename = NULL, fill=c("#2E2E2E", "#b05446","#9dcfdd"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                        euler.d = TRUE, scaled = TRUE,cat.fontfamily = "sans",cat.cex = 0.6,
                        print.mode = c("raw"),fontfamily = "sans",
                        category.names = c("Gill" , "Gill_warm" , "Gill_cold"), 
                        #cat.just=list(c(1,-40) , c(0.4,-40)),
                        lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                        cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(Gill, file="Gill_indval_shared.pdf", device = "pdf", height=3, width=3)

Tunic <- venn.diagram(x = list(IndvalT$ID, IndvalTW$ID, IndvalTC$ID), 
                     filename = NULL, fill=c("#E67E22", "#b05446","#9dcfdd"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                     euler.d = TRUE, scaled = TRUE,cat.fontfamily = "sans",cat.cex = 0.6,
                     print.mode = c("raw"),fontfamily = "sans",
                     category.names = c("Tunic" , "Tunic_warm" , "Tunic_cold"), 
                     #cat.just=list(c(1,-40) , c(0.4,-40)),
                     lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                     cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(Tunic, file="Tunic_indval_shared.pdf", device = "pdf", height=3, width=3)

Water <- venn.diagram(x = list(IndvalW$ID, IndvalWW$ID, IndvalWC$ID), 
                        filename = NULL, fill=c("#0B00FB", "#b05446","#9dcfdd"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                        euler.d = TRUE, scaled = TRUE,cat.fontfamily = "sans",cat.cex = 0.6,
                        print.mode = c("raw"),fontfamily = "sans",
                        category.names = c("Water" , "Water_warm" , "Water_cold"), 
                        #cat.just=list(c(1,-40) , c(0.4,-40)),
                        lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                        cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(Water, file="Water_indval_shared.pdf", device = "pdf", height=3, width=3)


WARM <- venn.diagram(x = list(IndvalGW$ID, IndvalWW$ID, IndvalTW$ID), 
                      filename = NULL, fill=c("#2E2E2E", "#0B00FB","#E67E22"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                      euler.d = TRUE, scaled = TRUE,cat.fontfamily = "sans",cat.cex = 0.6,
                      print.mode = c("raw"),fontfamily = "sans",
                      category.names = c("Gill_warm" , "Water_warm" , "Tunic_warm"), 
                      #cat.just=list(c(1,-40) , c(0.4,-40)),
                      lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                      cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(WARM, file="Warm_indval_shared.pdf", device = "pdf", height=3, width=3)

COLD <- venn.diagram(x = list(IndvalGC$ID, IndvalWC$ID, IndvalTC$ID), 
                      filename = NULL, fill=c("#2E2E2E", "#0B00FB","#E67E22"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                      euler.d = TRUE, scaled = TRUE,cat.fontfamily = "sans",cat.cex = 0.6,
                      print.mode = c("raw"),fontfamily = "sans",
                      category.names = c("Gill_cold" , "Water_cold" , "Tunic_cold"), 
                      #cat.just=list(c(1,-40) , c(0.4,-40)),
                      lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                      cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(COLD, file="Cold_indval_shared.pdf", device = "pdf", height=3, width=3)


Gill8 <- venn.diagram(x = list(indvalG8$ID, indvalGW8$ID, indvalGC8$ID), 
                     filename = NULL, fill=c("#2E2E2E", "#b05446","#9dcfdd"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                     euler.d = FALSE, scaled = FALSE,cat.fontfamily = "sans",cat.cex = 0.6,
                     print.mode = c("raw"),fontfamily = "sans",
                     category.names = c("" , "" , ""), 
                     #cat.just=list(c(1,-40) , c(0.4,-40)),
                     lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                     cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(Gill8, file="Gill8_indval_shared.pdf", device = "pdf", height=2.25, width=2.25)

Tunic8 <- venn.diagram(x = list(indvalT8$ID, indvalTW8$ID, indvalTC8$ID), 
                      filename = NULL, fill=c("#E67E22", "#b05446","#9dcfdd"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                      euler.d = FALSE, scaled = FALSE,cat.fontfamily = "sans",cat.cex = 0.6,
                      print.mode = c("raw"),fontfamily = "sans",
                      category.names = c("" , "" , ""), 
                      #cat.just=list(c(1,-40) , c(0.4,-40)),
                      lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                      cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(Tunic8, file="Tunic8_indval_shared.pdf", device = "pdf", height=2.25, width=2.25)

Water8 <- venn.diagram(x = list(indvalW8$ID, indvalWW8$ID, indvalWC8$ID), 
                      filename = NULL, fill=c("#0B00FB", "#b05446","#9dcfdd"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                      euler.d = FALSE, scaled = FALSE,cat.fontfamily = "sans",cat.cex = 0.6,
                      print.mode = c("raw"),fontfamily = "sans",
                      category.names = c("" , "" , ""), 
                      #cat.just=list(c(1,-40) , c(0.4,-40)),
                      lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                      cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(Water8, file="Water8_indval_shared.pdf", device = "pdf", height=2.25, width=2.25)


WARM8 <- venn.diagram(x = list( indvalWW8$ID, indvalTW8$ID, indvalGW8$ID), 
                     filename = NULL, fill=c("#0B00FB","#E67E22","#2E2E2E"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                     euler.d = FALSE, scaled = FALSE,cat.fontfamily = "sans",cat.cex = 0.6,
                     print.mode = c("raw"),fontfamily = "sans",
                     category.names = c("" , "" , ""), 
                     #cat.just=list(c(1,-40) , c(0.4,-40)),
                     lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                     cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(WARM8, file="Warm8_indval_shared.pdf", device = "pdf", height=2.25, width=2.25)

COLD8 <- venn.diagram(x = list( indvalWC8$ID, indvalTC8$ID, indvalGC8$ID), 
                     filename = NULL, fill=c( "#0B00FB","#E67E22","#2E2E2E"), alpha=c(0.5,0.5,0.5), disable.logging = TRUE,
                     euler.d = FALSE, scaled = FALSE,cat.fontfamily = "sans",cat.cex = 0.6,
                     print.mode = c("raw"),fontfamily = "sans",  inverted = FALSE,
                     category.names = c("" , "" , ""), 
                     #cat.just=list(c(1,-40) , c(0.4,-40)),
                     lwd = 2,lty = 'blank', cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
                     cex = 0.9,width = 10,height = 10, units = "in", cat.dist = c(0.055, 0.055, 0.085))
ggsave(COLD8, file="Cold8_indval_shared.pdf", device = "pdf", height=2.25, width=2.25)

