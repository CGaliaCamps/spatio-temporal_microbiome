library("devtools")
# install_github("microbiome/microbiome")
# install_github("microsud/jeevanuDB")
library(microbiome)
library(jeevanuDB)
library(RColorBrewer)
library(viridis)

sample <-read.table("asv_taxo_temporal_bo.txt",stringsAsFactors = F)
rownames(sample) <- sample$ASV

# Delete non-numeric samples
sample <- sample[,1:210]
# Relative Frequencies
sample<-as.matrix(sample)
sample <- as.data.frame(sample)

# Metadata time
codes <- read.table("metadata_temporal_bo.txt",header=T,stringsAsFactors = F)

#grab your samples
sample_T<-sample[,which((codes$TISSUE=="TUNIC"))]
sample_G<-sample[,which((codes$TISSUE=="GILL"))]
sample_W<-sample[,which((codes$TISSUE=="WATER"))]

sample_TW<-sample[,which((codes$TISSUE=="TUNIC"&codes$INDVAL=="WARM"))]
sample_TC<-sample[,which((codes$TISSUE=="TUNIC"&codes$INDVAL=="COLD"))]
sample_GW<-sample[,which((codes$TISSUE=="GILL"&codes$INDVAL=="WARM"))]
sample_GC<-sample[,which((codes$TISSUE=="GILL"&codes$INDVAL=="COLD"))]
sample_WW<-sample[,which((codes$TISSUE=="WATER"&codes$INDVAL=="WARM"))]
sample_WC<-sample[,which((codes$TISSUE=="WATER"&codes$INDVAL=="COLD"))]

# Remove ASVs with 0 reads

sample_T<-sample_T[which(rowSums(sample_T)>0),]
sample_G<-sample_G[which(rowSums(sample_G)>0),]
sample_W<-sample_W[which(rowSums(sample_W)>0),]

sample_TW<-sample_TW[which(rowSums(sample_TW)>0),]
sample_GW<-sample_GW[which(rowSums(sample_GW)>0),]
sample_WW<-sample_WW[which(rowSums(sample_WW)>0),]

sample_TC<-sample_TC[which(rowSums(sample_TC)>0),]
sample_GC<-sample_GC[which(rowSums(sample_GC)>0),]
sample_WC<-sample_WC[which(rowSums(sample_WC)>0),]

#Calculate relative abundances
sample_tot <- microbiome::transform(sample, "compositional")
sample_T <- microbiome::transform(sample_T, "compositional")
sample_G <- microbiome::transform(sample_G, "compositional")
sample_W <- microbiome::transform(sample_W, "compositional")

sample_TW <- microbiome::transform(sample_TW, "compositional")
sample_GW <- microbiome::transform(sample_GW, "compositional")
sample_WW <- microbiome::transform(sample_WW, "compositional")
sample_TC <- microbiome::transform(sample_TC, "compositional")
sample_GC <- microbiome::transform(sample_GC, "compositional")
sample_WC <- microbiome::transform(sample_WC, "compositional")

#Calculate 95% core
core_tot_9 <- core_members(sample_tot, detection = 1/10000000000000, prevalence = 0.95)
core_T_9 <- core_members(sample_T, detection = 1/10000000000000, prevalence = 0.95)
core_G_9 <- core_members(sample_G, detection = 1/10000000000000, prevalence = 0.95)
core_W_9 <- core_members(sample_W, detection = 1/10000000000000, prevalence = 0.95)

core_TW_9 <- core_members(sample_TW, detection = 1/10000000000000, prevalence = 0.95)
core_GW_9 <- core_members(sample_GW, detection = 1/10000000000000, prevalence = 0.95)
core_WW_9 <- core_members(sample_WW, detection = 1/10000000000000, prevalence = 0.95)
core_TC_9 <- core_members(sample_TC, detection = 1/10000000000000, prevalence = 0.95)
core_GC_9 <- core_members(sample_GC, detection = 1/10000000000000, prevalence = 0.95)
core_WC_9 <- core_members(sample_WC, detection = 1/10000000000000, prevalence = 0.95)

#Calculate 100% core
core_tot_1 <- core_members(sample_tot, detection = 1/10000000000000, prevalence = 1)
core_T_1 <- core_members(sample_T, detection = 1/10000000000000, prevalence = 1)
core_G_1 <- core_members(sample_G, detection = 1/10000000000000, prevalence = 1)
core_W_1 <- core_members(sample_W, detection = 1/10000000000000, prevalence = 1)

core_TW_1 <- core_members(sample_TW, detection = 1/10000000000000, prevalence = 1)
core_GW_1 <- core_members(sample_GW, detection = 1/10000000000000, prevalence = 1)
core_WW_1 <- core_members(sample_WW, detection = 1/10000000000000, prevalence = 1)
core_TC_1 <- core_members(sample_TC, detection = 1/10000000000000, prevalence = 1)
core_GC_1 <- core_members(sample_GC, detection = 1/10000000000000, prevalence = 1)
core_WC_1 <- core_members(sample_WC, detection = 1/10000000000000, prevalence = 1)


detections <- round(10^seq(log10(1e-5), log10(0.5), length = 20), 4)

#plot data
pdf("sample_tot_9.pdf",width=9, height=2.5)
plot_core(sample_tot,
          plot.type = "heatmap",
          colours = c("white","purple4"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.00001, length = 1000), 5), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold (%)") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_T_9.pdf",width=3, height=2)
plot_core(sample_T,
                plot.type = "heatmap",
                colours = c("white","#E67E22"),
                prevalences = c(1/10000000000000,1),
                detections = round(seq(1/10000000000000, 0.5, length = 100), 3), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_G_9.pdf",width=3, height=2)
plot_core(sample_G,
          plot.type = "heatmap",
          colours = c("white","#2E2E2E"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.2, length = 100), 3), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_W_9.pdf",width=3, height=5.2)
plot_core(sample_W,
          plot.type = "heatmap",
          colours = c("white","#0B00FB"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.5, length = 100), 3), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_TW_9.pdf",width=3, height=2)
plot_core(sample_TW,
          plot.type = "heatmap",
          colours = c("white","#E67E22"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.2, length = 100), 3), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_GW_9.pdf",width=3, height=2)
plot_core(sample_GW,
          plot.type = "heatmap",
          colours = c("white","#2E2E2E"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.2, length = 100), 3), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_WW_9.pdf",width=3, height=5.2)
plot_core(sample_WW,
          plot.type = "heatmap",
          colours = c("white","#0B00FB"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.2, length = 100), 3),
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text(size = 6))
dev.off()


pdf("sample_TC_9.pdf",width=3, height=2)
plot_core(sample_TC,
          plot.type = "heatmap",
          colours = c("white","#E67E22"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.5, length = 100), 3), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_GC_9.pdf",width=3, height=2)
plot_core(sample_GC,
          plot.type = "heatmap",
          colours = c("white","#2E2E2E"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.2, length = 100), 3), 
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text())
dev.off()

pdf("sample_WC_9.pdf",width=3, height=5.2)
plot_core(sample_WC,
          plot.type = "heatmap",
          colours = c("white","#0B00FB"),
          prevalences = c(1/10000000000000,1),
          detections = round(seq(1/10000000000000, 0.5, length = 100), 3),
          min.prevalence = 0.95) +
  xlab("Abundance threshold") +
  theme_bw() + ylab("ASVs") + theme(axis.text.x = element_text(angle=90, vjust=0.5),
                                    axis.text.y = element_text(size = 6))
dev.off()


#generate an upset plot for core ASV
#install.packages("splitstackshape")
library(splitstackshape)
library(reshape2)
library(UpSetR)
list_core <- list(All= core_tot_9, 
                  Tunic= core_T_9, Gill=core_G_9, Water=core_W_9,
                  Warm_tunic= core_TW_9, Warm_gill=core_GW_9, Warm_water=core_WW_9,
                  Cold_tunic= core_TC_9, Cold_gill=core_GC_9, Cold_water=core_WC_9)
list_plot <- t(splitstackshape:::charMat(list_core, fill = 0, mode="binary"))
colnames(list_plot) <- c("All", "Tunic", "Gill", "Seawater","Tunic_warm","Gill_warm","Seawater_warm","Tunic_cold","Gill_cold","Seawater_cold") 
list_plot <- list_plot[,c(1,4,7,10,2,5,8,3,6,9)]
list_plot <- rev(as.data.frame(list_plot))
sets <-    colnames(list_plot)              


pdf("core_summary_bo_name.pdf", height=2.5, width=5.7)
UpSetR::upset(list_plot, sets = sets, mb.ratio = c(0.55, 0.45), 
              sets.bar.color = c("#2E2E2E","#2E2E2E","black", "#E67E22","#E67E22","#E67E22", "#0B00FB","#0B00FB","#0B00FB", "#000000"), 
              order.by = c("freq"), keep.order=T, line.size = 0.5
              ,sets.x.label = "Number of core ASV",
              mainbar.y.label = "Shared ASV",
              show.numbers = "yes"
              )
dev.off()


write.table(list_plot, "presence-absence_core.txt")

