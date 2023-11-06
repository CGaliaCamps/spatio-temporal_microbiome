library(splitstackshape)
library(reshape2)
library(UpSetR)
library(ggplot2)
library(dplyr)

setwd("E:/styela/Molecular/Microbioma/mirobioma_temporal/Coses Liam")
sample <-read.table("asv_taxo_temporal_bo.txt", stringsAsFactors = F, header=T, sep = "\t")
rownames(sample) <- sample$ASV
# Remove non-numeroc columns
sample <- sample[,1:210]

# Relative Frequencies
sample<-as.matrix(sample)
sample<-prop.table(sample,2)
sample <- as.data.frame(sample)

# Metadata time
codes <- read.table("metadata_temporal_bo.txt",header=T,stringsAsFactors = F)
rownames(codes)<-codes$ID
codes <- codes[order(codes$ID),]
codes$SAMPLING<-as.character(codes$SAMPLING)
codes_tunic <- codes[which(codes$TISSUE!="GILL"),]
codes_gill <- codes[which(codes$TISSUE!="TUNIC"),]


####Trace Element selection####
setwd("G:/styela/Molecular/Microbioma/mirobioma_temporal")

PollutG <- read.table("pollut_gill_candidate_1axis.txt", header=T, sep="\t")
EnvG <- read.table("enviro_gill_candidate_1axis.txt", header=T, sep="\t")
PollutT <- read.table("pollut_tunic_candidate_1axis.txt", header=T, sep="\t")

polG <- sample[match(PollutG$ASV, rownames(sample)), which(codes$TISSUE!="TUNIC")]
envG <- sample[match(EnvG$ASV, rownames(sample)), which(codes$TISSUE!="TUNIC")]
polT <- sample[match(PollutT$ASV, rownames(sample)), which(codes$TISSUE!="GILL")] 

polG <- t(polG)
envG <- t(envG)
polT <-t(polT)


lines_polG<-aggregate(polG[,]~codes_gill$SAMPLING*codes_gill$TISSUE,polG,mean)
lines_envG<-aggregate(envG[,]~codes_gill$SAMPLING*codes_gill$TISSUE,envG,mean)
lines_polT<-aggregate(polT[,]~codes_tunic$SAMPLING*codes_tunic$TISSUE,polT,mean)

colnames(lines_polG)[1:2] <- colnames(lines_envG)[1:2] <- colnames(lines_polT)[1:2] <- c("Sampling","Tissue")

lines_polG$Sampling<-as.character(lines_polG$Sampling)
lines_envG$Sampling<-as.character(lines_envG$Sampling)
lines_polT$Sampling<-as.character(lines_polT$Sampling)

library(reshape2)
lines_polG<-melt(lines_polG)
lines_envG<-melt(lines_envG)
lines_polT<-melt(lines_polT)

lines_polG$Comp_ASV <- paste(lines_polG$Tissue, lines_polG$variable,  sep ="_")
lines_envG$Comp_ASV <- paste(lines_envG$Tissue, lines_envG$variable,  sep ="_")
lines_polT$Comp_ASV <- paste(lines_polT$Tissue, lines_polT$variable,  sep ="_")

lines_polG <- lines_polG[,c(1,5,4)]
lines_envG <- lines_envG[,c(1,5,4)]
lines_polT <- lines_polT[,c(1,5,4)]

lines_polG <- dcast(data = lines_polG,formula = Sampling~Comp_ASV,fun.aggregate = sum,value.var = "value")
lines_envG <- dcast(data = lines_envG,formula = Sampling~Comp_ASV,fun.aggregate = sum,value.var = "value")
lines_polT <- dcast(data = lines_polT,formula = Sampling~Comp_ASV,fun.aggregate = sum,value.var = "value")

####Environmental data####
codes_plot_e <- melt(codes[,c(4:6)])
codes_plot_p <- melt(codes[,c(4,7:15)])


codes_plot_e_mean <- dcast(data = codes_plot_e,formula = SAMPLING~variable,fun.aggregate = mean,value.var = "value")
codes_plot_e_sd <- dcast(data = codes_plot_e,formula = SAMPLING~variable,fun.aggregate = sd,value.var = "value")

codes_plot_p_mean <- dcast(data = codes_plot_p,formula = SAMPLING~variable,fun.aggregate = mean,value.var = "value")
codes_plot_p_sd <- dcast(data = codes_plot_p,formula = SAMPLING~variable,fun.aggregate = sd,value.var = "value")

codes_plot_e_mean <- melt(codes_plot_e_mean)
codes_plot_e_sd <- melt(codes_plot_e_sd)

codes_plot_p_mean <- melt(codes_plot_p_mean)
codes_plot_p_sd <- melt(codes_plot_p_sd)

desvest_E <- codes_plot_e_sd$value
desvest_P <- codes_plot_p_sd$value



codes_plot_p_vil <- melt(codes[which(codes$POP=="VIL"),c(4,7:15)])
codes_plot_p_bla <- melt(codes[which(codes$POP=="BLA"),c(4,7:15)])

####Indicator ASV####
setwd("E:/styela/Molecular/Microbioma/mirobioma_temporal/Coses Liam")
indvalTW8 <- read.table("IndvalTW7.txt", header=T)
indvalTC8 <- read.table("IndvalTC7.txt")
indvalGW8 <- read.table("IndvalGW7.txt")
indvalGC8 <- read.table("IndvalGC7.txt")
indvalWW8 <- read.table("IndvalWW7.txt")
indvalWC8 <- read.table("IndvalWC7.txt")

sample$ID <- rownames(sample)

IndvalTW8_freq<-sample %>% semi_join(indvalTW8, by ="ID")
IndvalTC8_freq<-sample %>% semi_join(indvalTC8, by ="ID")
IndvalGW8_freq<-sample %>% semi_join(indvalGW8, by ="ID")
IndvalGC8_freq<-sample %>% semi_join(indvalGC8, by ="ID")
IndvalWW8_freq<-sample %>% semi_join(indvalWW8, by ="ID")
IndvalWC8_freq<-sample %>% semi_join(indvalWC8, by ="ID")

tuniW <- IndvalTW8_freq[, which(codes$TISSUE!="GILL")] 
tuniC <- IndvalTC8_freq[, which(codes$TISSUE!="GILL")]
gillW <- IndvalGW8_freq[, which(codes$TISSUE!="TUNIC")] 
gillC <- IndvalGC8_freq[, which(codes$TISSUE!="TUNIC")]


tuniW <- t(tuniW)
tuniC <- t(tuniC)
gillW <- t(gillW)
gillC <- t(gillC)


lines_tuniW<-aggregate(tuniW[,]~codes_tunic$SAMPLING*codes_tunic$TISSUE,tuniW,mean)
lines_tuniC<-aggregate(tuniC[,]~codes_tunic$SAMPLING*codes_tunic$TISSUE,tuniC,mean)
lines_gillW<-aggregate(gillW[,]~codes_gill$SAMPLING*codes_gill$TISSUE,gillW,mean)
lines_gillC<-aggregate(gillC[,]~codes_gill$SAMPLING*codes_gill$TISSUE,gillC,mean)

colnames(lines_tuniW)[1:2] <- c("Sampling","Tissue")
colnames(lines_tuniC)[1:2] <- c("Sampling","Tissue")
colnames(lines_gillW)[1:2] <- c("Sampling","Tissue")
colnames(lines_gillC)[1:2] <- c("Sampling","Tissue")

lines_tuniW$Sampling<-as.character(lines_tuniW$Sampling)
lines_tuniC$Sampling<-as.character(lines_tuniC$Sampling)
lines_gillW$Sampling<-as.character(lines_gillW$Sampling)
lines_gillC$Sampling<-as.character(lines_gillC$Sampling)

lines_tuniW<-melt(lines_tuniW)
lines_tuniC<-melt(lines_tuniC)
lines_gillW<-melt(lines_gillW)
lines_gillC<-melt(lines_gillC)

lines_tuniW$Comp_ASV <- paste(lines_tuniW$Tissue, lines_tuniW$variable,  sep ="_")
lines_tuniC$Comp_ASV <- paste(lines_tuniC$Tissue, lines_tuniC$variable,  sep ="_")
lines_gillW$Comp_ASV <- paste(lines_gillW$Tissue, lines_gillW$variable,  sep ="_")
lines_gillC$Comp_ASV <- paste(lines_gillC$Tissue, lines_gillC$variable,  sep ="_")

lines_tuniW <- lines_tuniW[,c(1,5,4)]
lines_tuniC <- lines_tuniC[,c(1,5,4)]
lines_gillW <- lines_gillW[,c(1,5,4)]
lines_gillC <- lines_gillC[,c(1,5,4)]

lines_tuniW <- dcast(data = lines_tuniW,formula = Sampling~Comp_ASV,fun.aggregate = sum,value.var = "value")
lines_tuniC <- dcast(data = lines_tuniC,formula = Sampling~Comp_ASV,fun.aggregate = sum,value.var = "value")
lines_gillW <- dcast(data = lines_gillW,formula = Sampling~Comp_ASV,fun.aggregate = sum,value.var = "value")
lines_gillC <- dcast(data = lines_gillC,formula = Sampling~Comp_ASV,fun.aggregate = sum,value.var = "value")


#### plot environmental concentrations #### 
library(ggplot2)
library(scales)
pdf("tmp_sal_temporal_sd.pdf",width=6, height = 2)
ggplot(codes_plot_e_mean, aes(x=as.numeric(SAMPLING), y=value, fill=variable))+
geom_bar(stat="identity", position= "dodge")+
  scale_fill_manual(values = c("indianred","steelblue"))+
  geom_errorbar(aes(ymin = value -  desvest_E,
                    ymax = value +  desvest_E), 
                stat="identity", position=position_dodge(0.9), width=0.5, size=0.7)+
  scale_y_continuous(name="ºC / Salinity")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("TE_temporal_bla.pdf",width=6, height = 2)
ggplot(codes_plot_p_vil, aes(x=as.numeric(SAMPLING), y=value, fill=variable))+
  geom_bar(stat="identity", position= "dodge")+
  scale_fill_manual(values=c("#ff5400","#ff8e00","#ffd200","#81e650","#00d267","#00c0ff","#8b48fe","#ca41fc","#ff46fb"))+
  scale_y_continuous(name="TE concentration (ppb)", trans=scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 10, 100,1000,10000), labels = scales::trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("TE_temporal_VIL.pdf",width=6, height = 2)
ggplot(codes_plot_p_bla, aes(x=as.numeric(SAMPLING), y=value, fill=variable))+
  geom_bar(stat="identity", position= "dodge")+
  scale_fill_manual(values=c("#ff5400","#ff8e00","#ffd200","#81e650","#00d267","#00c0ff","#8b48fe","#ca41fc","#ff46fb"))+
  scale_y_continuous(name="TE concentration (ppb)", trans=scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 10, 100,1000,10000), labels = scales::trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("TE_temporal_sd.pdf",width=6, height = 2)
ggplot(codes_plot_p_mean, aes(x=as.numeric(SAMPLING), y=value, fill=variable))+
  geom_bar(stat="identity", position= "dodge")+
  scale_fill_manual(values=c("#ff5400","#ff8e00","#ffd200","#81e650","#00d267","#00c0ff","#8b48fe","#ca41fc","#ff46fb"))+
  geom_errorbar(aes(ymin = value -  desvest_P,
                    ymax = value +  desvest_P), 
                stat="identity", position=position_dodge(0.9), width=0.5, size=0.7)+
  scale_y_continuous(name="TE concentration (ppb)", trans=scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 10, 100,1000,10000), labels = scales::trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

#### plot indicator ####  
codes_plot_env_indic <- codes_plot_e_mean[1:8,1:3]
background <- data.frame(xstart = seq(0.5,7.5,1), xend = seq(1.5,8.5,1), col = letters[1:8])

pdf("tunic_W_indic_64_2.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_tuniW, aes(x=as.numeric(Sampling), y=TUNIC_ASV64), color="black", lwd=1) +
  geom_line(data=lines_tuniW, aes(x=as.numeric(Sampling), y=WATER_ASV64), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV64") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_W_indic_7.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_tuniW, aes(x=as.numeric(Sampling), y=TUNIC_ASV7), color="#000000", lwd=1) +
  geom_line(data=lines_tuniW, aes(x=as.numeric(Sampling), y=WATER_ASV7), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV7") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_C_indic_17.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=TUNIC_ASV17), color="#000000", lwd=1) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=WATER_ASV17), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV17") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_C_indic_37.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=TUNIC_ASV37), color="#000000", lwd=1) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=WATER_ASV37), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV37") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_C_indic_52.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=TUNIC_ASV52), color="#000000", lwd=1) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=WATER_ASV52), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV52") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_C_indic_772.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=TUNIC_ASV772), color="#000000", lwd=1) +
  geom_line(data=lines_tuniC, aes(x=as.numeric(Sampling), y=WATER_ASV772), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV772") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_W_indic_11.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=GILL_ASV11), color="#000000", lwd=1) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=WATER_ASV11), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV11") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_W_indic_117.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=GILL_ASV117), color="#000000", lwd=1) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=WATER_ASV117), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV117") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_W_indic_374.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=GILL_ASV374), color="#000000", lwd=1) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=WATER_ASV374), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV374") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_W_indic_50.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=GILL_ASV50), color="#000000", lwd=1) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=WATER_ASV50), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV50") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_W_indic_98.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=GILL_ASV98), color="#000000", lwd=1) +
  geom_line(data=lines_gillW, aes(x=as.numeric(Sampling), y=WATER_ASV98), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV98") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_C_indic_37.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=GILL_ASV37), color="#000000", lwd=1) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=WATER_ASV37), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV37") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_C_indic_4.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=GILL_ASV4), color="#000000", lwd=1) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=WATER_ASV4), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV4") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_C_indic_7.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=GILL_ASV7), color="#000000", lwd=1) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=WATER_ASV7), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV7") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_C_indic_8.pdf",width=3, height = 2.5)
ggplot() +
  geom_rect(data = background, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = c("#b05446","#b05446","#9dcfdd","#9dcfdd","#b05446","#b05446","#9dcfdd","#9dcfdd")), alpha = 0.7) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=GILL_ASV8), color="#000000", lwd=1) +
  geom_line(data=lines_gillC, aes(x=as.numeric(Sampling), y=WATER_ASV8), color="#000000", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="ASV relative abundance")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV8") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

#### plot environmental gill ####  
pdf("gill_env_0.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="SALINITY"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("steelblue"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV0*200), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV0*200), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Salinity (ppm)",sec.axis = sec_axis(~./200, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV0") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
)
dev.off()

pdf("gill_env_1.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="TEMPERATURE"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("indianred"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV1*100), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV1*100), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Temperature (ºC)",sec.axis = sec_axis(~./100, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV1") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_2.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="SALINITY"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("steelblue"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV2*100), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV2*100), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Salinity (ppm)",sec.axis = sec_axis(~./100, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV2") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_3.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="SALINITY"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("steelblue"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV3*100), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV3*100), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Salinity (ppm)",sec.axis = sec_axis(~./100, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV3") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_4.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="SALINITY"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("steelblue"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV4*500), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV4*500), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Salinity (ppm)",sec.axis = sec_axis(~./500, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV4") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_5.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="TEMPERATURE"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("indianred"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV5*300), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV5*300), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Temperature (ºC)",sec.axis = sec_axis(~./300, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV5") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_6.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="TEMPERATURE"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("indianred"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV6*300), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV6*300), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Temperature (ºC)",sec.axis = sec_axis(~./300, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV6") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_9.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="SALINITY"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("steelblue"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV9*400), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV9*400), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Salinity (ppm)",sec.axis = sec_axis(~./400, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV9") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_13.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="TEMPERATURE"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("indianred"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV13*1000), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV13*1000), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Temperature (ºC)",sec.axis = sec_axis(~./1000, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV13") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_14.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="TEMPERATURE"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("indianred"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV14*1000), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV14*1000), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Temperature (ºC)",sec.axis = sec_axis(~./1000, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV14") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_35.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="SALINITY"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("steelblue"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV35*1000), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV35*1000), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Salinity (ppm)",sec.axis = sec_axis(~./1000, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV35") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_env_40.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_e_mean[which(codes_plot_e_mean$variable=="TEMPERATURE"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values = c("indianred"))+
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=GILL_ASV40*1000), color="black", lwd=1) +
  geom_line(data=lines_envG, aes(x=as.numeric(Sampling), y=WATER_ASV40*1000), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="Temperature (ºC)",sec.axis = sec_axis(~./1000, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV40") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()


#### plot pollut gill #### 

pdf("gill_pol_0.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Fe"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#00d267"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV0*3), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV0*3), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./3, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV0") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_1.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Se"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ca41fc"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV1), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV1), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~., name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV1") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_2.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="B"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ff46fb"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV2*2.5), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV2*2.5), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./2.5, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV2") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', colour = "grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_3.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="B"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ff46fb"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV3*2.5), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV3*2.5), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./2.5, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV3") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_4.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Al"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ff5400"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV4*5), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV4*5), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./5, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV4") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_5.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV5/5), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV5/5), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*5, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV5") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_6.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="V"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ffd200"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV6), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV6), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~., name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV6") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_9.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV9/2), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV9/2), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*2, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV9") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_10.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="B"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ff46fb"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV10*50), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV10*50), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./50, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV10") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_11.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="As"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ff8e00"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV11/2), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV11/2), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*2, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV11") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_13.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Zn"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#00c0ff"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV13*30), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV13*30), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./30, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV13") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_14.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Cu"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#8b48fe"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV14*2), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV14*2), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./2, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV14") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_22.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV22*4), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV22*4), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./4, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV22") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_25.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Zn"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#00c0ff"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV25*100), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV25*100), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./100, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV25") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_35.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV35), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV35), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~., name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV35") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_40.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV40), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV40), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~., name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV40") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_50.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Zn"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#00c0ff"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV50*50), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV50*50), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./50, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV50") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_51.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV51), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV51), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~., name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV51") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("gill_pol_95.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Zn"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#00c0ff"))+
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=GILL_ASV95*50), color="black", lwd=1) +
  geom_line(data=lines_polG, aes(x=as.numeric(Sampling), y=WATER_ASV95*50), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./50, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV95") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

#### plot pollut tunic #### 
pdf("tunic_pol_0.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="V"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ffd200"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV0/10), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV0/10), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*10, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV0") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_1.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV1/5), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV1/5), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*5, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV1") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_2.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV2/2), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV2/2), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*2, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV2") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_3.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV3/2), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV3/2), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*2, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV3") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_4.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="V"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ffd200"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV4*1), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV4*1), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./1, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV4") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_5.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV5/5), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV5/5), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*5, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV5") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_9.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Pb"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#81e650"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV9/2), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV9/2), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~.*2, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV9") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_10.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="Zn"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#00c0ff"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV10*30), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV10*30), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./30, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV10") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_23.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="V"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ffd200"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV23*5), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV23*5), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./5, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV23") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

pdf("tunic_pol_108.pdf",width=3, height = 2.5)
ggplot() +
  geom_bar(data=codes_plot_p_mean[which(codes_plot_p_mean$variable=="B"),], aes(x=as.numeric(SAMPLING), y=value, fill=variable), 
           stat="identity", position="dodge")+
  scale_fill_manual(values=c("#ff46fb"))+
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=TUNIC_ASV108*50), color="black", lwd=1) +
  geom_line(data=lines_polT, aes(x=as.numeric(Sampling), y=WATER_ASV108*50), color="black", lwd=0.5, linetype="solid") +
  scale_y_continuous(name="TE concentration (ppm)", sec.axis = sec_axis(~./50, name=" ASV relative abundance"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  ggtitle("ASV108") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",colour = "white", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', color="grey"), 
        panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid')
  )
dev.off()

