library(dplyr)
library(reshape2)
install.packages("PairedData")
library(PairedData)
setwd("G:/styela/Molecular/Microbioma/mirobioma_temporal")
TE <- read.delim("metadata_pollutants_aigua.txt", header = TRUE, sep = "\t", dec = ".")

TE <- TE[,c(2,4,7:15)]
TE$SAMPLING <- as.character(TE$SAMPLING)

TE_ <- melt(TE)

pdf("paired_test.pdf", width=4.5, height=4)
ggplot(TE_, aes(x = POP, y = value)) + 
 # geom_boxplot(aes(fill = POP), alpha = .2) +
  geom_line(aes(group = SAMPLING, color=SAMPLING), alpha=0.5) + 
  scale_color_manual(values=c("#2980B9", "#27AE60", "#C0392B", "#F1C40F","#2980B9", "#27AE60", "#C0392B", "#F1C40F"))+
  geom_point(size = 2.5, alpha=0.5, aes(color=SAMPLING), pch=4, stroke=1.2) + 
  facet_wrap(~ variable, scales="free")+
  theme(legend.position = "null")
dev.off()

vils <- TE_[which(TE_$POP=="VIL"),]
blas <- TE_[which(TE_$POP=="BLA"),]

vil <- as.vector(vils$value)
bla <- as.vector(blas$value)

t.test(vil, bla, paired = TRUE, alternative="two.sided")


vilsAl <- as.vector(vils[which(vils$variable=="Al"),])
vilsAs <- as.vector(vils[which(vils$variable=="As"),])
vilsV <- as.vector(vils[which(vils$variable=="V"),])
vilsPb <- as.vector(vils[which(vils$variable=="Pb"),])
vilsFe <- as.vector(vils[which(vils$variable=="Fe"),])
vilsZn <- as.vector(vils[which(vils$variable=="Zn"),])
vilsCu <- as.vector(vils[which(vils$variable=="Cu"),])
vilsSe <- as.vector(vils[which(vils$variable=="Se"),])
vilsB <- as.vector(vils[which(vils$variable=="B"),])

blasAl <- as.vector(blas[which(blas$variable=="Al"),])
blasAs <- as.vector(blas[which(blas$variable=="As"),])
blasV <- as.vector(blas[which(blas$variable=="V"),])
blasPb <- as.vector(blas[which(blas$variable=="Pb"),])
blasFe <- as.vector(blas[which(blas$variable=="Fe"),])
blasZn <- as.vector(blas[which(blas$variable=="Zn"),])
blasCu <- as.vector(blas[which(blas$variable=="Cu"),])
blasSe <- as.vector(blas[which(blas$variable=="Se"),])
blasB <- as.vector(blas[which(blas$variable=="B"),])


t.test(vilsAl$value, blasAl$value, paired = TRUE)
t.test(vilsAs$value, blasAs$value, paired = TRUE)
t.test(vilsV$value, blasV$value, paired = TRUE)
t.test(vilsPb$value, blasPb$value, paired = TRUE)
t.test(vilsFe$value, blasFe$value, paired = TRUE)
t.test(vilsZn$value, blasZn$value, paired = TRUE)
t.test(vilsCu$value, blasCu$value, paired = TRUE)
t.test(vilsSe$value, blasSe$value, paired = TRUE)
t.test(vilsB$value, blasB$value, paired = TRUE)
