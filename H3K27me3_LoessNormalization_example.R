setwd("/Users/mgonzalez/Research/01_PhD/01_ChIPseq/05_LoessNormalization")

library(IDPmisc)
library(affy)
library(MASS)
library(ggplot2)
library(grid)
library(viridis)

# Load data
c01<-read.table("03_ChIPlevelsWholeGenome/mm10_2000_bin_H3K27me3_chr19_mESC_recoverChIPlevels/PEAKsignal_mm10_2000_bin_H3K27me3_chr19_mESC.bed")
c02<-read.table("03_ChIPlevelsWholeGenome/mm10_2000_bin_H3K27me3_chr19_MES_recoverChIPlevels/PEAKsignal_mm10_2000_bin_H3K27me3_chr19_MES.bed")
c03<-read.table("03_ChIPlevelsWholeGenome/mm10_2000_bin_H3K27me3_chr19_CP_recoverChIPlevels/PEAKsignal_mm10_2000_bin_H3K27me3_chr19_CP.bed")
c04<-read.table("03_ChIPlevelsWholeGenome/mm10_2000_bin_H3K27me3_chr19_CM_recoverChIPlevels/PEAKsignal_mm10_2000_bin_H3K27me3_chr19_CM.bed")
c05<-read.table("03_ChIPlevelsWholeGenome/mm10_2000_bin_H3K27me3_chr19_NPC_recoverChIPlevels/PEAKsignal_mm10_2000_bin_H3K27me3_chr19_NPC.bed")
c06<-read.table("03_ChIPlevelsWholeGenome/mm10_2000_bin_H3K27me3_chr19_CN_recoverChIPlevels/PEAKsignal_mm10_2000_bin_H3K27me3_chr19_CN.bed")

c07<-read.table("04_ChIPlevelsRegions/mESC_9_BP_clean_H3K27me3_chr19_mESC_recoverChIPlevels/PEAKsignal_mESC_9_BP_clean_H3K27me3_chr19_mESC.bed")
c08<-read.table("04_ChIPlevelsRegions/mESC_9_BP_clean_H3K27me3_chr19_MES_recoverChIPlevels/PEAKsignal_mESC_9_BP_clean_H3K27me3_chr19_MES.bed")
c09<-read.table("04_ChIPlevelsRegions/mESC_9_BP_clean_H3K27me3_chr19_CP_recoverChIPlevels/PEAKsignal_mESC_9_BP_clean_H3K27me3_chr19_CP.bed")
c10<-read.table("04_ChIPlevelsRegions/mESC_9_BP_clean_H3K27me3_chr19_CM_recoverChIPlevels/PEAKsignal_mESC_9_BP_clean_H3K27me3_chr19_CM.bed")
c11<-read.table("04_ChIPlevelsRegions/mESC_9_BP_clean_H3K27me3_chr19_NPC_recoverChIPlevels/PEAKsignal_mESC_9_BP_clean_H3K27me3_chr19_NPC.bed")
c12<-read.table("04_ChIPlevelsRegions/mESC_9_BP_clean_H3K27me3_chr19_CN_recoverChIPlevels/PEAKsignal_mESC_9_BP_clean_H3K27me3_chr19_CN.bed")

c13<-read.table("04_ChIPlevelsRegions/mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_mESC_recoverChIPlevels/PEAKsignal_mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_mESC.bed")
c14<-read.table("04_ChIPlevelsRegions/mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_MES_recoverChIPlevels/PEAKsignal_mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_MES.bed")
c15<-read.table("04_ChIPlevelsRegions/mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_CP_recoverChIPlevels/PEAKsignal_mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_CP.bed")
c16<-read.table("04_ChIPlevelsRegions/mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_CM_recoverChIPlevels/PEAKsignal_mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_CM.bed")
c17<-read.table("04_ChIPlevelsRegions/mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_NPC_recoverChIPlevels/PEAKsignal_mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_NPC.bed")
c18<-read.table("04_ChIPlevelsRegions/mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_CN_recoverChIPlevels/PEAKsignal_mESC_9_PE_clean_HiC_FDR0_P100_H3K27me3_chr19_CN.bed")


# Obtain data frame for whole genome
H3K27me3_ESC<-c01[,5]
H3K27me3_MES<-c02[,5]
H3K27me3_CP<-c03[,5]
H3K27me3_CM<-c04[,5]
H3K27me3_NPC<-c05[,5]
H3K27me3_CN<-c06[,5]

regions<-paste(c01[,1],c01[,2],c01[,3])
H3K27me3<-data.frame(H3K27me3_ESC,H3K27me3_MES,H3K27me3_CP,H3K27me3_CM,H3K27me3_NPC,H3K27me3_CN)
row.names(H3K27me3)<-regions
H3K27me3<-H3K27me3[H3K27me3_ESC!=0 & H3K27me3_MES!=0 & H3K27me3_CP!=0 & H3K27me3_CM!=0 & H3K27me3_NPC!=0 & H3K27me3_CN!=0,]
H3K27me3<-NaRV.omit(H3K27me3)
m1<-as.matrix(H3K27me3)
m1<-NaRV.omit(m1)
row.m1<-nrow(m1)


# Obtain data frame for promoters
H3K27me3_ESC<-c07[,5]
H3K27me3_MES<-c08[,5]
H3K27me3_CP<-c09[,5]
H3K27me3_CM<-c10[,5]
H3K27me3_NPC<-c11[,5]
H3K27me3_CN<-c12[,5]

regions<-paste(c12[,1],c12[,2],c12[,3])
H3K27me3<-data.frame(H3K27me3_ESC,H3K27me3_MES,H3K27me3_CP,H3K27me3_CM,H3K27me3_NPC,H3K27me3_CN)
row.names(H3K27me3)<-regions
H3K27me3<-H3K27me3[H3K27me3_ESC!=0 & H3K27me3_MES!=0 & H3K27me3_CP!=0 & H3K27me3_CM!=0 & H3K27me3_NPC!=0 & H3K27me3_CN!=0,]
H3K27me3<-NaRV.omit(H3K27me3)
m2<-as.matrix(H3K27me3)
m2<-NaRV.omit(m2)
row.m2<-nrow(m2)

# Obtain data frame for enhancers
H3K27me3_ESC<-c13[,5]
H3K27me3_MES<-c14[,5]
H3K27me3_CP<-c15[,5]
H3K27me3_CM<-c16[,5]
H3K27me3_NPC<-c17[,5]
H3K27me3_CN<-c18[,5]

regions<-paste(c13[,1],c13[,2],c13[,3])
H3K27me3<-data.frame(H3K27me3_ESC,H3K27me3_MES,H3K27me3_CP,H3K27me3_CM,H3K27me3_NPC,H3K27me3_CN)
row.names(H3K27me3)<-regions
H3K27me3<-H3K27me3[H3K27me3_ESC!=0 & H3K27me3_MES!=0 & H3K27me3_CP!=0 & H3K27me3_CM!=0 & H3K27me3_NPC!=0 & H3K27me3_CN!=0,]
H3K27me3<-NaRV.omit(H3K27me3)
m3<-as.matrix(H3K27me3)
m3<-NaRV.omit(m3)
row.m3<-nrow(m3)

# Normalization
m<-rbind(m1,m2,m3)
row.m<-nrow(m)
s <- seq(1,row.m1)
#s <- row.names(m1)
set.seed(123)
mn <- normalize.loess(m,subset=s)

mn1<-mn[s,]
p<-seq(row.m1+1,row.m1+row.m2)
mn2<-mn[p,]
e<-seq(row.m1+row.m2+1,row.m)
mn3<-mn[e,]

write.table(mn1,"NormFiles/normH3K27me3_chr19_bins.txt")
write.table(mn2,"NormFiles/normH3K27me3_chr19_BP.txt")
write.table(mn3,"NormFiles/normH3K27me3_chr19_PE.txt")


# MA plot before norm ggplot

pdf("plots/MAplot_H3K27me3_chr19_HiC_FDR0_P100_beforeNorm-ggplot.pdf",width = 15, height = 10)
M<-log2(c((m1[,2]+0.1)/(m1[,1]+0.1),(m1[,3]+0.1)/(m1[,1]+0.1),(m1[,4]+0.1)/(m1[,1]+0.1),(m1[,5]+0.1)/(m1[,1]+0.1),(m1[,6]+0.1)/(m1[,1]+0.1)))
A<-log2(c((m1[,2]+0.1)*(m1[,1]+0.1)/2,(m1[,3]+0.1)*(m1[,1]+0.1)/2,(m1[,4]+0.1)*(m1[,1]+0.1)/2,(m1[,5]+0.1)*(m1[,1]+0.1)/2,(m1[,6]+0.1)*(m1[,1]+0.1)/2))
Cell<-c(rep("MES",nrow(m1)),rep("CP",nrow(m1)),rep("CM",nrow(m1)),rep("NPC",nrow(m1)),rep("CN",nrow(m1)))
c<-data.frame(M,A,Cell)
c$Cell = factor(c$Cell, levels = c("MES","CP","CM","NPC","CN"))
ggplot(c) +
  geom_hex(aes(A, M), bins = 30) +
  scale_fill_gradientn("", colours = rev(viridis(300)))+
  geom_smooth(aes(A, M),method = "loess", level=0.5)+
  geom_hline(yintercept = 0,linetype="dashed")+
  labs(title="H3K27me3 MA plot before normalization",x="A", y = "M") +
  theme_bw() +
  theme(legend.position="right",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  facet_wrap(~ Cell, scales="free")
dev.off()



# MA plot after norm ggplot

pdf("plots/MAplot_H3K27me3_chr19_HiC_FDR0_P100_afterNorm-ggplot.pdf",width = 15, height = 10)
M<-round(log2(c((mn1[,2]+0.1)/(mn1[,1]+0.1),(mn1[,3]+0.1)/(mn1[,1]+0.1),(mn1[,4]+0.1)/(mn1[,1]+0.1),(mn1[,5]+0.1)/(mn1[,1]+0.1),(mn1[,6]+0.1)/(mn1[,1]+0.1))),2)
A<-round(log2(c((mn1[,2]+0.1)*(mn1[,1]+0.1)/2,(mn1[,3]+0.1)*(mn1[,1]+0.1)/2,(mn1[,4]+0.1)*(mn1[,1]+0.1)/2,(mn1[,5]+0.1)*(mn1[,1]+0.1)/2,(mn1[,6]+0.1)*(mn1[,1]+0.1)/2)),2)
Cell<-c(rep("MES",nrow(mn1)),rep("CP",nrow(mn1)),rep("CM",nrow(mn1)),rep("NPC",nrow(mn1)),rep("CN",nrow(mn1)))
c<-data.frame(M,A,Cell)
c$Cell = factor(c$Cell, levels = c("MES","CP","CM","NPC","CN"))
ggplot(c) +
  geom_hex(aes(A, M), bins = 30) +
  scale_fill_gradientn("", colours = rev(viridis(300)))+
  geom_smooth(aes(A, M),method = "loess", level=0.5)+
  geom_hline(yintercept = 0,linetype="dashed")+
  labs(title="H3K27me3 MA plot after normalization",x="A", y = "M") +
  theme_bw() +
  theme(legend.position="right",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  facet_wrap(~ Cell, scales="free")
dev.off()




# Boxplots

# Bins
Levels<-c(log2(m1[,1]),log2(m1[,2]),log2(m1[,3]),log2(m1[,4]),log2(m1[,5]),log2(m1[,6]),
          log2(mn1[,1]),log2(mn1[,2]),log2(mn1[,3]),log2(mn1[,4]),log2(mn1[,5]),log2(mn1[,6]))
Cell<-c(rep("ESC",nrow(m1)),rep("MES",nrow(m1)),rep("CP",nrow(m1)),rep("CM",nrow(m1)),rep("NPC",nrow(m1)),rep("CN",nrow(m1)),
        rep("ESC",nrow(mn1)),rep("MES",nrow(mn1)),rep("CP",nrow(mn1)),rep("CM",nrow(mn1)),rep("NPC",nrow(mn1)),rep("CN",nrow(mn1)))
Norm<-c(rep("Before normalization",nrow(m1)*6),rep("After normalization",nrow(mn1)*6))
c<-data.frame(Levels,Cell,Norm)
c$Cell = factor(c$Cell, levels = c("ESC","MES","CP","CM","NPC","CN"))
c$Norm = factor(c$Norm, levels = c("Before normalization","After normalization"))

pdf("plots/Histones_normH3K27me3_chr19_bins.pdf",width = 12, height = 10)
ggplot(c, aes(x=Cell, y=Levels,fill=Cell)) + 
  #geom_violin(position=position_dodge(1),lwd=1) +
  scale_fill_manual(values=plasma(n=6)) +
  geom_boxplot(width=0.4,position=position_dodge(1),outlier.size=-1,show.legend = FALSE,lwd=1,colour="black") +
  labs(title="Histone mark levels before and after normalization (bins)",x="", y = "log2(Histone mark + 0.1)") +
  theme(legend.position="bottom",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  facet_wrap(~ Norm, scales="free")
dev.off()

# Promoters
Levels<-c(log2(m2[,1]),log2(m2[,2]),log2(m2[,3]),log2(m2[,4]),log2(m2[,5]),log2(m2[,6]),
          log2(mn2[,1]),log2(mn2[,2]),log2(mn2[,3]),log2(mn2[,4]),log2(mn2[,5]),log2(mn2[,6]))
Cell<-c(rep("ESC",nrow(m2)),rep("MES",nrow(m2)),rep("CP",nrow(m2)),rep("CM",nrow(m2)),rep("NPC",nrow(m2)),rep("CN",nrow(m2)),
        rep("ESC",nrow(mn2)),rep("MES",nrow(mn2)),rep("CP",nrow(mn2)),rep("CM",nrow(mn2)),rep("NPC",nrow(mn2)),rep("CN",nrow(mn2)))
Norm<-c(rep("Before normalization",nrow(m2)*6),rep("After normalization",nrow(mn2)*6))
c<-data.frame(Levels,Cell,Norm)
c$Cell = factor(c$Cell, levels = c("ESC","MES","CP","CM","NPC","CN"))
c$Norm = factor(c$Norm, levels = c("Before normalization","After normalization"))

pdf("plots/Histones_normH3K27me3_chr19_BP.pdf",width = 12, height = 10)
ggplot(c, aes(x=Cell, y=Levels,fill=Cell)) + 
  #geom_violin(position=position_dodge(1),lwd=1) +
  scale_fill_manual(values=plasma(n=6)) +
  geom_boxplot(width=0.4,position=position_dodge(1),outlier.size=-1,show.legend = FALSE,lwd=1,colour="black") +
  labs(title="Histone mark levels before and after normalization (BP)",x="", y = "log2(Histone mark + 0.1)") +
  theme(legend.position="bottom",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  facet_wrap(~ Norm, scales="free")
dev.off()

# Enhancers
Levels<-c(log2(m3[,1]),log2(m3[,2]),log2(m3[,3]),log2(m3[,4]),log2(m3[,5]),log2(m3[,6]),
          log2(mn3[,1]),log2(mn3[,2]),log2(mn3[,3]),log2(mn3[,4]),log2(mn3[,5]),log2(mn3[,6]))
Cell<-c(rep("ESC",nrow(m3)),rep("MES",nrow(m3)),rep("CP",nrow(m3)),rep("CM",nrow(m3)),rep("NPC",nrow(m3)),rep("CN",nrow(m3)),
        rep("ESC",nrow(mn3)),rep("MES",nrow(mn3)),rep("CP",nrow(mn3)),rep("CM",nrow(mn3)),rep("NPC",nrow(mn3)),rep("CN",nrow(mn3)))
Norm<-c(rep("Before normalization",nrow(m3)*6),rep("After normalization",nrow(mn3)*6))
c<-data.frame(Levels,Cell,Norm)
c$Cell = factor(c$Cell, levels = c("ESC","MES","CP","CM","NPC","CN"))
c$Norm = factor(c$Norm, levels = c("Before normalization","After normalization"))

pdf("plots/Histones_normH3K27me3_chr19_PE.pdf",width = 12, height = 10)
ggplot(c, aes(x=Cell, y=Levels,fill=Cell)) + 
  #geom_violin(position=position_dodge(1),lwd=1) +
  scale_fill_manual(values=plasma(n=6)) +
  geom_boxplot(width=0.4,position=position_dodge(1),outlier.size=-1,show.legend = FALSE,lwd=1,colour="black") +
  labs(title="Histone mark levels before and after normalization (PE)",x="", y = "log2(Histone mark + 0.1)") +
  theme(legend.position="bottom",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  facet_wrap(~ Norm, scales="free")
dev.off()

