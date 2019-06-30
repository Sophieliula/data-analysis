setwd("W://计算机肿瘤/CPTAC/TCGA-HCC/")
library(ggplot2)
library(reshape2)
variant <- read.table("cbioportal_TCGA_366.txt", header = T, stringsAsFactors = F, na.strings = "NA")
chcc <- read.table("E://work/进行中/CPTAC/script/bardat.txt",header=T,stringsAsFactors = F,sep="\t")
clin <- read.delim("TCGA_CELL_paper_clinic_data.txt", header = T, stringsAsFactors = F)

clin$Barcode <- gsub("A$","",clin$Barcode)
clin$Barcode <- gsub("-",".",clin$Barcode)
variant <- melt(variant,id.vars = c("GENE_ID","COMMON"))
variant <- variant[which(variant$value != "NaN"),]
variant$variable <- as.character(variant$variable)
variant <- variant[,2:3]
variant <- unique(variant)


hbv <- clin$Barcode[which(clin$HBV_consensus=="pos")]
hcv <- clin$Barcode[which(clin$HCV_consensus=="pos")]


tcga_hbv <- variant[which(variant$variable %in% hbv), ]
tcga_hcv <- variant[which(variant$variable %in% hcv), ]
tcga_hbv <- data.frame(table(tcga_hbv$COMMON)/length(intersect(variant$variable,hbv))*100)
tcga_hcv <- data.frame(table(tcga_hcv$COMMON)/length(intersect(variant$variable,hcv))*100)
tcga_hbv$Var1 <- as.character(tcga_hbv$Var1)
tcga_hcv$Var1 <- as.character(tcga_hcv$Var1)

length(intersect(variant$variable,hbv)) #37
length(intersect(variant$variable,hcv)) #31

tcga_hbv <- tcga_hbv[which(tcga_hbv$Var1 %in% chcc$Var1),]
tcga_hcv <- tcga_hcv[which(tcga_hcv$Var1 %in% chcc$Var1),]



plotdat <- rbind(chcc,tcga_hbv,tcga_hcv)
plotdat$Type <- rep(c("CHCC-HBV (n=159)","TCGA-HBV (n=37)","TCGA-HCV (n=31)"), times=c(nrow(chcc),nrow(tcga_hbv),nrow(tcga_hcv)))
plotdat$Var1 <- factor(plotdat$Var1, levels = c("TP53","CTNNB1","AXIN1","ARID1A","APOB","ALB","KMT2C","KEAP1","TSC2","RB1","ATRX","CPS1","SMARCA2"))

dat <- rep(c("KMT2C","TSC2","ATRX","CPS1","SMARCA2"), times=2)
dat <- data.frame(dat)
dat$Freq <- 0
dat$Type <- rep(c("TCGA-HBV (n=37)","TCGA-HCV (n=31)"),times=c(5,5))
names(dat)[1] <- "Var1"
plotdat <- rbind(plotdat,dat)


ggplot(plotdat,aes(Var1,Freq,fill=Type)) +
  geom_bar(stat = "identity", position=position_dodge(), width = 0.7) +
  xlab("") +  ylab("Mutation Frequency(%)") +
  scale_fill_manual(values=c("CHCC-HBV (n=159)"="#CD3333","TCGA-HBV (n=37)"="#FF6A6A","TCGA-HCV (n=31)"="#00BFFF")) +
  theme_bw() + scale_y_continuous(expand = c(0,0), limits =c(0,80))+
  theme(panel.grid=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color="black", size=12),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black",size = 12),
        legend.position = c(0.9,0.8),
        axis.title.y = element_text(color="black", size=13))


ggsave('W://计算机肿瘤/CPTAC/CPTAC_20181229/TCGA_barplot.pdf',width = 12, height = 4)
