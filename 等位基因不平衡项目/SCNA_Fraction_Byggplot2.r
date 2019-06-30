library(ggplot2)
library(reshape2)


Data <- read.delim("AI_SNP_SCNA.txt", head = F, stringsAsFactors = F)
Data1 <- Data[,1:5]
colnames(Data1) <- c("Pi.threshold","Tumor",  "Group","nonAI","AI")
temp <- data.frame(rbind(Data1[,1:3], Data1[,1:3]), Fraction = c(Data1[,4], Data1[,5]), Class = c(rep("NonAI",nrow(Data1)), rep("AI",nrow(Data1))))
#temp$Fraction <- signif(temp$Fraction,3)
temp$Pi.threshold <- factor(temp$Pi.threshold,levels =c("5%Pi","10%Pi","20%Pi"))
temp$Group <- factor(temp$Group, levels = c("30%-100%","20%-30%","10%-20%","0-10%"))

#temp1 <- temp[which(temp$Tumor == Tumor),]
temp1 <- temp
p = ggplot(temp1, aes(x=Group, y=Fraction, fill = Class)) + 
  geom_bar(width = 0.5, position="dodge", stat="identity") +
  scale_fill_manual(values = alpha(c("#ec7014","#225ea8"),alpha = 0.8)) +
#  scale_fill_brewer(palette="Pastel1") +
  labs(x = "", y = "")

p + coord_flip() + facet_grid(Tumor~Pi.threshold)

pdf(paste(Tumor,"_SCNA_Enrichment.pdf",sep = ""))  
p + coord_flip() + facet_grid(Pi.threshold~.)
dev.off()

