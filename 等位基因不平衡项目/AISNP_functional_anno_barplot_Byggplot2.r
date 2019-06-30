library(ggplot2)
library(reshape2)


#anno <- read.delim("Permut_AI_SNP_ANNO.3.6.txt", head = F, stringsAsFactors = F)
anno <- read.delim("Permut_AI_SNP_ANNO.3.6.txt", head = F, stringsAsFactors = F)
anno1 <- anno[,c(2:4,9:11,12)]
#anno1 <- data.frame(anno1,nonAI = anno[,7]/(anno[,5] + anno[,7]), AI = anno[,8]/(anno[,6] + anno[,8]))
colnames(anno1) <- c("Tumor", "Pi.threshold", "ANNO","nonAI","AI","OR","pval")
temp <- data.frame(rbind(anno1[,c(1:3, 6:7)], anno1[,c(1:3, 6:7)]), Fraction = c(anno1[,4], anno1[,5]), Class = c(rep("NonAI",nrow(anno1)), rep("AI",nrow(anno1))))
#temp$Fraction <- signif(temp$Fraction,3)

#temp1 <- temp[which(temp$Tumor == "BRCA" & temp$ANNO == "Polyphen2_HVAR_pred"),]
temp1 <- temp[which(temp$Tumor == "BRCA" & temp$Pi.threshold =="5%Pi"),]

#temp1$Pi.threshold <- as.factor(temp1$Pi.threshold)
temp1$ANNO <- as.factor(temp1$ANNO)

#temp1 <- temp[which(temp$ANNO == "Polyphen2_HVAR_pred"),]

myANNO.plot = function(DB = NA){
temp1 <- temp[which(temp$ANNO == DB),]
p = ggplot(temp1, aes(x=Pi.threshold, y=Fraction, fill=Class)) + 
  geom_bar(width = 0.4, position="dodge", stat="identity") +
  scale_fill_manual(values = alpha(c("#fc4e2a","#4393c3"),alpha = 0.8)) +
  scale_x_discrete(limits = c("20%Pi","10%Pi", "5%Pi")) +
#  geom_text(aes(label=Fraction), vjust= 0.25, hjust = - 0.2, colour="black", position=position_dodge(.9), size=3) +
  ylim(0,max(temp1$Fraction) + 0.01) +
  labs(x = "Pi Threshold", y = "Functional SNP Fraction")
p + coord_flip() + facet_wrap(~Tumor)

}

for(x in c("CADD","SIFT","LRT","Polyphen2_HVAR_pred")){
  pdf(paste(x,"_ANNO.pdf", sep = ""),width = 11)  
  myANNO.plot(x)
  dev.off()
}


pdf("Polyphen2_HVAR_pred_ANNO.pdf")  
p + coord_flip() + facet_grid(Tumor~Tumor)
dev.off()
 

"#4393c3"
c("#fc4e2a", "#4393c3")

myANNO.plot = function(PI = NA){
    temp1 <- temp[which(temp$Pi.threshold == PI),]
  #  temp1 <- temp[which(temp$Tumor == PI),]
  p = ggplot(temp1, aes(x=ANNO, y=Fraction, fill=Class)) + 
    geom_bar(width = 0.4, position="dodge", stat="identity") +
    scale_fill_manual(values = alpha(c("#fc4e2a", "#4393c3"),alpha = 0.8)) +
    scale_x_discrete(limits = c("LRT","Polyphen2_HVAR_pred", "SIFT", "CADD")) +
    #  geom_text(aes(label=Fraction), vjust= 0.25, hjust = - 0.2, colour="black", position=position_dodge(.9), size=3) +
    ylim(0,max(temp1$Fraction) + 0.01) +
    labs(x = "Database", y = "Fraction")
 
  p + coord_flip() + facet_wrap(~Tumor,nrow = 1)
  
#  p + theme_set(theme_bw()) 
#  p +theme(panel.grid.major=element_line(colour=NA))
}

pdf("5percent.pdf")
myANNO.plot("5%Pi")
dev.off()


temp1 <- temp
#  temp1 <- temp[which(temp$Tumor == PI),]
p = ggplot(temp1, aes(x=ANNO, y=Fraction, fill=Class)) + 
  geom_bar(width = 0.4, position="dodge", stat="identity") +
  scale_fill_manual(values = alpha(c("#fc4e2a", "#4393c3"),alpha = 0.8)) +
  scale_x_discrete(limits = c("LRT","Polyphen2_HVAR_pred", "SIFT", "CADD")) +
 # geom_text(aes(label=ANNO), angle = 30) +
  ylim(0,max(temp1$Fraction) + 0.01) +
  labs(x = "Database", y = "Fraction") 

p + coord_flip() + facet_grid(Pi.threshold~Tumor)


