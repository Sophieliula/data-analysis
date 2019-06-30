library(ggplot2)
library(dplyr)

TMB_value <- read.csv("neoantigen_TMB/TMB.csv",header = F,stringsAsFactors = F,quote = "")
colnames(TMB_value) <- c("sample","TMB")

result_type <- c("nonsy_results","nonsy+sy_results","nonsy+sy+utr_results")

for(type in result_type){
  
  
  AA_df <- read.table(paste("AA+ABF1_results/",type,"/AA.xls",sep = ""),header = T,stringsAsFactors = F,sep = "\t")
  AA_df$group <- "AA"
  AA_df$group[which(AA_df$q..0.05==0)]="Non-AA"
  
  AA_df <- AA_df %>% left_join(TMB_value)
  wilcox_test <- wilcox.test(AA_df$TMB[which(AA_df$group=="AA")],AA_df$TMB[which(AA_df$group=="Non-AA")])
  
  
  
  # boxplot
  p <- ggplot(AA_df ,aes(x=group, y=TMB)) +geom_boxplot(fill=c('#FA8072', '#63B8FF'),width=0.4,outlier.colour = "#363636", outlier.size = 1.1) + geom_jitter(width = 0.03,size=1.1,color="#363636")
  p <- p + xlab("") + ylab("TMB (mutation per Mb)")
  if(wilcox_test$p.value < 0.001){
    p <- p + annotate(x=1.5, y=25, label="p < 0.001",geom="text", size=10,color="black")
  }
  p <- p + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(color = "black"))
  p <- p + theme(axis.text=element_text(size=24,color="black"),axis.title=element_text(size=24,color="black")) + theme(plot.title = element_text(size=30,hjust = 0.5))
  p <- p + scale_x_discrete(labels=c("AA" = "AA", "Non-AA" = "Non-AA"))
  ggsave(paste("W://计算机肿瘤/CPTAC/CPTAC_20181229/TMB_boxplot_v3_20181127_",type,".pdf",sep = ""), p,width = 6,height = 6)
  print(median(AA_df$TMB[which(AA_df$group=="AA")],na.rm = T))
  print(median(AA_df$TMB[which(AA_df$group=="Non-AA")],na.rm = T))
}
