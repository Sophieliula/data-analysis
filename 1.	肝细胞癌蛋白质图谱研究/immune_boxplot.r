setwd("E://project/部门内部科研/20190401/CPTAC/")
library(ggplot2)
load("AA_signature_associated_boxplot.RData")
genes = unique(TMEmarker$ImmunCell)
plot_gene = c("TReg","PDL1","IFN-gama","Eosinophils")
pos = match(plot_gene,genes)

for (i in pos){
  dat = RNA_TMEbox[[i]][[1]]
  aa = genes[i]
  wilcox_test = wilcox.test(dat$expre[dat$AA=="AA"],dat$expre[dat$AA=="Non-AA"])
  p = ggplot(dat ,aes(x=AA, y=expre)) +
    geom_boxplot(fill=c('#FA8072', '#63B8FF'),width=0.4,outlier.colour = "#363636", outlier.size = 1.1) + geom_jitter(width = 0.03,size=1.1, color="#363636") +
    labs(y="Abundance\n", size=8) +
    ggtitle(aa) +annotate(x=1.5, y=max(dat$expre), label= signif(wilcox_test$p.value,2),geom="text", size=6,color="black") +
    scale_x_discrete(labels=c("AA-1","AA-0")) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 18,  color = 'black'),
          legend.title = element_blank(),
          axis.title.y = element_text(size = 18,  color = 'black'),
          legend.position = "",
          plot.title = element_text(hjust = 0.5,size = 21),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black')
    )
  ggsave(paste0("E://project/部门内部科研/20190401/CPTAC/","RNA_RNA_TME_",aa,".pdf"),p, width = 6, height = 6)
  print(p)
}

pos2 = pos[3:4]

for (i in pos2){
  dat = Protein_TMEbox[[i]][[1]]
  aa = genes[i]
  wilcox_test = wilcox.test(dat$expre[dat$AA=="AA"],dat$expre[dat$AA=="Non-AA"])
  p = ggplot(dat ,aes(x=AA, y=expre)) +
    geom_boxplot(fill=c('#FA8072', '#63B8FF'),width=0.4,outlier.colour = "#363636", outlier.size = 1.1) + geom_jitter(width = 0.03,size=1.1, color="#363636") +
    labs(y="Abundance\n", size=8) +
    ggtitle(aa) +annotate(x=1.5, y=max(dat$expre), label= signif(wilcox_test$p.value,2),geom="text", size=6,color="black") +
    scale_x_discrete(labels=c("AA-1","AA-0")) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 18,  color = 'black'),
          legend.title = element_blank(),
          axis.title.y = element_text(size = 18,  color = 'black'),
          legend.position = "",
          plot.title = element_text(hjust = 0.5,size = 21),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black')
    )
  ggsave(paste0("E://project/部门内部科研/20190401/CPTAC/","Protein_RNA_TME_",aa,".pdf"),p, width = 6, height = 6)
  print(p)
}




