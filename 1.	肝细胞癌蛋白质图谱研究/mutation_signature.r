library(ggplot2)
library(scales)
setwd("E://")
dat <- read.table("Documents/20181230_CPTAC_mutational/prof96_plot_Mutational_Profile.txt", 
                  header = T,stringsAsFactors = F)
d = data.frame(x1=c(0.2,16.2,32.2,48.2,64.2,80.2),x2=c(15.8,31.8,47.8,63.8,79.8,96),
               y1=rep(0.124,6),y2=rep(0.14,6),r=unique(dat$Substitution_Type))
d1 = data.frame(x=c(0,97),y=c(0,0))

pdf("mutation.pdf", height = 3, width = 14)
ggplot()+
  geom_bar(data=dat,mapping=aes(1:96,V1,fill=Substitution_Type),stat = "identity", width = 0.5) + 
  xlab('') + ylab('Mutation Type\nProbability') +
  geom_rect(data = d,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=r))+
  geom_text(data=d, aes(x=x1+(x2-x1)/2, y=0.15, label=r), size=7, fontface="bold")+
  scale_y_continuous(labels = percent, expand = c(0,0), breaks = seq(0.03,0.12,0.03),limits = c(0,0.16)) +
  scale_x_continuous(breaks = 1:96,labels = dat$Trinucleotide,expand = c(0,0)) +
  theme_bw() + geom_hline(yintercept = 0,size = 1,colour = 'grey')+
  geom_segment(data = d1,aes(x=x,y=y,xend=x,yend=y+0.12),size = 0.50,colour = 'grey90')+
  scale_fill_manual(values = c("DeepSkyBlue","black","grey","Firebrick2","YellowGreen","pink")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 16,face = "bold", colour = "black"),
        axis.text.x.bottom = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 20, face = "bold", colour = "black"),
        legend.position = "",
        axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))
dev.off()




