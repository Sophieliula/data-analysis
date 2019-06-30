#setwd("E://work/进行中/CPTAC/CPTAC_mutational_landscape/")
library(readxl)
library(reshape2)
library(ComplexHeatmap)


clin <- as.data.frame(read_excel("随访数据20181228 - zhw 修改(1).xlsx"))
clin <- clin[,c("sample","Age","TNM grouped","Gender","Tumor thrombus","Medical History of Liver Cirrhosis")]
purity <- as.data.frame(read_excel("Table tumor purity.xlsx"))
purity$Tumor <- paste0("T",purity$Tumor,sep="")
names(purity)[1:2] <- c("sample","Purity")
names(clin)[6] <- "Liver cirrhosis"
names(clin)[3] <- "TNM stage"
names(clin)[5] <- "Tumor thrombus"
clin <- merge(clin,purity[,1:2],by = "sample")
clin$Age <- cut(clin$Age,c(19,40,61,82))
clin$Age <- gsub(",","-",clin$Age)
clin$Purity <- cut(clin$Purity,c(0.4,0.6,0.8,1))
rownames(clin) <- clin[,1]
clin <- clin[,-1]


mutation <- read.delim("CPTAC_159_mutational_type_filter.IGV.txt",header = F, stringsAsFactors = F)
mutation <- mutation[which(mutation$V1 %in% rownames(clin)), ]
landmatrix <- dcast(mutation,V1 ~ V2, fun.aggregate = function(x){paste(x,collapse = ":")})


rownames(landmatrix) <- landmatrix[,1]
landmatrix <- t(as.matrix(landmatrix[,-1]))
com1 <- intersect(rownames(clin),colnames(landmatrix))
landmatrix <- landmatrix[,match(com1,colnames(landmatrix))]
diffdat <- data.frame(matrix("",13,20))
colnames(diffdat) = setdiff(rownames(clin),mutation$V1)
rownames(diffdat) = rownames(landmatrix)
landmatrix <- cbind(landmatrix,diffdat)


nosy <- read.table("CPTAC_159_sys_nosys.xls", header = T, stringsAsFactors = F)
rownames(nosy) <- nosy[,1]
nosy <- nosy[,-1]
nosy <- nosy[,c(2,1)]


sampleid <- read.table("sampleid.txt", header = F, stringsAsFactors = F)


landmatrix <- landmatrix[,match(sampleid$V1,colnames(landmatrix))]
nosy <- nosy[match(sampleid$V1,rownames(nosy)),]
clin <- clin[match(sampleid$V1, rownames(clin)),]
clin <- clin[,c(1,6,2:5)]
clin$`Tumor thrombus` <- ifelse(clin$`Tumor thrombus`==0,"No","Yes")
clin$`Tumor thrombus` <- factor(clin$`Tumor thrombus`,levels = c("Yes","No"))
clin$`Liver cirrhosis` <- ifelse(clin$`Liver cirrhosis`==0,"No","Yes")
clin$`Liver cirrhosis` <- factor(clin$`Liver cirrhosis`,levels = c("Yes","No"))
clin$Gender <- factor(clin$Gender)

cnv <- read.table("20181230_CPTAC_CNV/CNV_0.8.xls", header = T, stringsAsFactors = F, sep = "\t")
rownames(cnv) <- cnv$protein
cnv1 <- cnv[2:5,]
cnv1[cnv1==1] = "CNA Loss"
cnv2 <- cnv[1,]
cnv2[cnv2==1] = "CNA Gain"
cnv = rbind(cnv1,cnv2)
cnv[cnv==0] = ""
rownames(cnv) = c("CDKN2A del", "RB1 del" , "AXIN1 del" , "TP53 del", "CCND1 amp")
cnv <- cnv[,-1]
cnv <- cnv[,match(sampleid$V1,colnames(cnv))]
cnv_count <- read.table("cnv_count_0.8.txt", header = F, stringsAsFactors = F)
cnv_count$Freq <- cnv_count$V3/159*100
cnv_count$Var1 <- c("CDKN2A del", "RB1 del" , "AXIN1 del" , "TP53 del", "CCND1 amp")
landmatrix <- rbind(landmatrix,cnv)

r_orders = c("TP53", "CTNNB1","AXIN1","ARID1A", "APOB", "ALB", "KMT2C", "KEAP1","TSC2","RB1","ATRX","CPS1", "SMARCA2","CDKN2A del", "RB1 del" , "AXIN1 del" , "TP53 del", "CCND1 amp")

bardat <- unique(mutation[,1:2])
bardat <- data.frame(table(bardat$V2)/159*100)
bardat <- rbind(bardat,cnv_count[,c(5,4)])
bardat$Var1 <- as.character(bardat$Var1)

bardat <- bardat[match(rownames(landmatrix),bardat$Var1),]



oncocol <- c("frameshift deletion" = "#EFC000FF",
             "frameshift insertion" = "#006400",
             "nonframeshift deletion" = "#7AA6DCFF",
             "nonframeshift substitution" = "#EEEE00",
             "nonsynonymous SNV" = "#1F77B4FF",
             "splicing" = "#F4A460",
             "stopgain" = "#8E388E",
             "CNA Loss" = "#00008B",
             "CNA Gain" = "#B22222")




annoTMB <- HeatmapAnnotation(barplot=anno_barplot(nosy, baseline = 0, bar_width = 0.70,
                                              gp = gpar(fill = c("#36648B","#3CB371"),col=NA),
                                              axis = T,axis_side = 'left', border = F),
                             show_annotation_name = F,
                             annotation_name_gp = gpar(fontsize=8),
                             annotation_height = unit.c(unit(4, 'cm'), unit(0.2, 'cm')),
                             gap = unit(0.2, 'cm'), na_col = 'grey')

annorow <- rowAnnotation(barplot = row_anno_barplot(bardat$Freq, baseline = 0, ylim = c(0,100),
                         gp = gpar(fill = "#104E8B", col=NA),
                         axis = TRUE, axis_side = "top", border = F),
                         show_annotation_name = F,
                         annotation_name_gp = gpar(fontsize=8),
                         gap = unit(0.2, 'cm'), na_col = 'grey',
                         width = unit(4, "cm"))


alter_fun = function(x,y,w,h,v){
  n=sum(v)
  h = h*0.72
  w = w*0.80
  grid.rect(x, y, w, h, gp = gpar(fill = "grey90", col = NA))
  if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")}




annoother <- HeatmapAnnotation(df=clin,
                               col = list(Gender = c(female = 'black', male = 'grey90'),
                                          Age = c(`(19-40]` = '#B0E2FF', `(40-61]` = '#32CD32',
                                                  `(61-82]` = '#228B22'),
                                          `Tumor thrombus` = c(`Yes` = 'black', `No` = 'grey90'),
                                          `TNM stage` = c(`I` = '#7AA6DCFF',`II` = '#1F77B4FF',
                                                          `III+IV` = '#CD5555'),
                                          Purity = c(`(0.4,0.6]` = '#9ad7ff', `(0.6,0.8]` = '#818dff',
                                                     `(0.8,1]` = '#a567ff'),
                                          `Liver cirrhosis` = c(`Yes` = 'black',`No` = 'grey90')),
                               show_annotation_name = T, annotation_name_side = 'left',
                               annotation_legend_param = list(
                                 title_gp=gpar(fontsize=10),
                                 labels_gp=gpar(fontsize=10)
                               ),
                               annotation_name_gp = gpar(fontsize=10),
                               gap = unit(0.1, 'cm'), na_col = 'grey90')



pdf("mutational_landscape_20190104.pdf", width = 13, height = 6)
onco <- oncoPrint(landmatrix, get_type = function(x) unique(strsplit(x, ":")[[1]]),
                  alter_fun = alter_fun, col = oncocol,
                  heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=9,
                                              title_gp=gpar(fontsize=10),
                                              labels_gp=gpar(fontsize=10)),
                  column_names_gp = gpar(fontsize = 9),
                  row_names_gp = gpar(fontsize = 10),
                  row_names_side='left',
                  row_title_side='left',
                  row_title_gp=gpar(fontsize=10),
                  top_annotation = annoTMB,
                  bottom_annotation = annoother,
                  column_order = rownames(clin),
                  row_order = r_orders,
                  show_pct = T,
                  column_title_side='top',
                  column_title_gp=gpar(col='white'),
                  show_column_names = F,
                  show_row_barplot = FALSE,
                  pct_gp=gpar(fontsize=10, col='black'))

onco <- onco + annorow
draw(onco, heatmap_legend_side='right')
dev.off()

