draw_landscape <- function(variant,width,height,res,filename) {
  library(ComplexHeatmap)
 
  wd = as.numeric(width)
  ht = as.numeric(height)
  res = as.numeric(res)
  fn = filename
  variant = read.delim(variant,header = T,stringsAsFactors = F)
  variant = variant[which(variant$REPORT_OR_VUS=="report"), ]
 

  msi = unique(variant[,c('ORDER_ID','MSI')])
  msi[msi$MSI=="","MSI"] = "Unsure"
  msi[! grepl("MSS|MSI-H",msi$MSI),"MSI"] = "Unsure"


  land <- variant[ ,c("ORDER_ID", "GENE", "VAR_TYPE_SX")]
  land <- unique(land)


  freq = land[,c("ORDER_ID","GENE")]
  freq = unique(freq)
  freq = as.data.frame(table(freq$GENE)/length(unique(freq$ORDER_ID))*100, stringsAsFactors = F)
  freq = freq[order(freq$Freq,decreasing = T),]

  if(nrow(freq) > 30){
    genes = freq$Var1[1:30]
    freq = freq[match(genes,freq$Var1), ]
  } else{
    genes = freq$Var1
  }

  TMB = unique(variant[,c("ORDER_ID","TMB")])
  TMB$TMB[grepl('[A-Za-z]',TMB$TMB)] = 0
  TMB$TMB[is.na(TMB$TMB)] = 0
  TMB$TMB = as.numeric(TMB$TMB)
  TMB = TMB[order(TMB$TMB, decreasing = T), ]
  TMB = TMB[! duplicated(TMB$ORDER_ID),]
  
  msi = msi[match(TMB$ORDER_ID,msi$ORDER_ID),]
  maxTMB = ceiling(max(TMB$TMB,na.rm = T))
  freq_max = ceiling(max(freq$Freq))

  landmatrix = reshape2::dcast(land,ORDER_ID ~ GENE, fun.aggre = function(x){paste(x,collapse = ":")})
  rownames(landmatrix) = landmatrix$ORDER_ID


  if(length(genes)==1){
    landmatrix = landmatrix[match(TMB$ORDER_ID,landmatrix$ORDER_ID),]
    landmat = t(landmatrix[,-1])
    colnames(landmat) = landmatrix$ORDER_ID
    rownames(landmat) = colnames(landmatrix)[2:ncol(landmatrix)]
  } else{
    landmat = t(landmatrix[,-1])
    landmat = landmat[match(genes,rownames(landmat)), ] 
    landmat = landmat[,match(TMB$ORDER_ID,colnames(landmat))]
  }



  oncocol <- c("Substitution/Indel" = "#228B22",        ##green
             "Gene Amplification" = "#EE0000",        ##red
             "Gene Homozygous Deletion" = "#0000EE",  ##blue
             "Fusion/Rearrangement" = "#EEEE00",      ##yellow
             "Truncation" = "#8E388E")                ##purple


  alter_fun = function(x,y,w,h,v){
    n=sum(v)
    h = h*0.72
    w = w*0.70
    grid.rect(x, y, w, h, gp = gpar(fill = "grey90", col = NA))
    if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
  }



  top_annotation = HeatmapAnnotation(TMB=anno_barplot(TMB$TMB, height = unit(3, "cm"),
                                                    gp = gpar(fill = '#5A4A75', col=NA),
                                                    which = 'column',
                                                    gap = unit(0.2, 'cm'), na_col = 'grey',
                                                    border = F,
                                                    axis_param = list(at = c(0,maxTMB),
                                                                      labels = c(0, maxTMB))),
                                   MSI=msi$MSI,
                                   annotation_name_gp = gpar(fontsize=10),
                                   col = list(MSI = c(`MSI-H` = "#EE0000", MSS = "#228B22", Unsure = "#EE9A49")),
                                   show_annotation_name = T, annotation_name_side = 'left',
                                   show_legend = F)




  annorow <- rowAnnotation('Mutation \npercentage(%)' = anno_barplot(freq$Freq, ylim = c(0,freq_max),
                                                    gp = gpar(fill = "#104E8B", col=NA),
                                                    border = F,
                                                    axis_param = list(side = "top")),
                         annotation_name_gp = gpar(fontsize=9),
                         annotation_name_side = 'top', 
                         gap = unit(0.5, 'cm'), 
                         na_col = 'grey',
                         width = unit(3, "cm"))



  ld1 = Legend(labels = c("MSI-H","MSS","Unsure"), title = "MSI", 
             legend_gp = gpar(fill = c("#EE0000", "#228B22","#EE9A49")),
             title_gp = gpar(fontsize = 10))

  ld2 = Legend(labels = c("Substitution/Indel","Gene Amplification","Gene Homozygous Deletion",
                        "Fusion/Rearrangement","Truncation"), title = "Alterations", 
             legend_gp = gpar(fill = c("#228B22", "#EE0000","#0000EE", "#EEEE00", "#8E388E")),
             title_gp = gpar(fontsize = 10))




  jpeg(file=paste0('/home/dev/Pictures/landscape/',fn,'.jpeg'), width = wd, height=ht,res=res)
  onco = oncoPrint(landmat, get_type = function(x) unique(strsplit(x, ":")[[1]]),
                 alter_fun = alter_fun, col = oncocol,
                 column_names_gp = gpar(fontsize = 9),
                 row_names_gp = gpar(fontsize = 9),
                 row_names_side='left',
                 row_title_side='left',
                 row_title_gp=gpar(fontsize=9),
                 top_annotation = top_annotation,
                 column_order = TMB$ORDER_ID,
                 row_order = freq$Var1,
                 show_column_names = F,
                 show_heatmap_legend = F,
                 show_pct = F,
                 right_annotation = annorow)
  draw(onco, heatmap_legend_side='right',heatmap_legend_list = list(ld1,ld2))
  invisible(dev.off())
}


     



