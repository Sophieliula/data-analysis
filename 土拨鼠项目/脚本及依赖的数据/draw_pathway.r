draw_pathway <- function(variant,width,height,res,filename,pathway) {
  library(ComplexHeatmap)
  library(dplyr)
  wd = as.numeric(width)
  ht = as.numeric(height)
  res = as.numeric(res)
  fn = filename
  pw = pathway
  variant = read.delim(variant,header = T,stringsAsFactors = F)
  variant = variant[which(variant$REPORT_OR_VUS=="report"), ]
  load("/home/dev/Documents/tuboshu/pathway_info.RData")
  land <- variant[ ,c("ORDER_ID", "GENE", "VAR_TYPE_SX")]
  land <- unique(land)
  
  freq = land[,c("ORDER_ID","GENE")]
  freq = unique(freq)
  freq = as.data.frame(table(freq$GENE)/length(unique(freq$ORDER_ID))*100, stringsAsFactors = F)
  colnames(freq) = c("gene","freq")
  
  pw = pathway_name[match(pw,pathway_name$value),"geneGroup"]
  allpw = as.data.frame(pathway_genes[,pw])
  colnames(allpw) = pw
  
  allgene = reshape2::melt(allpw,id = 0)
  allgene = na.omit(allgene)
  colnames(allgene) = c("pathway","gene")
  
  freq = merge(freq,allgene,by="gene",all = F)

  gene_freq = freq %>% group_by(pathway) %>% top_n(n = 5, wt = freq)
  gene_freq = as.data.frame(gene_freq)
  
  landmat = land[which(land$GENE %in% gene_freq$gene),]
  landmat = unique(landmat)
  landmat = reshape2::dcast(landmat,ORDER_ID ~ GENE, fun.aggregate = function(x){paste(x, collapse = ':')})
  rownames(landmat) = landmat$ORDER_ID
  landmat = t(as.matrix(landmat[,-1]))
  
  gene_freq = gene_freq[match(rownames(landmat),gene_freq$gene),]
  gene_freq = gene_freq[order(gene_freq$pathway,gene_freq$freq, decreasing = T),]
  landmat = landmat[match(gene_freq$gene,rownames(landmat)),]
  
  TMB = unique(variant[,c("ORDER_ID","TMB")])
  TMB = TMB[which(TMB$ORDER_ID %in% colnames(landmat)), ]
  TMB$TMB[grepl('[A-Za-z]',TMB$TMB)] = 0
  TMB$TMB[is.na(TMB$TMB)] = 0
  TMB$TMB = as.numeric(TMB$TMB)
  TMB = TMB[order(TMB$TMB, decreasing = T), ]
  TMB = TMB[! duplicated(TMB$ORDER_ID),]
  
  msi = unique(variant[,c('ORDER_ID','MSI')])
  msi[msi$MSI=="","MSI"] = "Unsure"
  msi[! grepl("MSS|MSI-H",msi$MSI),"MSI"] = "Unsure"
  msi = msi[match(TMB$ORDER_ID,msi$ORDER_ID),]
  
  landmat = landmat[,match(TMB$ORDER_ID,colnames(landmat))]
  
  
  maxTMB = ceiling(max(TMB$TMB,na.rm = T))
  

  
  # plot
  
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
  
  
  
  ld1 = Legend(labels = c("MSI-H","MSS","Unsure"), title = "MSI", 
               legend_gp = gpar(fill = c("#EE0000", "#228B22","#EE9A49")),
               title_gp = gpar(fontsize = 10))
  
  ld2 = Legend(labels = c("Substitution/Indel","Gene Amplification","Gene Homozygous Deletion",
                          "Fusion/Rearrangement","Truncation"), title = "Alterations", 
               legend_gp = gpar(fill = c("#228B22", "#EE0000","#0000EE", "#EEEE00", "#8E388E")),
               title_gp = gpar(fontsize = 10))
  
  
  jpeg(file=paste0('/home/dev/Pictures/pathway/',fn,'.jpeg'), width = wd, height=ht,res=res)
  
  onco <- oncoPrint(landmat, get_type = function(x) unique(strsplit(x, ":")[[1]]),
                    alter_fun = alter_fun, col = oncocol,
                    row_names_gp = gpar(fontsize = 10),
                    row_names_side='left',
                    row_title_side='left',
                    row_title_gp=gpar(fontsize=10),
                    column_order = TMB$ORDER_ID,
                    top_annotation = top_annotation,
                    show_column_names = F,
                    show_heatmap_legend = F,
                    show_pct = F,
                    row_split = gene_freq$pathway, gap=unit(1, 'cm'))
  
  draw(onco, heatmap_legend_side='right', heatmap_legend_list = list(ld1,ld2))
  invisible(dev.off())
  
}


