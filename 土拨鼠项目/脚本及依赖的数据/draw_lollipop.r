draw_lollipop <- function(maf,width,height,res,filename,gene){ 
  wd = as.numeric(width)
  ht = as.numeric(height)
  res = as.numeric(res)
  fn = filename
  gene = gene
  source("/home/dev/Documents/tuboshu/lollipop.r")
  load("/home/dev/Documents/tuboshu/gff.RData")
  maf = read.delim(maf, header = T,stringsAsFactors = F)
  maf = maf[which(maf$REPORT_OR_VUS=="report"),]
  maf = maf[which(maf$GENE == gene), ]
  #refID = maf$TRANSCRIPT_ID[grepl('^[NM_0-9]+$',maf$TRANSCRIPT_ID)]
  refID = sapply(strsplit(maf$TRANSCRIPT_ID,"；"),'[',1)
  refID = names(sort(table(refID),decreasing = T))[1]
  jpeg(file=paste0('/home/dev/Pictures/lollipop/',fn,'.jpeg'), width = wd, height=ht,res=res)
  #dev.control(displaylist = "enable")
  lollipop(maf = maf, gene = gene, gff=gff, AACol = 'Protein_Change', refSeqID = refID)
  #lolplot = recordPlot()
  #invisible(dev.off())
  dev.off()
}

# args <- commandArgs(trailingOnly = TRUE)
# maf <- args[1]
# width <- args[2]
# height <- args[3]
# res <- args[4]
# filename <- args[5]
