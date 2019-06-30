#####function####

MYPi = function(Path = NA, SampleNumber = NA, phi.down = NA, phi.up = NA, delta = NA){
		Pi <- read.delim(Path, head = F, stringsAsFactor = F)	
		colnames(Pi) <- c("gene","rsid","Sample_Number","theta", "pi", "phi", "delta","convergence","message")
#		rownames(Pi) <- Pi$rsid
		Pi <- data.frame(Pi[,1:2], ANNO[Pi$rsid,c(2,4,5:8,10)], Pi[,3:9])
		n1 <- grep("^OR", Pi$gene)
		if(length(n1) > 0) Pi <- Pi[-n1,]
		n2 <- grep("^HLA", Pi$gene)
		if(length(n2) > 0) Pi <- Pi[-n2,]
		Pi <- Pi[which(Pi[,"convergence"] == 0),]
		Pi <- Pi[which(Pi[,"delta"] < delta & Pi[,"phi"] > phi.down & Pi[,"phi"] < phi.up),]
		Pi <- Pi[which(Pi$Sample_Number >= SampleNumber),]
		Pi
}

PAR_PLOT = function(Pi = NA, Group = NA){
	hist(Pi[,"pi"], br = 80, xlab = "pi", main = Group)
	hist(Pi[,"phi"], br = 80, xlab = "phi", main = Group)
	hist(Pi[,"theta"], br = 80, xlab = "theta", main = Group)
	hist(Pi[,"delta"], br = 80, xlab = "delta", main = Group)
}

####plot ####

#====================  Pi  =================================

Pi <- MYPi(llresult, 10, 0.4, 0.6, 0.1)
exonic.snp <- anno$snp144[which(anno$Func.refGene %in% c("splicing;splicing","splicing","exonic;splicing","exonic;exonic","UTR3", "UTR3;UTR3","UTR5","UTR5;UTR3","exonic"))]
Pi <- Pi[which(Pi$rsid %in% exonic.snp),]
Pi  <-  Pi[order(Pi$pi),]
Pi <- Pi[which(Pi$chrom != "chrX"),]

Permu.Pi <- MYPi(ll.permut.result, 10, 0.4, 0.6, 0.1)
#exonic.snp <- anno$snp144[which(anno$Func.refGene %in% c("splicing;splicing","splicing","exonic;splicing","exonic;exonic","UTR3", "UTR3;UTR3","UTR5","UTR5;UTR3","exonic"))]
Permu.Pi <- Permu.Pi[which(Permu.Pi$rsid %in% exonic.snp),]
Permu.Pi  <-  Permu.Pi[order(Permu.Pi$pi),]
Permu.Pi <- Permu.Pi[which(Permu.Pi$chrom != "chrX"),]



Permut.AI.Pi <- Pi[which(Pi$pi <= quantile(Permu.Pi$pi, ai.pi.threshold[1]) | Pi$pi >= quantile(Permu.Pi$pi, ai.pi.threshold[2])),]

pi_result <- c(date(),Tumor,Pi.percent,nrow(Pi), nrow(Permut.AI.Pi), quantile(Permu.Pi$pi, ai.pi.threshold))

write.table(t(pi_result), file = "AllelicImbalance/new/Permut_AI_SNP_Pi_result.1.txt", row.names = F, col.names = F, quote = F, sep = "\t", append = T)

AI.Pi <- Permut.AI.Pi


#######

#write.table(Pi[,1:14], file = paste(Tumor, "_SNP_Pi.1.csv", sep = ""),  col.names = T,  row.names = F,  quote = F,sep = ",")
#write.table(AI.Pi[,1:14], file = paste(Tumor, "_AISNP_Pi.1.csv", sep = ""), col.names = T, row.names = F, quote = F,sep = ",")

##
write.table(unique(AI.Pi$gene), file = paste(Tumor,gsub("%","percent",Pi.percent),"genelist.1.txt", sep = "_"),col.names = F,  row.names = F,  quote = F)

