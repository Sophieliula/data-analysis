############## Package

library(VennDiagram)
library(ggplot2)


############### Data
load("snp144_15w_ANNO_variant_function.RData")
load("REFSEQ.hg19.RData")
load("snp144_annotation.RData")

############### Read Pi

### Tumor Type, real Pi, permut Pi 
Tumor = "BRCA"
llresult = "BRCA_SNP_VarAI.5.txt"
ll.permut.result = "BRCA_SNP_VarAI.permutation.4times.txt"

### at three confidence levels: 1%, 5%, 10% 

perc = list(c(0.025, 0.975), c(0.05, 0.95), c(0.1, 0.9)) 
pi.perc <- c("5%Pi","10%Pi","20%Pi")

for(i in 1:3){
print(i)
ai.pi.threshold <- perc[[i]]
Pi.percent = pi.perc[i]
## read pi
source("AllelicImbalance/scripts/AI/ReadPi.r")
## AI SNP functional enrichment 
AI_SNP_ANNO_file = "AllelicImbalance/new/Permut_AI_SNP_ANNO.1.txt"
source("AllelicImbalance/scripts/AI/AI_SNP_ANNO_Pval.r")
}

############## AI at different SCNA levels 

load("AllelicImbalance/BRCA_ERpos_AI/SCNA/BRCA_snp144_SCNA.RData")
load("AllelicImbalance/LUAD_AI/SCNA/LUAD_snp144_SCNA.RData")
load("AllelicImbalance/COAD_AI/SCNA/COAD_snp144_SCNA.RData")

library(parallel)
Sample <- intersect(colnames(SCNA), unique(X$subject))
mySCNACount = function(SNP, down = NA, up = NA){
 	ans <- mclapply(SNP, function(x){
		#	print(x)
			Sample <- intersect(X[which(X$Rsid == x),"subject"],colnames(SCNA))
			scna <- na.omit(SCNA[x,Sample])
			snp1 <- which(scna >= up)
			snp2 <- which(scna <= down)

			c(length(snp2), length(snp1), length(snp1) + length(snp2), length(scna))
	}, mc.cores = 35)
	ans
	do.call(rbind, ans)
}

test <- mySCNACount(rownames(SCNA), -0.5, 0.5)
scna_count <- test
colnames(scna_count) <- c("loss", "gain","Sum","SampleNumber")
rownames(scna_count) <- rownames(SCNA)
save(scna_count, file = paste(Tumor,"_scna_count.RData",sep = ""))

for(i in 1:3){
print(i)
ai.pi.threshold <- perc[[i]]
Pi.percent = pi.perc[i]
Permut.AI.Pi <- Pi[which(Pi$pi <= quantile(Permu.Pi$pi, ai.pi.threshold[1]) | Pi$pi >= quantile(Permu.Pi$pi, ai.pi.threshold[2])),]
AI.Pi <- Permut.AI.Pi
source("AllelicImbalance/scripts/AI/SCNA_Enrichment.r")




