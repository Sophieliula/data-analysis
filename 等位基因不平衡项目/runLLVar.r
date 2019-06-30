
Tumor = "BRCA"
load("AllelicImbalance/new/BRCA/BRCA_LLinput_filter_X.RData")
X <- get("X")

temp1 <- sample(as.vector(t(X[,7:10])))
Random.x <- matrix(temp1, ncol = 4)
colnames(Random.x) <- c("y0","y1","yt0","yt1")
X <- data.frame(X[,1:6], Random.x)


purity <- get("purity")
source("AllelicImbalance/scripts/AI/VarAI.r")
source("AllelicImbalance/scripts/AI/Ki_Ki0.r")

#load("AllelicImbalance/new/PRAD/PRAD_Ki_Ki0.RData")
#Ki <- get("Ki")
#Ki0 <- get("Ki0")
#temp <- read.delim("AllelicImbalance/new/PRAD/PRAD_SNP_VarAI.4.txt", head = F, stringsAsFactor = F)
#snp <- unique(X$Rsid)
#unfinish <- snp[-match(temp[,2], snp)]

#LL_SNP(unique(X$Rsid), "AllelicImbalance/new/BRCA/BRCA_SNP_VarAI.permutation.2")
LL_SNP(unique(X$Rsid), "AllelicImbalance/new/BRCA/BRCA_SNP_VarAI.permutation.4")
