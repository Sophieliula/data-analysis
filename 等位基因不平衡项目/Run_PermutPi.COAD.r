Tumor = "COAD"
load("AllelicImbalance/new/COAD/COAD_LLinput_filter_X.RData")

X <- get("X")
purity <- get("purity")
load("AllelicImbalance/new/Cytoband_hg19.RData")
load("AllelicImbalance/new/snp144_chrom_arm_info.RData")
Cytoband <- get("Cytoband")
arm <- get("arm")

X <- data.frame(X, arm = arm[match(X$Rsid, arm[,1]),2])
X$arm <- as.character(X$arm)
X <- X[which(!is.na(X$arm)),]

arm_by_sample <- unique(paste(X$subject, X$arm, sep = "_"))
random_arm_by_sample <- sample(0:1, size = length(arm_by_sample), replace = T)
names(random_arm_by_sample) <- arm_by_sample

Permut.X <- X
Permut.X$arm_bysample <- paste(Permut.X$subject, Permut.X$arm, sep = "_")
Permut.X$random_arm_by_sample <- random_arm_by_sample[Permut.X$arm_bysample]

Permut.X$y0[which(Permut.X$random_arm_by_sample == 0)] <-  X$y1[which(Permut.X$random_arm_by_sample == 0)]
Permut.X$y1[which(Permut.X$random_arm_by_sample == 0)] <-  X$y0[which(Permut.X$random_arm_by_sample == 0)]
Permut.X$yt0[which(Permut.X$random_arm_by_sample == 0)] <-  X$yt1[which(Permut.X$random_arm_by_sample == 0)]
Permut.X$yt1[which(Permut.X$random_arm_by_sample == 0)] <-  X$yt0[which(Permut.X$random_arm_by_sample == 0)]

Permut.X <- Permut.X[,-c(11:13)]
X <- Permut.X

source("AllelicImbalance/scripts/AI/VarAI.1.r")
source("AllelicImbalance/scripts/AI/Ki_Ki0.r")


LL_SNP(unique(X$Rsid), "AllelicImbalance/new/COAD/COAD_SNP_VarAI.permutation.new4")
