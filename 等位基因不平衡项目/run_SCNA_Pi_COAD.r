
Tumor = "COAD"
path = paste("AllelicImbalance/new/", Tumor,"/",sep = "")
load("AllelicImbalance/new/COAD/COAD_LLinput_filter_X.RData")
SCNAfile <- "AllelicImbalance/COAD_AI/SCNA/COAD_snp144_SCNA.RData"

#=========
X <- get("X")
purity <- get("purity")
source("AllelicImbalance/scripts/AI/LLinput_filterX_bySCNA.r")

#==================
X <- gain_X
llresult = paste(path, Tumor,"_SNP_VarAI_Gain.1",sep = "")
source("AllelicImbalance/scripts/AI/Ki_Ki0.r")
source("AllelicImbalance/scripts/AI/VarAI.r")
LL_SNP(unique(X$Rsid), llresult)
#===================
X <- loss_X
llresult = paste(path, Tumor,"_SNP_VarAI_loss.1", sep = "")
source("AllelicImbalance/scripts/AI/Ki_Ki0.r")
source("AllelicImbalance/scripts/AI/VarAI.r")
LL_SNP(unique(X$Rsid), llresult)
#==================
X <- neutral_X
llresult = paste(path, Tumor, "_SNP_VarAI_neutral.1", sep = "")
source("AllelicImbalance/scripts/AI/Ki_Ki0.r")
source("AllelicImbalance/scripts/AI/VarAI.r")
LL_SNP(unique(X$Rsid), llresult)