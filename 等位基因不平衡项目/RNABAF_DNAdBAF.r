#========RNA BAF ~ DNA dBAF ===========
load("AllelicImbalance/BRCA_ERpos_AI/ASE/BRCA_AD_Rawdata.RData")
load("AllelicImbalance/BRCA_ERpos_AI/BRCA_Exome_snp_forLL.RData")
X <- Data
X <- X[which(X$y0 >= 10 & X$y1 >= 10),]
X <- X[which(X$subject %in% Purity[,1]),]

AD.dep <- AD[,,1] + AD[,,2]
AD.dep[AD.dep < 20] <- NA
BAF <- AD[,,2]/AD.dep

AI.Pi <- Permut.AI.Pi
AI.Pi <- AI.Pi[order(AI.Pi$pi),]
AISNP.BAF <- BAF[AI.Pi$rsid,]
AISNP.BAF <- AISNP.BAF[,intersect(colnames(AISNP.BAF), unique(X$subject))]
n <- apply(AISNP.BAF, 1, function(x){length(na.omit(x))})
AISNP.BAF <- AISNP.BAF[which(n >0),]

ans <- mclapply(rownames(AISNP.BAF),function(snp){
		print(snp)
		N <- which(X$Rsid == snp)
		SUJ <- X[N,"subject"]
		tbaf <- X[N,"yt1"]/(X[N,"yt0"] + X[N,"yt1"])
		nbaf <- X[N,"y1"]/(X[N,"y0"] + X[N,"y1"])
		dbaf <- tbaf - nbaf
		names(tbaf) <- SUJ
		names(nbaf) <- SUJ
		names(dbaf) <- SUJ
		baf <- na.omit(BAF[snp, ])
		if(length(baf) < 1) return(rep(NA,8))

		SUJ <- intersect(SUJ, names(baf))
		if(length(SUJ) < 1) return(rep(NA,8))

		baf <- baf[SUJ]
		dbaf <- dbaf[SUJ]
				
		dep <- AD[snp, SUJ,]
		if(length(SUJ) > 1){
		cbind(baf, dbaf, dep, nbaf[SUJ], tbaf[SUJ], rep(snp, length(SUJ)),SUJ)
		}else{
		c(baf, dbaf, dep, nbaf[SUJ], tbaf[SUJ], snp,SUJ)
		}

#		info <- cbind(rep(paste("hs",AI.Pi[snp,"chrom"],sep = ""), length(baf)), rep(AI.Pi[snp,"start"], length(baf)), rep(AI.Pi[snp,"end"], length(baf)))
#		write.table(cbind(info, baf), file = "BRCA_RNA_BAF_forCircos.txt",quote = F, col.names = F, row.names = F,append = T)
#		write.table(cbind(info, tbaf[names(baf)]), file = "BRCA_DNA_TBAF_forCircos.txt",quote = F, col.names = F, row.names = F,append = T)

},mc.cores = 30)
dat <- do.call(rbind, ans)
dat <- as.data.frame(dat, stringsAsFactors = F)
dat <- dat[which(!is.na(dat[,7])),]
dat[,1] <- as.numeric(dat[,1])
dat[,2] <- as.numeric(dat[,2])
dat[,3] <- as.numeric(dat[,3])
dat[,4] <- as.numeric(dat[,4])
dat[,5] <- as.numeric(dat[,5])
dat[,6] <- as.numeric(dat[,6])
dat$gene <- ANNO[dat[,7],"gene"]

write.table(cbind(paste("hs", AI.Pi$chrom, sep = ""), AI.Pi$start, AI.Pi$end, AI.Pi$pi), file = "BRCA_AISNP_forCircos.txt", quote = F, col.names = F, row.names = F)

dat.ref <- do.call(rbind, ans)
dat.alt <- do.call(rbind, ans)
summary(dat.ref[which(dat.ref[,1] < 0.1 & dat.ref[,2] > -0.1 & dat.ref[,2] < 0.1),3])
colnames(dat) <- c("rna.baf", "dna.baf")

pdf("BRCA_RNABAF_DNAdBAF_3.22.pdf")
par(pty = "s", las = 1)
smoothScatter(dat[,2], dat[,1], xlab = "DNA dBAF", ylab = "RNA Tumor BAF")
dev.off()

pdf("BRCA_RNABAF_DNAdBAF_scatter_3.22.pdf")
par(pty = "s", las = 1)
plot(dat[,2], dat[,1],pch = ".", xlab = "DNA dBAF", ylab = "RNA Tumor BAF")
abline(h=0.15,col = "red")
dev.off()


temp <- dat[which(dat[,1] < 0.18 & dat[,2] > -0.1 & dat[,2] < 0.1),]
extreme_snp <- names(which(table(temp[,7]) > 5))
temp_table <- data.frame(snp = extreme_snp,Sample_number = table(temp[,7])[extreme_snp],gene = ANNO[extreme_snp,"gene"],func = ANNO[extreme_snp,"func"])
temp_table <- data.frame(temp_table, pp= anno[extreme_snp,"Polyphen2_HDIV_pred"],lrt=anno[extreme_snp,"LRT_pred"],sift = anno[extreme_snp,"SIFT_pred"],CADD=anno[extreme_snp,"CADD_phred"])
temp_table$snp <- as.character(temp_table$snp)
temp_table$gene <- as.character(temp_table$gene)


pdf("test.pdf")
par(mfrow = c(1,3), las = 1)
plot(0,0, xlim = c(0,0.18), ylim = c(-1,1),pch = "",xlab = "RNA BAF",ylab = "DNA dBAF")
points(temp[,1], temp[,2], pch = ".")
plot(0,0, xlim = c(0,0.18), ylim = c(0,1),xlab = "RNA BAF",ylab = "DNA Normal BAF")
points(temp[,1], temp[,5], pch = ".", col = "green")
plot(0,0, xlim = c(0,0.18), ylim = c(0,1), xlab = "RNA BAF",ylab = "DNA Tumor BAF")
points(temp[,1], temp[,6], pch = ".", col = "red")
dev.off()


dat$name <- paste(dat[,7],dat[,8],sep = "_")
rownames(dat) <- dat$name
depth.temp1 <- temp[,3] + temp[,4]
depth.temp2 <- dat[,3] + dat[,4]
depth.temp2 <- depth.temp2[-which(dat$name %in% temp$name)]
depth.temp <- data.frame(depth = c(depth.temp1,depth.temp2),class = c(rep("Abnormal",length(depth.temp1)),rep("Normal",length(depth.temp2))))

pdf("Extreme_value_depth_boxplot.pdf")
par(pty = "s", las = 1)
boxplot(depth.temp[,1] ~ depth.temp[,2], xlab = "Group",ylab = "RNA depth", main = "AI SNP")
dev.off()


pdf("Extreme_value_depth.pdf")
par(mfrow = c(1,2),pty = "s", las = 1)
hist(dat[,3] + dat[,4], br = 40, xlab = "RNA depth", main = "AI SNP")
hist(depth.temp, br = 40, xlab = "RNA depth",ylab = "", main = "Extreme value")
dev.off()

pdf("Extreme_value_depth.pdf")
par(mfrow = c(1,2),pty = "s", las = 1)
boxplot(dat[,3] + dat[,4], br = 40, xlab = "RNA depth", main = "AI SNP")
hist(depth.temp, br = 40, xlab = "RNA depth",ylab = "", main = "Extreme value")
dev.off()

