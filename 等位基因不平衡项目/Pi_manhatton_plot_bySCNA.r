###=================  Manhatten plot ===================================================


Manhaton.plot.one = function(pi.neu = NA, pi.loss = NA, pi.gain = NA, chr=NA){
		chr = paste("chr",chr,sep = "")
        snp1 <- pi.neu[which(pi.neu[,"chrom"] == chr), "rsid"]
        pos1 <- pi.neu[snp1,"end"]/1E6
        pi1 <- abs(pi.neu[snp1,"pi"] - 0.5)
        plot(0, 0, type = "p", pch = ".", col = "black", xlim = c(min(pos1), max(pos1)),ylim = c(0,0.35), xlab = "Position (mb)", ylab = "Pi", main = chr)

		snp.loss <- pi.loss[which(pi.loss[,"chrom"] == chr), "rsid"]
		if(length(snp.loss) > 0){
			pos.loss <- pi.loss[snp.loss,"end"]/1E6
			pi.loss <- abs(pi.loss[snp.loss,"pi"] - 0.5)
			points(pos.loss, pi.loss, type = "p", pch = ".", col = "green")
		}

		snp.gain <- pi.gain[which(pi.gain[,"chrom"] == chr), "rsid"]
		if(length(snp.gain) > 0){
			pos.gain <- pi.gain[snp.gain,"end"]/1E6
			pi.gain <- abs(pi.gain[snp.gain,"pi"] - 0.5)
			points(pos.gain, pi.gain, type = "p", pch = ".", col = "red")
		}

        points(pos1, pi1, type = "p", pch = ".", col = "grey")

#	abline(h = 0.5, col = "black")
}



pdf(paste(Tumor,"_SNP_Pi_Manhattanplot.1.pdf",sep = ""), paper="a4")

sapply(seq(1,19,3), function(x){
print(x)
par(mfrow = c(3,1), las = 1)
Manhaton.plot.one(Neutral.Pi, Loss.Pi, Gain.Pi, x)
Manhaton.plot.one(Neutral.Pi, Loss.Pi, Gain.Pi, x+1)
Manhaton.plot.one(Neutral.Pi, Loss.Pi, Gain.Pi, x+2)
})

par(mfrow = c(3,1), las = 1)
Manhaton.plot.one(Neutral.Pi, Loss.Pi, Gain.Pi, 22)
dev.off()

