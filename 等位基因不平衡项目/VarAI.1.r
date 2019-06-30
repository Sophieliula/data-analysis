library(emdbook)
library(parallel)

LL.varai <- function(par, YY, Ki, Ki0, H, H0, L, purity, track = FALSE) {
	if (track) cat(par, "\t")
	theta <- par[1]   
	pi0   <- 0.5
	pi1   <- par[2]      
	phi   <- par[3]  
	rho   <- purity
	delta <- par[4]  
	kk <- 2.0 
	ww <- 0.2
	sigma2 <- 1E4
	subj <- unique(YY[, "subject"])
	N <- nrow(YY)
	pr <- -dbeta(pi1, shape1 = 10, shape2 = 10, log = TRUE) -
	dbeta(delta, shape1 = 1.01, shape2 = 1.99, log = TRUE) - 
	dgamma(theta, shape = kk / 2, rate = ww / 2, log = TRUE) -
	dbeta(phi, shape1 = 1, shape2 = 1, log = TRUE)
	lli <- sapply(subj, function(i) {
	indx <- which(YY[, "subject"] == i & L == 1)
	if (length(indx) == 0) return(c(-999))
	Ymat <- as.matrix(YY[indx, -c(1:4)])
	l <- L[indx]
	h <- H[indx]
	h0 <- H0[indx]
	pi_ <- rho[i] * pi1 + (1 - rho[i]) * pi0
	QD <- matrix(c(1 - pi_, pi_), nr = 1)
	QDadj <- c(2 * (1 - phi) * ((1 - delta) * QD[1] + delta * QD[2]), 2 * phi * ((1 - delta) * QD[2] + delta * QD[1]))
	K_il <- cbind(h * Ki[i] * QDadj[1], h * Ki[i] * QDadj[2])
	Kappa_il <- K_il[, 1] + K_il[, 2]
	prob <-  K_il[, 1] / Kappa_il
	pl <- sum(sapply(1:length(Kappa_il), function(j) {
	dbetabinom(x = Ymat[j, "yt0"], size = Ymat[j, "yt0"] + Ymat[j, "yt1"], prob = prob[j], theta = theta, log = TRUE )
	}))
	if (all(c("y0", "y1") %in% colnames(YY)) & all(!is.na(Ki0))) {
	QD <- c(pi0, 1 - pi0)
	QDadj <- c(2 * (1 - phi) * ((1 - delta) * QD[1] + delta * QD[2]), 2 * phi * ((1 - delta) * QD[2] + delta * QD[1]))
	K_il <- cbind(h0 * Ki0[i] * QDadj[1], h0 * Ki0[i] * QDadj[2])
	Kappa_il <- K_il[, 1] + K_il[, 2]
	prob <-  K_il[, 1] / Kappa_il
	p0 <- sum(sapply(1:length(Kappa_il), function(j) {
	dbetabinom(x = Ymat[j, "y0"], size = Ymat[j, "y0"] + Ymat[j, "y1"], prob = prob[j], theta = theta, log = TRUE )
	}))
	return(-p0 -pl)
	}
	return(-pl)
	})
	liklihd <- sum(lli, na.rm = T) + pr
	if (track) cat(liklihd, "\n")
	return(liklihd)
}

LL_SNP <- function(SNP,filename){
			mclapply(SNP, function(g) {
			cat(g, "\n")

			Y <- X[which(X[, "Rsid"] == g & 
				!is.na(X[, "yt0"]) & 
				!is.na(X[, "yt1"])), ]
			sid <- unique(Y[, "subject"])
			vid <- unique(Y[, "Rsid"])
			N <- length(sid)
			Y <- Y[order(Y[, "i"]), ]
### get L, for het L = 1; for homo L = 2 ###		
			L <- rep(1, nrow(Y))

			H <- rep(1, length(sid))

			H0 <- rep(1, length(sid))
			
		m =	optim(par = c(0.3, 0.4, 0.4, 0.1), fn = LL.varai, 
					YY = Y[, -c(2:3)], Ki = Ki[sid], 
					Ki0 = Ki0[sid], H = H, H0 = H0, L = L,
					 purity = purity[sid], track = F,
			method = "L-BFGS-B",
			lower = rep(1E-8, 4),
			upper = c(Inf, 1 - 1E-8, 1 - 1E-8, 0.1))
			output <- c(INFO[g,1],g, N, m$par, m$convergence, m$message)
			write.table(t(output), file = paste(filename,"txt",sep = "."), quote = F, col.names = F, row.names = F, append = T, sep = "\t")
},mc.cores = 32)
}


INFO <- X[,2:3]
INFO <- INFO[which(!duplicated(INFO[,2])),]
rownames(INFO) <- INFO[,2]



