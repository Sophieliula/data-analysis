#===============SCNA Density================
Gain.Pi <- MYPi(gain_pi_file, 10, 0.4, 0.6, 0.1)

Loss.Pi <- MYPi(loss_pi_file, 10, 0.4, 0.6, 0.1)

Neutral.Pi <- MYPi(neutral_pi_file, 10, 0.4, 0.6, 0.1)


plot.multi.dens <- function(s)
{
    junk.x = NULL
    junk.y = NULL
    for(i in 1:length(s)) {
        junk.x = c(junk.x, density(s[[i]])$x)
        junk.y = c(junk.y, density(s[[i]])$y)
    }
    xr <- range(junk.x)
    yr <- range(junk.y)
    plot(density(s[[1]]), xlim = xr, ylim = yr, main = "")
    for(i in 1:length(s)) {
        lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
    }
}

# the input of the following function MUST be a numeric list

pdf(paste(Tumor, "_SNP_Pi_GroupBySCNA_density.1.pdf",sep = ""))
par(pty = "s", las = 1)
plot.multi.dens( list(Neutral.Pi$pi, Gain.Pi$pi, Loss.Pi$pi))
legend("topright",pch=c(-1,-1,-1),lty=c(1,1,1),col=c("black","green","red"),bty = "n",legend=c("Neutral","Loss","Gain"))

dev.off()

