load("FourTumor_SCNA_Pi.RData")

Gain <- rbind(brca.gain, coad.gain, luad.gain, prad.gain)
p = ggplot(Gain)+geom_boxplot(aes(x=Tumor, y=Pi,fill=Group)) +
  labs(x = "Tumor", y = "Absolute value of Pi-0.5") +
  scale_fill_manual(values = alpha(c("#ef6548", "#f7fbff"),alpha = 0.8))


Loss <- rbind(brca.loss, coad.loss, luad.loss, prad.loss)
p = ggplot(Loss)+geom_boxplot(aes(x=Tumor, y=Pi,fill=Group)) +
  labs(x = "Tumor", y = "Absolute value of Pi-0.5") +
  scale_fill_manual(values = alpha(c("green", "#f7fbff"),alpha = 0.8))



#ef6548
#41ab5d