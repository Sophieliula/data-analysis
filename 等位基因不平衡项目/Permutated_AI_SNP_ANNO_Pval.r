annovar_temp <- data.frame(anno[permu.pi$rsid,], Class = vector(length = nrow(permu.pi)))
annovar_temp$Class <- "NonAI"
annovar_temp$Class[match(permu.ai.pi$rsid, annovar_temp$snp144)] <- "AI"

ns <- c("stopgain;stopgain","nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV","stopgain","stoploss")


c("LRT_pred", "Polyphen2_HVAR_pred", "SIFT_pred","CADD_phred")

##LRT_pred=====================================================
a <- annovar_temp[which(annovar_temp$LRT_pred != "."),]
Fit <- glm(as.numeric(a$LRT_pred == "D") ~ I(a$Class=="AI"),family = binomial("logit"),na.action = na.omit)

Table <- table(as.numeric(a$LRT_pred == "D") , I(a$Class=="AI"))
p <- summary(Fit)$coefficient[2, 4]
or <- exp(coef(Fit))[2]
hd.p <- phyper(Table[2,2], sum(Table[2,]), sum(Table[1,]), sum(Table[,2]), lower.tail = FALSE, log.p = FALSE)
Table <- as.vector(Table)
anno_result <- c(date(),Tumor, Pi.percent,"LRT",Table,signif(Table[2]/c(Table[1] + Table[2]),3),signif(Table[4]/c(Table[3] + Table[4]),3), signif(or,3), signif(p,3), signif(hd.p,3))

write.table(t(anno_result), file = AI_SNP_ANNO_file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)

#Polyphen======================================================

a <- annovar_temp[which(annovar_temp$Polyphen2_HVAR_pred != "."),]
Fit <- glm(as.numeric(a$Polyphen2_HVAR_pred == "D") ~ I(a$Class=="AI"),family = binomial("logit"),na.action = na.omit)

Table <- table(as.numeric(a$Polyphen2_HVAR_pred == "D") , I(a$Class=="AI"))
p <- summary(Fit)$coefficient[2, 4]
or <- exp(coef(Fit))[2]
hd.p <- phyper(Table[2,2], sum(Table[2,]), sum(Table[1,]), sum(Table[,2]), lower.tail = FALSE, log.p = FALSE)
Table <- as.vector(Table)
anno_result <- c(date(),Tumor,Pi.percent,"Polyphen2_HVAR_pred",Table,signif(Table[2]/c(Table[1] + Table[2]),3),signif(Table[4]/c(Table[3] + Table[4]),3), signif(or,3), signif(p,3), signif(hd.p,3))

write.table(t(anno_result), file = AI_SNP_ANNO_file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)

###SIFT========================================================
a <- annovar_temp[which(annovar_temp$SIFT_pred != "."),]
Fit <- glm(as.numeric(a$SIFT_pred == "D") ~ I(a$Class=="AI"),family = binomial("logit"),na.action = na.omit)

Table <- table(as.numeric(a$SIFT_pred == "D") , I(a$Class=="AI"))
p <- summary(Fit)$coefficient[2, 4]
or <- exp(coef(Fit))[2]
hd.p <- phyper(Table[2,2], sum(Table[2,]), sum(Table[1,]), sum(Table[,2]), lower.tail = FALSE, log.p = FALSE)
Table <- as.vector(Table)
anno_result <- c(date(),Tumor,Pi.percent,"SIFT",Table,signif(Table[2]/c(Table[1] + Table[2]),3),signif(Table[4]/c(Table[3] + Table[4]),3), signif(or,3), signif(p,3), signif(hd.p,3))

write.table(t(anno_result), file = AI_SNP_ANNO_file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)

##CADD=========================================================

a <- annovar_temp[which(annovar_temp$CADD_phred != "."),]
a$CADD_phred <- as.numeric(a$CADD_phred)
Fit <- glm(as.numeric(a$CADD_phred > 15) ~ I(a$Class=="AI"),family = binomial("logit"),na.action = na.omit)

Table <- table(as.numeric(a$CADD_phred > 15) , I(a$Class=="AI"))
p <- summary(Fit)$coefficient[2, 4]
or <- exp(coef(Fit))[2]
hd.p <- phyper(Table[2,2], sum(Table[2,]), sum(Table[1,]), sum(Table[,2]), lower.tail = FALSE, log.p = FALSE)
Table <- as.vector(Table)
anno_result <- c(date(),Tumor,Pi.percent,"CADD",Table,signif(Table[2]/c(Table[1] + Table[2]),3),signif(Table[4]/c(Table[3] + Table[4]),3), signif(or,3), signif(p,3), signif(hd.p,3))

write.table(t(anno_result), file = AI_SNP_ANNO_file, row.names = F, col.names = F, quote = F, sep = "\t", append = T)


