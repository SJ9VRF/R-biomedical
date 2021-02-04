#package for analysing biological issues with linear models
library(limma)
#package for getting dataset in our test
library(GEOquery) 
#downloading data
gset <- getGEO("GSE59259", GSEMatrix =TRUE)
#put 8 control samples from non-neoplastic liver tissue and
# 8 samples of hepatocellular carcinoma in each appropriate together.
type <- c(rep("N",8),rep("H",8))
# defining the design matrix for fitting linear model
des <- model.matrix(~ type+0)
# Put expresions of gens in log scale
Data <- log(exprs(gset[[1]])+0.001)
# Fitting linear model 
fit <- lmFit(Data, des) 
# setup bayes on fitted model to moderate the standard errors
fit <- eBayes(fit)
# ADJUSTING P-VALUE
pVal_ad<-p.adjust(fit$p.value, method = p.adjust.methods)
pVal_ad <- as.array(pVal_ad)
colnam <- colnames(fit$p.value)
fit$p.value<- matrix( pVal_ad, nrow = nrow(fit$p.value), ncol=2)
colnames(fit$p.value) <- colnam
# classifyTests for adjusted p-value thereshold
results3 <- decideTests(fit, p.value=1e-17)
results4 <- decideTests(fit, p.value=0.01)
#top 5 de gens
top5DEgens <- topTable(fit, number=5, fit$p.value = 0.01, adjust.method ="BH")
#heatmap plot, save to "heatmap.png"
pheatmap::pheatmap(top5DEgens, filename = "heatmap.png")
# save de gens with adjusted p.value<0.01 to "de.csv"
gen_names <- row.names(Data)
colnam_re <- c("H", "N")
colnames(results4) <- colnam_re
DE <- cbind(gen_names, fit$p.value, results4)
DE <- data.frame(DE)
DE <- subset(DE, N == 1 & H==1, select = c(2,3))
write.csv(DE, file = "de.csv")
