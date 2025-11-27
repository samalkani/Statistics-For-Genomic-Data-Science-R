# Batch effects in R - Part A

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "sva", "bladderbatch", "snpStats"
# BiocManager::install(c("sva", "bladderbatch", "snpStats"), type = "source", force = TRUE)


# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load Libraries
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

# 5. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 6. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 7. Making the circles on the plots filled solid
par(pch=19)

# 8. Load "bladderdata"
data(bladderdata)
ls()

# 9. Set up the data
pheno = pData(bladderEset)
head(pheno, n = 3)
edata = exprs(bladderEset)
head(edata, n = 3)

# 10. Dimensions of expression object
dim(edata)

# 11. Dimensions for the phenotype object
dim(pheno)

# 12. Information about phenotype variables
head(pheno)

# 13. Adjusting for batch effects with a linear model

# 13. A. Adjusting for batch effects with a linear model – Model 1
mod = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
head(mod)
fit = lm.fit(mod,t(edata))
hist(fit$coefficients[2,],col=2,breaks=100)

# 13. B. Checking for correlation between the outcome variable and the batch effect – Model 1
table(pheno$cancer,pheno$batch)

# 14. Adjusting for batch effects using Combat

# 14. A. Adjusting for batch effects using Combat – Model 2
batch = pheno$batch

# 14. B. Construct model matrix with intercept term, and pheno data
modcombat = model.matrix(~1, data=pheno)

# 14. C. Add the Cancer (outcome) term to the model
modcancer = model.matrix(~cancer, data=pheno)

# 14. D. Run combat and clean expression data add in batch variable, plots
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# 14. E. Fit the adjusted model on the cleaned combat data
combat_fit = lm.fit(modcancer,t(combat_edata))

# 14. F. Plot histogram of coefficients from Combat fit
hist(combat_fit$coefficients[2,],col=2,breaks=100)

# 14. G. Comparing the Combat model with the linear adjustment
plot(fit$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="Linear Model",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)

# 15. Adjusting for batch effects, when the batch variable is unknown using SVA

# 15. A. Model with known adjustment variables + variable we care about (cancer) - Model 3
mod = model.matrix(~cancer, data=pheno)

# 15. B. Model with unknown adjustment variables to discover batch effect
mod0 = model.matrix(~1, data=pheno)

# 15. C. Create some covariates (2 surrogate batch effects)
sva1 = sva(edata,mod,mod0,n.sv=2)

# 15. D. Display new covariates that have been created
names(sva1)

# 15. E. Dimensions
dim(sva1$sv)

# 15. F. Correlate the potential batch effects with the actual observed batch variable
summary(lm(sva1$sv ~ pheno$batch))

# 15. G. Boxplot comparing surrogate batch variable with observed batch variable
boxplot(sva1$sv[,2] ~ pheno$batch)

# 15. H. Boxplot comparing surrogate batch variable with observed batch variable
boxplot(sva1$sv[,2] ~ pheno$batch)
points(sva1$sv[,2] ~ jitter(as.numeric(pheno$batch)),col=as.numeric(pheno$batch))

# 16. Construct adjusted model with surrogate batch effects and observed batch effects

# 16. A. Construct adjusted model with surrogate batch effects and observed batch effects – Model 4
modsv = cbind(mod,sva1$sv)
head(modsv)

# 16. B. Construct adjusted model with surrogate batch effects and observed batch effects
modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(edata))

# 16. C. Plotting SVA vs Combat adjusted values
par(mfrow=c(1,2))
plot(fitsv$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="SVA",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)

# 16. D. Plotting SVA vs Combat adjusted values
par(mfrow=c(1,2))
plot(fitsv$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="SVA",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
plot(fitsv$coefficients[2,], fit$coefficients[2,],col=2,
     xlab="SVA",ylab="linear model",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)


















