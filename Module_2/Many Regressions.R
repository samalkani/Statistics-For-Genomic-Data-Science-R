# Many Regressions

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing edge
# BiocManager::install(c("edge"), type = "source", force = TRUE)

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load Libraries
library(devtools)
library(Biobase)
library(limma)
library(edge)

# 5. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 6. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 7. Making the circles on the plots filled solid
par(pch=19)

# 8. Load Data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
head(bot)
pdata=pData(bot)
head(pdata)
edata=as.matrix(exprs(bot))
head(edata)
fdata = fData(bot)
head(fdata)
ls()

# 9. Remove the lowly expressed genes
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

# 10. Dimensions of new matrix
dim(edata)

# 11. Construct a model matrix
mod = model.matrix(~ pdata$strain)
mod

# 12. Fit a Model using the lm() function

# 12. A. Fit the model using the lm() function - Model 1
fit = lm.fit(mod,t(edata))
names(fit)

# 12. B. Coefficients from the first model fit - Model 1
fit$coefficients[,1]

# 12. C. Coefficients matrix is very big (first 10 models) - Model 1
fit$coefficients[,1:10]

# 13. Compare first model in the model matrix to fitting a model on individual data sets

# 13. A. 1st model in the model matrix (Fits many models simultaneously) - Model 1
fit$coefficients[,1]

# 13. B. Fitting model to the 1st data set (Fits one model at a time - Model 1
library(broom)
tidy(lm(as.numeric(edata[1, ]) ~ pdata$strain))

# 13. C. Display distributions of the many regression models - Model 1
par(mfrow=c(1,2))
hist(fit$coefficients[1,],breaks=100,col=2,xlab="Intercept")
hist(fit$coefficients[2,],breaks=100,col=2,xlab="Strain")
abline(v=0,lwd=3,col=1)

# 13. D. Plot residuals from the 1st and 2nd data set - Model 1
par(mfrow=c(1,2))
plot(fit$residuals[,1],col=2)
plot(fit$residuals[,2],col=2)

# 13. E. Fitting a model with adjustments - Model 1
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_adj = lm.fit(mod_adj,t(edata))
fit_adj$coefficients[,1]

# 14. Using the Limma package to fit models

# 14. A. Fitting multiple models using the Limma package – Model 2
fit_limma = lmFit(edata,mod_adj)
names(fit_limma)

# 14. B. Display Coefficients – Model 2
fit_limma$coefficients[1,]

# 14. C. Coefficients for the lm.fit() package appear in the 1st column – Model 1
fit$coefficients[,1]

# 15. Using the Edge Package to fit Models

# 15. A. Using the Edge package to fit models (no knowledge of model matrices required) – Model 3
edge_study = build_study(data=edata,grp=pdata$strain,adj.var=as.factor(pdata$lane.number))
fit_edge = fit_models(edge_study)
summary(fit_edge)

# 15. B. Display coefficients for the 1st gene (using the Edge package) – Model 3
fit_edge@beta.coef[1,]

# 15. C. Display coefficients for the 1st gene (using the Limma Package) - Model 2
fit_limma$coefficients[1,]





