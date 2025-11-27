# Calculating Statistics in R

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing genefilter
# BiocManager::install(c("genefilter"), type = "source", force = TRUE)

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load Libraries
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)

# 5. Plotting Parameters
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
par(pch=19)

# 6. Load Data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
bot
pdata=pData(bot)
head(pdata)
edata=as.matrix(exprs(bot))
head(edata)
fdata = fData(bot)
head(fdata)
ls()

# 7. Log2 Transformation of data
edata = log2(as.matrix(edata) + 1)

# 8. Remove lowly expressed gene data (gene counts below 10)
edata = edata[rowMeans(edata) > 10, ]

# 9. Performing row t-tests to compare between groups (different strains)
tstats_obj = rowttests(edata,pdata$strain)
names(tstats_obj)

# 10. Histogram of T-statistic for every gene in the edata for the two strains
hist(tstats_obj$statistic,col=2)

# 11. Display the levels of the lanes variable
table(pdata$lane.number)

# 12. Performing row F-tests to compare between multiple groups (e.g lane.number)
fstats_obj = rowFtests(edata,as.factor(pdata$lane.number))
names(fstats_obj)

# 13. F-statistic for the lane.number variable for all genes in edata
hist(fstats_obj$statistic,col=2)

# 14. Fit a non-adjusted model of the strain variable using the limma package
# which gives moderated results (T-statistic is not inflated by small variance)
mod = model.matrix(~ pdata$strain)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)

# 15. Display Test statistics for the intercept and the strain variable
head(ebayes_limma$t)

# 16. Plot of moderated t-statistic from the limma package vs t-statistic (unmoderated)
plot(ebayes_limma$t[,2],-tstats_obj$statistic,col=4,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),col="darkgrey",lwd=3)

# 17. Use the Limma package to give an adjusted model and adjusted t-statistic
# This cannot be done with the Genefilter package
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_limma_adj = lmFit(edata,mod_adj)
ebayes_limma_adj = eBayes(fit_limma_adj)
head(ebayes_limma_adj$t)

# 18. Plot adjusted, moderated t-stat from Limma package vs unadjusted, unmoderated t-stat
plot(ebayes_limma_adj$t[,2],-tstats_obj$statistic,col=3,
     xlab="Moderated, Adjusted T-stat",ylab="T-stat")
abline(c(0,1),lwd=3,col="darkgrey")

# 19. Using Limma package for multiple level variable lane.number
mod_lane = model.matrix(~ as.factor(pdata$lane.number))
fit_limma_lane = lmFit(edata,mod_lane)
ebayes_limma_lane = eBayes(fit_limma_lane) 
head(ebayes_limma_lane$t)

# 20. Find the most significant level in the lane variable for each gene
top_lane = topTable(ebayes_limma_lane, coef=2:7,number=dim(edata)[1],sort.by="none")
head(top_lane)

# 21. Plot of the moderated version of the F-statistic we saw earlier
# Remember F-statistic is for multiple level variables
plot(top_lane$F,fstats_obj$statistic,
     xlab="Moderated F-statistic",ylab="F-statistic",col=3)

# 22. Using the Edge package for multi-level lane variable
edge_study = build_study(edata, grp = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
names(qval)

# 23. Plot of comparison of F-statistic from the gene filter and the edge packages
plot(qval$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter")

# 24. Use the Edge package to construct a new model with adjustment
edge_study2 = build_study(edata, grp = as.factor(pdata$lane.number),
                          adj.var=pdata$strain)
de_obj2 = lrt(edge_study2)
qval2 = qvalueObj(de_obj2)
plot(qval2$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter")










