# Permutation in R

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing genefilter
# BiocManager::install(c("genefilter"), type = "source", force = TRUE)

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load the Libraries
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)

# 5. Plotting Parameters

# 5. A. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 5. B. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 5. C. Making the circles on the plots filled solid
par(pch=19)

# 6. Load the Data
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

# 7. Log2 transform data and filter out the lowly expressed genes
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

# 8. Use the rowT-test on the strain variable
tstats_obj = rowttests(edata,pdata$strain)

# 9. Display observed statistics through histogram plot
hist(tstats_obj$statistic,col=2,xlim=c(-5,2))

# 10. Setting up the permutation
set.seed(135)                     # Random seed
strain = pdata$strain             # Original order
strain0 = sample(strain)          # Random order

# 11. Display original order
strain

# 12. Display random order
strain0

# 13. Re-perform row T-tests on strain0 variable (with randomly permuted labels)
tstats_obj0 = rowttests(edata,strain0)
hist(tstats_obj0$statistic,col=2,xlim=c(-5,2))

# 14. Display Quantiles of the permuted statistics
quantile(tstats_obj0$statistic)

# 15. Display quantiles of the observed statistics
quantile(tstats_obj$statistic)



















