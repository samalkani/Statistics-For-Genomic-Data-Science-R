# Quantile Normalisation in R

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing preprocessCore
BiocManager::install(c("preprocessCore"), force = TRUE)

# 3. Load Libraries
library(devtools)
library(Biobase)
library(preprocessCore)

# 4. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 5. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 6. Making the circles on the plots filled solid
par(pch=19)

# 7. Load Data and extract out the phenotype, genotype & feature data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
ls()

# 8. Log Transform and remove outlying values
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 3, ]

# 9. Dimensions of expression data
dim(edata)

# 10. Plot density for the first sample
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))

# 11. Plot for the remaining 19 samples (overlay plot)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}


# 12. Quantile Normalization & dimensions of normalized data (matrix)
norm_edata = normalize.quantiles(as.matrix(edata))
dim(norm_edata)

# 13. Plot density of Quantile normalized data
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.20))
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}

# 14. Plot of Normalized data for the first gene against study group
plot(norm_edata[1,],col=as.numeric(pdata$study))

# 15. Plot of SVD1 vs SVD2 by study group
svd1 = svd(norm_edata - rowMeans(norm_edata))
plot(svd1$v[,1],svd1$v[,2],xlab="SVD1",ylab="SVD2", col=as.numeric(pdata$study))













