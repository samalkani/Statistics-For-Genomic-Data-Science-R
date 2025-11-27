# Dimension reduction for genomics

# 1. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 2. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 3. Making the circles on the plots filled solid
par(pch=19)

# 4. Installing Libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")
# BiocManager::install(c("devtools", "Biobase"))
# install.packages("usethis")

# 5. Load the Libraries
library(devtools)
library(Biobase)

# 6. Load the data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
ls()

# 7. Subtract out the row means less than 100 to aid visualization
edata = edata[rowMeans(edata) > 100, ]
dim(edata)

# 8. Apply the Log transform, to make the scale of the data easier to work with
edata = log2(edata + 1)

# 9. Centre the data set remove row and column means to expose variation between samples and genes
edata_centered = edata - rowMeans(edata)

# 10. Apply singular value decomposition to centred data
svd1 = svd(edata_centered)
names(svd1)
length(svd1$d)
svd1$d


# 11. Dimensions of the data set from where the 'd' matrix derives its 129 values
dim(edata)

# 12. Dimensions of the 'u' and 'v' matrices also has 129 values

# 12. A. U component or matrix (variation across samples)
dim(svd1$u)


# 12. B. V component or matrix (variation across genes)
dim(svd1$v)

# 13. Plot the 'd' values the singular values
plot(svd1$d,ylab="Singular value",col=2)

# 14. Plotting the Variance explained
plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col=2)

# 15. Plotting 1st principal component (Right Singular Vector) and the 2nd
par(mfrow=c(1,2))
plot(svd1$v[,1],col=2,ylab="1st PC")
plot(svd1$v[,2],col=2,ylab="2nd PC")

# 16. Plotting PC1 against PC2
par(mfrow=c(1,1))
plot(svd1$v[,1],svd1$v[,2],col=2,ylab="2nd PC", xlab="1st PC")

# 17. Plotting PC1 against PC2 colouring the PC's by the particular study
plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC", xlab="1st PC", col=as.numeric(pdata$study))

# 18. Box plot of PC1 against the two different study groups
boxplot(svd1$v[,1] ~ pdata$study, border=c(1,2))
points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)),col=as.numeric(pdata$study))

# 19. Plot of PC1 against SVD1 (not the same thing) even though the terms are used interchangeably
pc1 = prcomp(edata)
plot(pc1$rotation[,1],svd1$v[,1])

# 20. Perform a second centring subtracting column means then plot PC1 vs SVD1
edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)
plot(pc1$rotation[,1],svd2$v[,1],col=2)

# 21. Plotting SVD1 with and without outlier
edata_outlier = edata_centered
edata_outlier[1,] = edata_centered[1,] * 10000
svd3 = svd(edata_outlier)
plot(svd1$v[,1], svd3$v[,1], xlab="Without outlier", ylab="With outlier")

# 22. Highly correlated plot of SVD1 + outlier vs outlier itself
par(mfrow=c(1,1))
plot(svd3$v[,1],edata_outlier[1,],col=4)















