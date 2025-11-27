# Data Transforms

# 1. Installing packages

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")


# BiocManager::install(c("devtools", "Biobase"))

# 2. Load Libraries
library(devtools)
library(Biobase)

# 3. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 4. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 5. Making the circles on the plots filled solid
par(pch=19)

# 6. Load Data from URL
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)

# 7. Body Map eset
bm = bodymap.eset

# 8. Phenotype Data
pdata=pData(bm)

# 9. Gene Expression Data
edata=as.data.frame(exprs(bm))

# 10. Feature Data
fdata = fData(bm)

# 11. List
ls()

# 12. Histogram of Normal data (simulated)
hist(rnorm(1000),col=2)

# 13. Realistically most genomic data is skewed (sample 1 from expression set)
hist(edata[,1],col=2,breaks=100)

# 14. Realistically most genomic data is skewed (sample 1 from exp set, 1000 bins)
hist(edata[,1],col=2,breaks=1000)

# 15. Log transform of skewed data for sample 1
hist(log(edata[,1]),col=2,breaks=100)

# 16. Log transform genomic data with mainly zeros gives -Inf
min(log(edata))
quantile(log(edata[,1]))

# 17. Add 1 to the zero counts so log transform is properly defined
hist(log(edata[,1] + 1),col=2)

# 18. Log2 transform of sample 1 genomic data
hist(log2(edata[,1] + 1),breaks=100,col=2)

# 19. Log2 transform of sample 1 genomic data with limits on X + Y axis
hist(log2(edata[,1] + 1),breaks=100,col=2,xlim=c(1,15),ylim=c(0,400))

# 20. Counting the number of row which have zero
hist(rowSums(edata==0),col=2)

# 21. Identifying the low genes
low_genes = rowMeans(edata) < 5
table(low_genes)

# 22. Filter out low gene values, counts < 5
library(dplyr)
filt_edata = filter(edata,!low_genes)
dim(filt_edata)

# 23. Summary of edata
summary(edata)

# 24. Identifying the low genes based on row Medians instead of row means
low_genes2 = rowMedians(as.matrix(edata)) < 5
table(low_genes2,low_genes)

# 25. Comparing filtering by mean value compared with median value
filt_edata2 = filter(edata,!low_genes2)

# 26. Dimensions filtering by median
dim(filt_edata2)

# 27. Dimensions filtering by mean
dim(filt_edata)

# 28. Plot with filtered data by mean
hist(log2(filt_edata[,1] + 1),col=2)


















