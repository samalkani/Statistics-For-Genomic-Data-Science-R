# Exploratory Analysis in R Part 2

##### Part 1 - Basic Checks

# 1. Installing packages

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")

# install.packages("rlang")
# BiocManager::install(c("org.Hs.eg.db", "gplots"))
# devtools::install_github('alyssafrazee/RSkittleBrewer')
# BiocManager::install(c("AnnotationDbi", "Seqinfo"))

# BiocManager::valid()

# BiocManager::install(c(
#   "affxparser", "affyio", "AnnotationHub", "Biobase", "BiocParallel", "Biostrings", "DelayedArray", "DESeq2", "edgeR",
#   "genefilter", "GenomicAlignments", "GenomicRanges", "h5mread", "illuminaio", "IRanges", "limma", "multtest", "oligo",
#   "org.Hs.eg.db", "preprocessCore", "promises", "pwalign", "rhdf5", "rhdf5filters", "Rhdf5lib", "Rhtslib", "Rsamtools",
#   "rtracklayer", "S4Arrays", "S4Vectors", "ShortRead", "SparseArray", "sparseMatrixStats", "UCSC.utils", "VariantAnnotation",
#   "xfun", "XML", "XVector" 
# ), update = TRUE, ask = FALSE, force = TRUE)

# 2. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 3. Use the palette command to direct R to use those colours outlined above
palette(tropical)

# 4. Making the circles on the plots filled solid
par(pch=19)

# 5. Load Libraries
library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
# library(org.Hs.eg.db)
# library(DBI)
# library(AnnotationDbi)

# 6. Load in data from a connection
con=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)

# 7. Close the connection
close(con)

# 8. Display body map expression set
ls()

# 9. Re-assign bodymap.eset to a smaller name
bm = bodymap.eset

# 10. Extract phenotype data from expression set
pdata=pData(bm)

# 11. Extract expression data from the expression set
edata=exprs(bm)

# 12. Extract feature data from the expression set
fdata = fData(bm)

# 13. Display new variables
ls()

# 14. Displaying gender type information from phenotype data
table(pdata$gender)

# 15. Cross-Tabulation of Gender and Race
table(pdata$gender,pdata$race)

# 16. Summary of distribution for each column
summary(edata)

# 17. Finding missing values
table(pdata$age)

# 18. Use option useNA to include NA's in table
table(pdata$age,useNA="ifany")

# 19. Check for other common missing names
sum(pdata$age==" ")

# 20. Check for other common missing names
sum(pdata$age==" ", na.rm = TRUE)

# 21. Check genomic data for NAs
is.na(edata)[1,]

# 22. Check genomic data for NAs sum all values
sum(is.na(edata))

# 23. Make the distribution of NA's by genes
gene_na = rowSums(is.na(edata))
table(gene_na)
gene_na

# 24. Make the distribution of NA's by samples
sample_na = colSums(is.na(edata))
table(sample_na)
sample_na

# 25. Check to see if the dimensions of the expression set match up

# 25. A. Feature Data
dim(fdata)

# 25. B. Phenotype Data
dim(pdata)

# 25. C. Expression data
dim(edata)

##### Part 2 - Plotting the Data

# 26. A. Box plot
boxplot(edata[,1])

# 26. B. Log transform Expression data box plot
boxplot(log2(edata[,1]+1))

# 26. C. Display a box plot of the entire matrix not just one column
boxplot(log2(edata+1),col=2,range=0)

# 27. Two Sample Histogram plots of expression data side by side
par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2)

# 28. One sample Density plot of expression data
par(mfrow=c(1,1))
plot(density(log2(edata[,1]+1)),col=2)

# 29. Two sample Density plot of expression data (overlay plots)
plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=3)

# 30. Two sample Q-Q plot of expression data (Overlay plots)
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3)

# 31. Two sample Q-Q plot of expression data (Overlay plots)
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3)

abline(c(0,1))

# 32. MA or Bland-Altman plot
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)

# 33. Use filtering to keep only data points with a mean > 1

# 33. A. Install dplyr package
# install.packages("dplyr")
library(dplyr)

# 33. B. Convert matrix to dataframe to do filtering
edata = as.data.frame(edata)

# 33. C. Filter data frame to keep means > 1
filt_edata = filter(edata,rowMeans(edata)>1)

# 33. D. Dimensions of this new data set
dim(filt_edata)

# 33. E. Convert filtered dataframe back to a matrix to create box plot
boxplot(as.matrix(log2(filt_edata+1)),col=2)











