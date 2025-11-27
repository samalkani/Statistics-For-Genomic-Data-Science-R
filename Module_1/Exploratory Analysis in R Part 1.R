# Exploratory Analysis in R Part 1

# 1. Installing packages

 if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install(version = "3.22")
 
 BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
 
 BiocManager::available("Seqinfo")
 
 if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 pkgs <- rownames(installed.packages())
 BiocManager::install(pkgs, type = "source", checkBuilt = TRUE)
 
# install.packages(c("devtools","gplots"))
# devtools::install_github('alyssafrazee/RSkittleBrewer')

# Cannot install because cannot install Seqinfo package
BiocManager::install(c("org.Hs.eg.db","AnnotationDbi"))

# Use the latest version of Bioconductor for your version of R
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# Troubleshooting with Biocoductor packages
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
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")

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



































