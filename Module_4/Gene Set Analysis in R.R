# Gene Set Analysis in R

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "goseq"
# BiocManager::install("goseq")

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load Libraries
library(goseq)

# 5. Display data
head(supportedGenomes())

# 6. Display data – Available cases
head(supportedGeneIDs())

# 7. Load data from goseq package
temp_data =read.table(system.file("extdata","Li_sum.txt",
                                  package="goseq"),sep="\t",
                      header=TRUE,
                      stringsAsFactors=FALSE)

# 8. First row of the data set
temp_data[1,]

# 9. Remove 1st column of data set and replace it with rownames = gene names
expr= temp_data[,-1]
rownames(expr) = temp_data[,1]

# 10. Remove lowly expressed genes
expr = expr[rowMeans(expr) > 5,]

# 11. Create a group variable, the data doesn't have a group variable
grp=factor(rep(c("Control","Treated"),times=c(4,3)))
grp

# 12. Create a pdata (phenotype) object with just the grp variable in it
pdata  = data.frame(grp)
pdata[1,]

# 13. Passing expression(expr), phenotype data(pdata) and the grp variable you want to model
library(DESeq2)
de = DESeqDataSetFromMatrix(expr, pdata, ~grp)

# 14. Use DESeq to identify differentially expressed genes
de_fit = DESeq(de)

# 15. Calculate the results using the results function
de_results = results(de_fit)

# 16. Display differentially expressed results
head(de_results)

# 17. Remove the genes where we can’t calculate a statistic

# 17. A. Filter genes that are differentially expressed < 0.05 (FDR)
genes = as.integer(de_results$padj < 0.05)

# 17. B. Remove data were there wasn't enough data
not_na = !is.na(genes)

# 17. C. Apply row names to match genes to genomic data
names(genes) = rownames(expr)

# 17. D. Remove genes we can't calculate a statistic
genes = genes[not_na]

# 18. We need to set up a weighting function for all the genes in that genome
pwf=nullp(genes,"hg19","ensGene")

# 19. Use goseq package to calculate the statistics for enrichment
GO.wall=goseq(pwf,"hg19","ensGene")

# 20. Display results
head(GO.wall)

# 21. Selecting a specific category of interest for your analysis
GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
head(GO.MF)











