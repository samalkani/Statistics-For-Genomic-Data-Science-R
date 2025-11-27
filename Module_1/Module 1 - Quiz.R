# Module 1 - Quiz

# Question 1

# Reproducibility is defined informally as the ability to recompute data analytic 
# results conditional on an observed data set and knowledge of the statistical pipeline 
# used to calculate them Peng 2011, Science. Replicability of a study is the chance that 
# a new experiment targeting the same scientific question will produce a consistent result 
# Asendorpf 2013 European Journal of Personality.

# Susan asks Joe for his data shared according to the data sharing plan discussed in
# the lectures. Which of the following are reasons the study may be reproducible, 
# but not replicable?
  
# o	Joe doesn't make the raw data accessible so Susan can't re-run his code.

# o	The code for processing the data is written as a collection of R scripts and 
#   not R markdown files, so it is difficult to understand what each step does.

# o	The processed data alone don't give enough information to replicate the entire 
#   study.

# o	The identified effect can be reproduced from Joe's code and data, but may be 
#   due only to random variation and not appear in future studies.

# Question 2

# Put the following code chunk at the top of an R markdown document called 
# test.Rmd but set eval=TRUE

# 2. A. Set eval option to True
knitr::opts_chunk$set(cache=TRUE)

# 2. B. Plot Graph
 x = rnorm(10)
 plot(x,pch=19,col="dodgerblue")

# 2. C. Draw Table
y = rbinom(20,size=1,prob=0.5)
table(y)

# o	The table is random each time you knit the document, but the plot is always the
# same after you knit it the first time.

# o	The plot and table are random every time you knit the document, except for the 
# last time.

# o	The plot is random the first time you knit the document. It is identical to the 
# first time the second time you knit the document. After removing the folders 
# test_cache and test_files they generate new random versions.

# o	The plot and table are random the first time you knit the document. They are 
# identical the second time you knit the document. After removing the folders 
# test_cache and test_files they are still identical.

# Question 3

# 3. A. Create a summarizedExperiment object and load libraries with the following code
library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
data(sample.ExpressionSet, package = "Biobase")
se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)

# 3. B. How do you access the genomic data for this object?
genomic_table = assay(se)
head(genomic_table)

# 3. C. How do you access the phenotype table?
phenotype_table = colData(se)
phenotype_table

# 3. D. How do you access the feature data?
feature_data = rowData(se)
feature_data

# 3. E. What is the unique additional information provided by rowRanges(se)?
additional_info = rowRanges(se)
head(additional_info)

# Look up the help files for summarizedExperiment with the code?summarizedExperiment. 
# How do you access the genomic data for this object? How do you access the phenotype 
# table? How do you access the feature data? What is the unique additional information 
# provided by rowRanges(se)?
  
# o	Get the genomic table with assay(se), get the phenotype table with colData(se), 
# get the feature data with rowData(se).rowRanges(se) gives information on the genomic 
# location and structure of the measured features.

# o	Get the genomic table with assay(se), get the phenotype table with colData(se), 
# get the feature data with rowRanges(se). rowRanges(se) gives the range of possible 
# values for the expression data.

# o	Get the genomic table with assay(se), get the phenotype table with pData(se), 
# get the feature data with rowData(se).rowRanges(se) gives information on the 
# genomic location and structure of the measured features.

# o	Get the genomic table with assay(se), get the phenotype table with colData(se),
# get the feature data with rowData(se).rowRanges(se) gives the range of possible 
# values for the expression data.

# Question 4

# Suppose that you have measured ChIP-Seq data from 10 healthy individuals and 
# 10 metastatic cancer patients. For each individual you split the sample into 
# two identical sub-samples and perform the ChIP-Seq experiment on each sub-sample. 
# How can you measure (a) biological variability, (b) technical variability and 
# (c) phenotype variability?
  
# o	(a) By looking at variation across samples from 10 different healthy individuals 

#   (b) By looking at variability between the measurements on the two sub-samples 
#       from the same sample and 

#   (c) by comparing the average measurements on the healthy individuals to the 
#       measurements on the individuals with cancer. 


# o	(a) By looking at variation across samples from 10 different individuals with cancer 

#   (b) By comparing the average variability in the cancer and normal individuals 

#   (c) By comparing the average measurements on the healthy individuals to the 
#       measurements on the individuals with cancer.

# o	(a) By looking at variation across replicate sub-samples within the normal individuals 

#   (b) By looking at variation across samples from 10 different healthy individuals 

#   (c) By comparing the average measurements on the healthy individuals to the 
#       measurements on the individuals with cancer.

# o	(a) & (b) By looking at variation across samples from 10 different healthy 
#       individuals. 

#   (c) by comparing the average measurements on the healthy individuals to the 
#       measurements on the individuals with cancer.

# Question 5

# 5. A. Load the Bottomly and the Bodymap data sets with the following code:
  
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
pdata_bot
  
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)
pdata_bm
  
#   Just considering the phenotype data what are some reasons that the Bottomly 
#   data set is likely a better experimental design than the Bodymap data? Imagine 
#   the question of interest in the Bottomly data is to compare strains and in the 
#   Bodymap data it is to compare tissues.
  
#   o	The Bottomly data has biological replicates for each group but the Bodymap 
#   data does not. 
  
#   o	The Bodymap data has measured more levels of the outcome of interest (tissues)
#   than the Bottomly data has measured (strains).
  
#   o	The Bottomly data set does not measure the age of the mice.
  
#   o	Most of the tissues in the Bodymap data have a consistent number of technical 
#   replicates (2).
  
#   The covariates in the Bottomly data set (experiment number, lane number) are 
#   balanced with respect to strain.
#   The covariates in the Bodymap data set (gender, age, number of technical replicates)
#   are not balanced with
#   respect to tissue.
  
# Question 6

# What are some reasons why this plot is not useful for comparing the number of 
# technical replicates by tissue (you may need to install the plotrix package).

# 6. A. Load Biobase library
library(Biobase)

# 6. B. Load Body Map data set
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

# 6. C. Plot the 3D Pie Chart
# install.packages('plotrix')
library(plotrix)
pie3D(pdata_bm$num.tech.reps,labels=pdata_bm$tissue.type)

# o	The plot is in 3-d so it makes it hard to compare the angles.

# o	There are a large number of data points underlying each wedge and you can't 
#   see them. 

# o	The plot would be much easier to see if the pie chart were rotated by 90 
#   degrees from its current position. 

# o	There is nothing wrong with the plot, it accurately shows how many replicates 
#   of each type there are. 

# The b

# Question 7

### 7. A.	Load the Bottomly data:

library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset
bm
edata = exprs(bm)
head(edata)

# Which of the following code chunks will make a heatmap of the 500 most highly 
# expressed genes (as defined by total count), without re-ordering due to clustering? 
# Are the highly expressed samples next to each other in sample order?
  
# 7. B. Code Chunk 1
  
row_sums = rowSums(edata)
edata = edata[order(-row_sums),]
index = 1:500
heatmap(edata[index,],Rowv=NA,Colv=NA)

# Yes they are

# 7. C. Code Chunk 2

# row_sums = rowSums(edata)
# index = which(rank(row_sums) < 500 )
# heatmap(edata[index,],Colv=NA)



# The highly expressed samples are not next to each other.

# 7. D. Code Chunk 3

# row_sums = rowSums(edata)
# edata = edata[order(row_sums),]
# index = which(rank(-row_sums) < 500 )
# heatmap(edata[index,],Rowv=NA,Colv=NA)

# No they are not next to each other.

# 7. E. Code Chunk 4

# row_sums = rowSums(edata)
# index = which(rank(-row_sums) < 500 )
# heatmap(edata[index,],Rowv=NA,Colv=NA)


# Question 8

# Load the Bodymap data using the following code:
  
# 8. A. Loading Body map data

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset
bm
pdata = pData(bm)
head(pdata)
edata = exprs(bm)
head(edata)

# 8. B. MA or Bland-Altman plot
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)

## 8. C. MA or Bland-Altman plot (using DESeq2 package)
library(DESeq2)
edata_log = rlog(edata + 1)
mm = edata_log[,1] - edata_log[,2]
aa = edata_log[,1] + edata_log[,2]
plot(aa,mm,col=2)

# Make an MA-plot of the first sample versus the second sample using the log2 
# transform (hint: you may have to add 1 first) and the rlog transform from the 
# DESeq2 package. How are the two MA-plots different? Which kind of genes appear 
# most different in each plot?
  
# The plots are very different as the log2 plot seems to shrink low abundance 
# genes more and the rlog plot seems to shrink high abundance genes more. The 
# genes in the middle of the distribution show the biggest differences.

# The plots look pretty similar, but there are two strong diagonal stripes** 
# (corresponding to the zero count genes) in the log2 plot. In both cases,** 
# the genes in the middle of the expression distribution show the biggest** 
# differences, but the low abundance genes seem to show smaller differences** 
# with the rlog transform.
  
# The plots look pretty similar, but there are two strong diagonal stripes 
# (corresponding to the zero count genes) in the rlog plot. In both cases, 
# the genes in the middle of the expression distribution show the biggest differences, 
# but the low abundance genes seem to show smaller differences with the log2 transform.

# The plots are nearly identical. Both transforms seem to deal with the low 
# abundance genes, including the zero genes the same way. The high-abundance 
# genes show the most differences.

# Question 9

# 9. A.	Load the Montgomery and Pickrell eSet

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
library(dendextend)

# Cluster the data in three ways:
  
# With no changes to the data

# After filtering all genes with rowMeans less than 100

# After taking the log2 transform of the data without filtering

# Color the samples by which study they came from (Hint: consider using the 
# function myplclust.R in the package rafalib available from CRAN and looking at
# the argument lab.col.) How do the methods compare in terms of how well they 
# cluster the data by study? Why do you think that is?
  
# 9. B. With no changes to the data
dist1 = dist(t(edata))
hclust1 = hclust(dist1)
par(mar=c(0, 4, 4, 2))
plot(hclust1, hang = -1, main="origin", labels=FALSE)

# 9. C. After filtering all genes with rowMeans less than 100
low_genes = rowMeans(edata) < 100
filter_edata = filter(edata, !low_genes)
f_dist1 = dist(t(filter_edata))
f_hclust1 = hclust(f_dist1)
par(mar=c(0, 4, 4, 2))
plot(f_hclust1, hang = -1, main="remove low expression", labels=FALSE)


# 9. D. After taking the log2 transform of the data without filtering
log_edata = log2(edata + 1)
l_dist1 = dist(t(log_edata))
l_hclust1 = hclust(l_dist1)
par(mar=c(0, 4, 4, 2))
plot(l_hclust1, hang=-1, main="perform log2 transform", labels=FALSE)


# Clustering with or without filtering is about the same. Clustering after the 
# log2 transform shows better clustering with respect to the study variable. The 
# likely reason is that the highly skewed distribution doesn't match the Euclidean 
# distance metric being used in the clustering example.

# Clustering is identical with all three approaches and they show equal clustering. 
# The distance is an average over all the dimensions so it doesn't change. 

# Clustering with or without log2 transform is about the same. Clustering after 
# filtering shows better clustering with respect to the study variable. The reason 
# is that it is just the lowly expressed genes that make the distance hard to calculate. 

# Clustering with or without log2 transform is about the same. Clustering after** 
# filtering shows better clustering with respect to the study variable. The reason** 
# is that the lowly expressed genes have some extreme outliers that skew the calculation.**
  
# Question 10
  
# 10. A. Loading Montpick data set
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# Cluster the samples using k-means clustering after applying the log2 transform (be sure to add 1). Set a seed for reproducible results (use set.seed(1235)). If you choose two clusters, do you get the same two clusters as you get if you use the cutree function to cluster the samples into two groups? Which cluster matches most closely to the study labels?
  
# 10. B. Log Transform
edata = log2(edata + 1)

# 10. C. Perform k-means clustering - two clusters
set.seed(1235)
k2 = kmeans(edata,centers=2)
matplot(t(k2$centers),col=1:2,type="l",lwd=3)

# 10. D. Hierarchicial Clustering
dist1 = dist(t(edata))
hclust1 = hclust(dist1)
tree = cutree(hclust1, 2)
par(mar=c(0, 4, 4, 2))
plot(hclust1, tree, main="cutree")


# They produce the same answers and match the study variable equally well. 

# They produce different clustering’s with hierarchical clustering more closely 
# matching the study variable. K-means clustering is too random to pick up the 
# study difference. 

# They produce different answers, with k-means clustering giving a much more 
# unbalanced clustering. The hierarchical clustering matches study better. 

# They produce different answers, with hierarchical clustering giving a much** 
# more unbalanced clustering. The k-means clustering matches study better.**



