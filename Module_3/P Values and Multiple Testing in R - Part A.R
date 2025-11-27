# P values and Multiple Testing in R: Part A

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing qvalue
# BiocManager::install(c("qvalue"), type = "source", force = TRUE)

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load Libraries
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)
library(qvalue)

# 5. Plotting Parameters
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
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

# 7. Log2 transform and removing lowly expressed genes
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

# 8. Calculating the F-statistic and p-values using the gene filter package
fstats_obj = rowFtests(edata,as.factor(pdata$strain))
hist(fstats_obj$p.value,col=2)

# 9. Construct an adjusted model using the edge package, Strain (variable of interest),
# Lane (adjustment variable), model is not moderated (may get inflated statistics)
edge_study = build_study(edata, grp = pdata$strain, adj.var = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
hist(qval$pvalues,col=3)

# 10. P-values for moderated statistics with Limma package
mod = model.matrix(~ pdata$strain + pdata$lane.number)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1])$P.Value
hist(limma_pvals,col=4)

# 11. Calculating empirical permutation p-values with edge

# 11. A. Setting randomisation seed
set.seed(3333)

# 11. B. Number of permutations
B = 1000

# 11. C. Row T-tests on observed data
tstats_obj = rowttests(edata,pdata$strain)

# 11. D. Empty Matrix ready to be populated by simulated data (permutations)
tstat0 = matrix(NA,nrow=dim(edata)[1],ncol=B)

# 11. E. Re-assign T-statistic from observed data
tstat = tstats_obj$statistic

# 11. F. Re-assign phenotype data for the strain variable
strain = pdata$strain

# 11. G. Populate empty matrix with simulated data (permutations)
for(i in 1:B){
  strain0 = sample(strain)
  tstat0[,i] = rowttests(edata,strain0)$statistic
}
# 11. H. Generate P-values from t-statistics from observed and simulated data
emp_pvals = empPvals(tstat,tstat0)
hist(emp_pvals,col=2)






