# Generalized linear models (GLMs) in R

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "MASS", "DESeq2"
# BiocManager::install(c("MASS", "DESeq2"), type = "source", force = TRUE)

# 3. Troubleshooting installing packages
BiocManager::valid()

# 4. Load Libraries
library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)

# 5. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 6. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 7. Making the circles on the plots filled solid
par(pch=19)

# 8. Load Data from snpStats package
data(for.exercise)
ls()

# 9. Subset data so the data is computationally easier to handle
use <- seq(1, ncol(snps.10), 10)
head(use)

# 10. Take the created subset to reduce the size of the data set
sub.10 <- snps.10[,use]
sub.10

# 11. Calculate the Principal Components (PC's)

# 11. A. Right transformation of the data to get the xx transpose
xxmat <- xxt(sub.10, correct.for.missing=FALSE)

# 11. B. Calculate the Eigen vectors 
evv <- eigen(xxmat, symmetric=TRUE)

# 11. C. Calculate the PC's
pcs <- evv$vectors[,1:5]
head(pcs)

# 11. D.Dimensions of PC's matrix
dim(pcs)

# 12. Subset out data into R format so it’s easy to process
snpdata = sub.10@.Data
dim(snpdata)

# 13. Case Control status
status = subject.support$cc
status

# 14. Take the 1st SNP and extract the data for the 1st SNP
snp1 = as.numeric(snpdata[,1])
length(snp1)
length(status)

# 15. Display the way the SNP data is coded
table(snp1)
snp1

# 16. Recode the zero values (missing values) as NA so R will know what to do with them
snp1[snp1==0] = NA

# 17. Fit a generalized linear model (GLM), Logistic Regression
glm1 = glm(status ~ snp1,family="binomial")
tidy(glm1)

# 18. Fitting the dominant model (AA) where only Aa or aa carries risk using a dummy variable
snp1_dom = (snp1 == 1)
glm1_dom = glm(status ~ snp1_dom,family="binomial")
tidy(glm1_dom)

# 19. Previous model using SNP variable for comparison with the dominant model
glm1 = glm(status ~ snp1,family="binomial")
tidy(glm1)

# 20. Adjust for population structure using principal components in the GLM model
glm2 = glm(status ~ snp1 + pcs[,1:5],family="binomial")
tidy(glm2)

# 21. Construct many GLM models, un-adjusted model
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)
slotNames(glm_all)

# 22. Chi-Squared statistics + plot
qq.chisq(chi.squared(glm_all),df=1)

# 23. Construct many GLMs, adjusted models with PC's added
glm_all_adj = snp.rhs.tests(status ~ pcs,snp.data=sub.10)
qq.chisq(chi.squared(glm_all_adj),df=1)

# 24. Load the Bottomly data to construct a Poisson model from count data
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

# 25. Filter out the lowly expressed genes
edata = edata[rowMeans(edata) > 10, ]

# 26. Fit a GLM, Poisson model on the Bottomly data
glm3 = glm(edata[1, ] ~ pdata$strain,family="poisson")
tidy(glm3)

# 27. Fit a GLM, Negative Binomial model, that offers more flexibility in terms 
# of variance than the Poisson model (where mean & variance are the same)
glm.nb1 = glm.nb(edata[1, ] ~ pdata$strain)
tidy(glm.nb1)

# 28. Fit a GLM, Negative Binomial model on many genes using the DESeq package

# 28. A. Build a DESeq data set
de = DESeqDataSetFromMatrix(edata, pdata, ~strain)

# 28. B. Construct the models
glm_all_nb = DESeq(de)

# 28. C. Accessing the results of the Negative Binomial Model on many genes
result_nb = results(glm_all_nb)

# 28. D. Histogram of the statistics
hist(result_nb$stat)








