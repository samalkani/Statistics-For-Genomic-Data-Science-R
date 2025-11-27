# Module 2 - Quiz

# Question 1

# 
library(devtools)
library(Biobase)

# 1. B. Load the Montgomery and Pickrell eSet:
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
ls()
pdata=pData(mp)
head(pdata)
edata=as.data.frame(exprs(mp))
head(edata)
fdata = fData(mp)
head(fdata)
  
# What percentage of variation is explained by the 1st principal component in the data set if you:
    
#   a)	Do no transformations?
  
#   b)	log2(data + 1) transform?
  
#   c)	log2(data + 1) transform and subtract row means?
  
# 1. C. No transformations
svd1 = svd(edata)
ori_pca = svd1$d^2/sum(svd1$d^2)
ori_pca[1]

# 1. D. log2 transform
edata_log2 = log2(edata + 1)
svd2 = svd(edata_log2)
log2_pca = svd2$d^2/sum(svd2$d^2)
log2_pca[1]

# 1. E. log2 transform, subtract row means
edata_centered = edata_log2 - rowMeans(edata_log2)
svd3 = svd(edata_centered)
centered_data_pca = svd3$d^2/sum(svd3$d^2)
centered_data_pca[1]

# o	a. 0.97 b. 0.97 c. 0.35

# o	a. 0.97 b. 0.97 c. 0.97

# o	a. 0.97 b. 0.89 c. 0.35

# o	a. 0.89 b. 0.97 c. 0.35

# Question 2

# 2. A. Load the Montgomery and Pickrell eSet:

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# Perform the log2(data + 1) transform and subtract row means from the samples. 
# Set the seed to 333 and use k-means to cluster the samples into two clusters. 
# Use svd to calculate the singular vectors. What is the correlation between the 
# first singular vector and the sample clustering indicator?

# 2. B. Log2 Transform and centre data by subtracting row means
edata_log2 = log2(edata + 1)
edata_centered = edata_log2 - rowMeans(edata_log2)

# 2. C. Use SVD to calculate the singular vectors
set.seed(333)
svd1 = svd(edata_centered)

# 2. D. Perform K-means Clustering and correlate 1st SVD and K-means cluster
edata_kmeans = kmeans(t(edata_centered), centers=2)
head(svd1$v[,1])
head(edata_kmeans$cluster)
tail(edata_kmeans$cluster)
cor.test(svd1$v[,1], edata_kmeans$cluster)

# o	-0.08

# o	0.84

# o	-0.52

# o	0.87

# Question 3

# 3. A. Load the Bodymap data with the following command

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset
edata = exprs(bm)
head(edata)
pdata_bm=pData(bm)
head(pdata_bm)

# Fit a linear model relating the first gene’s counts to the number of technical 
# replicates, treating the number of replicates as a factor. Plot the data for this 
# gene versus the covariate. Can you think of why this model might not fit well?

# 3. B. Fit linear model
lm1 = lm(edata[1,] ~ pdata_bm$num.tech.reps)
summary(lm1)

# 3. C. Plot the data
plot(pdata_bm$num.tech.reps,edata[1,])
abline(lm1, col=2, lwd=3)

# o	The difference between 2 and 5 technical replicates is not the same as the 
# difference between 5 and 6 technical replicates. 

# o	There are very few samples with more than 2 replicates so the estimates for 
# those values will not be very good. 

# o	The variable num.tech.reps is not a good measure of the technical variability 
# in the data.

# o	There may be different numbers of counts for different numbers of technical replicates.

# Question 4

# 4. A.	Load the Bodymap data with the following command

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset
head(bm)
edata = exprs(bm)
head(edata)
pdata_bm=pData(bm)
head(pdata_bm)

# Fit a linear model relating he first gene’s counts to the age of the person and 
# the sex of the samples. What is the value and interpretation of the coefficient 
# for age?
  
# 4. B. Fit linear model
lm2 = lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
summary(lm2)

# o	-23.25. This coefficient means that there is an average decrease of 23.91 in 
# the count variable per year within each gender. 

# o	2187.91. This means that for a person that is zero years old, the average count 
# is 2187.91

# o	-23.91. This coefficient means that there is an average decrease of 23.91 in 
# the count variable per year within each gender. 

# o	-207.26. This coefficient means that for each additional year of age, the 
# count goes down by an average of 207.26 for a fixed sex.

# Question 5

# 5. A.	Load the Montgomery and Pickrell eSet:
  
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
ls()
mp = montpick.eset
head(mp)
pdata=pData(mp)
head(pdata)
edata=as.data.frame(exprs(mp))
head(edata)
fdata = fData(mp)
head(fdata)

# Perform the log2(data + 1) transform. Then fit a regression model to each sample
# using population as the outcome. Do this using the lm.fit function (hint: don't 
# forget the intercept). What is the dimension of the residual matrix, the effects 
# matrix and the coefficients matrix?

# 5. B. Log2 transform the expression data
edata = log2(edata + 1)

# 5. C. fit a regression model to each sample, using population as the outcome
mod = model.matrix(~ pdata$population)
head(mod)
fit = lm.fit(mod, t(edata))

# 5. D. dimension of the residual matrix
dim(fit$residuals)

# 5. E. Dimension of the effects matrix
dim(fit$effects)

# 5. F. Dimension of the coefficients matrix
dim(fit$coefficients)

# *	Residual matrix: 52580 x 129   
# Effects matrix: 52580 x129  
# Coefficients matrix: 2 x 52580 

# *	Residual matrix: 52580 x 129   
# Effects matrix: 129 x 52580
# Coefficients matrix: 2 x 52580 

# *	Residual matrix:  129 x 52580 
# Effects matrix: 129 x 52580
# Coefficients matrix: 1 x 52580 

# *	Residual matrix:  129 x 52580 
# Effects matrix: 129 x 52580
# Coefficients matrix: 2 x 52580

# Question 6

# 6. A.	Load the Montgomery and Pickrell eSet:

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
ls()
mp = montpick.eset
head(mp)
pdata=pData(mp)
head(pdata)
edata=as.data.frame(exprs(mp))
head(edata)
fdata = fData(mp)
head(fdata)

# Perform the log2(data + 1) transform. Then fit a regression model to each sample 
# using population as the outcome. Do this using the lm.fit function (hint: don't 
# forget the intercept). What is the effects matrix?

# 6. B. Log2 transform the expression data
edata = log2(edata + 1)

# 6. C. Fit a regression model to each sample, using population as the outcome
mod = model.matrix(~ pdata$population)
fit = lm.fit(mod, t(edata))

# o	The estimated fitted values for all samples for each gene, with the values for
# each gene stored in the columns of the matrix. 

# o	The model coefficients for all samples for each gene, with the values for each
# gene stored in the rows of the matrix. 

# o	The model coefficients for all samples for each gene, with the values for each
# gene stored in the columns of the matrix. 

# o	The estimated fitted values for all samples for each gene, with the values for 
# each gene stored in the rows of the matrix.

# Question 7

# 7. A.	Load the Bodymap data with the following command

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset
head(bm)
edata = exprs(bm)
head(edata)
pdata_bm=pData(bm)
head(pdata_bm)

# Fit many regression models to the expression data where age is the outcome variable 
# using the lmFit function from the limma package (hint: you may have to subset the 
# expression data to the samples without missing values of age to get the model to fit). 
# What is the coefficient for age for the 1,000th gene? Make a plot of the data and 
# fitted values for this gene. Does the model fit well?
  
# 7. B. Load Libraries
library(devtools)
library(Biobase)
library(limma)
library(edge)

# 7. C. Subset the expression data to the samples without missing values of age
pdata_bm = na.omit(pdata_bm)
edata = edata[,rownames(pdata_bm), drop=FALSE]
head(edata)

# 7. D. Fit many regression models to the expression data where age is the outcome
mod_adj = model.matrix(~ pdata_bm$age)
mod_adj
fit_limma = lmFit(edata,mod_adj)
dim(fit_limma$coefficients)
fit_limma$coefficients[1000,]

# 7. E. make a plot of the 1,000th gene and fitted values
intercept = fit_limma$coefficients[1000,][1]
slope = fit_limma$coefficients[1000,][2]
x = edata[1000,]*slope+intercept
plot(x,pdata_bm$age)

# o	-27.61. The model fits well since there seems to be a flat trend in the counts. 

# o	-27.61. The model doesn't fit well since there appears to be a non-linear trend
# in the data.

# o	2469.87. The model doesn't fit well since there appears to be a non-linear trend 
# in the data. 

# o	-27.61. The model doesn't fit well since there are two large outlying values 
# and the rest of the values are near zero.

# Question 8

# 8. A. Load the Bodymap data with the following command

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset
head(bm)
edata = exprs(bm)
head(edata)
pdata_bm=pData(bm)
head(pdata_bm)

# Fit many regression models to the expression data where age is the outcome variable 
# and tissue.type is an adjustment variable using the lmFit function from the limma 
# package (hint: you may have to subset the expression data to the samples without 
# missing values of age to get the model to fit). What is wrong with this model?
  
# 8. B. Subset the expression data to the samples without missing values of age
pdata_bm = na.omit(pdata_bm)
edata = edata[,rownames(pdata_bm), drop=FALSE]

# 8. C. Fit many regression models to the expression data where age is the outcome 
# and tissue type is the adjustment variable
mod_adj = model.matrix(~ pdata_bm$age + as.factor(pdata_bm$tissue.type))
mod_adj
fit_limma = lmFit(edata,mod_adj)

# 8. D. There are 18 levels in the tissue type variable
pdata_bm$tissue.type

# 8. E. There are only 16 data points per gene in the expression data (edata)
dim(edata)


# o	This model doesn't fit well since tissue.type makes more sense as an outcome 
# variable.

# o	Since tissue.type is a factor variable with many levels, this model has more 
# coefficients to estimate per gene (18) than data points per gene (16).

# o	The model doesn't fit well since age should be treated as a factor variable.

# o	The model doesn't fit well because there are a large number of outliers for 
# the white blood cell tissue. 

# Question 9

# 9. A. Multicollinearity between the population and study variables
head(pdata)
pop = as.numeric(as.factor(pdata$population))
pop
study = as.numeric(as.factor(pdata$study))
study
cor.test(pop, study)

# Why is it difficult to distinguish the study effect from the population effect 
# in the Montgomery Pickrell dataset from ReCount?
  
# o	The effects are difficult to distinguish because the study variable and 
# population variable are perfectly correlated. 

# o	The study effects and population effects are difficult to distinguish because 
# of confounding, we can estimate the population effects by adjusting for study. 

# o	The study effects and population effects are difficult to distinguish because 
# the study effects are stronger. 

# o	The study effects and population effects are difficult to distinguish because 
# the study effects appear in the first principal component.

# Question 10

# 10.	A. Load the Bodymap data with the following command

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
ls()
bm = bodymap.eset
head(bm)
edata = exprs(bm)
head(edata)
pdata_bm=pData(bm)
head(pdata)

# Set the seed using the command set.seed(33353) then estimate a single surrogate 
# variable using the sva function after log2(data + 1) transforming the expression 
# data, removing rows with rowMeans less than 1, and treating age as the outcome 
# (hint: you may have to subset the expression data to the samples without missing 
# values of age to get the model to fit). What is the correlation between the 
# estimated surrogate for batch and age? Is the surrogate more highly correlated 
# with race or gender?
  
# 10. B. Loading Libraries
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

# 10. C. Preprocessing the data
set.seed(33353)
pheno = na.omit(pdata_bm)
edata = edata[,rownames(pheno), drop=FALSE]
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 1,]

# 10. D. Fit a sva model
mod = model.matrix(~age, data=pheno)
mod
mod0 = model.matrix(~1, data=pheno)
mod0
sva1 = sva(edata, mod,mod0, n.sv=1)

# 10. E. Correlation between surrogate for batch and age
cor(sva1$sv, pheno$age)

# 10. F. Correlation between surrogate for batch and race
cor(sva1$sv, as.numeric(pheno$race))

# 10. G. Correlation between surrogate for batch and gender
cor(sva1$sv, as.numeric(pheno$gender))

# o	Correlation with age: 0.99 More highly correlated with race.

# o	Correlation with age: 0.20 More highly correlated with gender. 

# o	Correlation with age: 0.99 More highly correlated with gender. 

# o	Correlation with age: 0.33 More highly correlated with gender. 






