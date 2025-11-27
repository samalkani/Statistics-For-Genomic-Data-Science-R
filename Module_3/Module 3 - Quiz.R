# Module 3 - Quiz

# Question 1

# 1. A.	Load the example SNP data with the following code:
  
library(snpStats)
library(broom)
ls()
data(for.exercise)
snps.10
use <- seq(1, ncol(snps.10), 10)
head(use)
sub.10 <- snps.10[,use]
sub.10
snpdata = sub.10@.Data
head(subject.support)
status = subject.support$cc
head(status)

# Fit a linear model and a logistic regression model to the data for the 3rd SNP. 
# What are the coefficients for the SNP variable? How are they interpreted? (Hint: 
# Don't forget to recode the 0 values to NA for the SNP data)

# 1. B. Recode 0 values to NA
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA

# 1. C. Fit a linear model
lm3 = lm(status ~ snp3)
tidy(lm3)


# 1. D. Fit a logistic regression model
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3)

# o	Linear Model =  -0.16

#   Logistic Model = -0.04

# Both models are fit on the additive scale. So in the linear model case, the 
# coefficient is the decrease in probability associated with each additional copy 
# of the minor allele. In the logistic regression case, it is the decrease in the 
# log odds ratio associated with each additional copy of the minor allele.

# o	Linear Model = -0.04

#   Logistic Model = -0.16

# Both models are fit on the additive scale. So in the linear model case, the 
# coefficient is the decrease in probability associated with each additional copy 
# of the minor allele. In the logistic regression case, it is the decrease in the 
# log odds ratio associated with each additional copy of the minor allele.

# o	Linear Model = 0.54

#   Logistic Model = 0.18

# Both models are fit on the additive scale. So in the linear model case, the 
# coefficient is the decrease in probability associated with each additional 
# copy of the minor allele. In the logistic regression case, it is the decrease 
# in the log odds ratio associated with each additional copy of the minor allele.

# Linear Model = 0.54

# Logistic Model = 0.18

# Both models are fit on the additive scale. So in both cases the coefficient is 
# the decrease in probability associated with each additional copy of the minor 
# allele.

# 2.	In the previous question why might the choice of logistic regression be better 
# than the choice of linear regression?
  
# o	The linear model only allows modeling relationships on the additive scale 
# but we might want to consider a dominant or recessive model. 

# o	The log odds is always more interpretable than a change in probability on the
# additive scale. 

# o	It is customary to use logistic regression for case-control data like those 
# obtained from genome-wide association studies.

# o	If you included more variables it would be possible to get negative estimates 
# for the probability of being a case from the linear model, but this would be 
# prevented with the logistic regression model.

# Question 3

# 3. A. Load the example SNP data with the following code:
  
library(snpStats)
library(broom)
ls()
data(for.exercise)
snps.10
use <- seq(1, ncol(snps.10), 10)
head(use)
sub.10 <- snps.10[,use]
sub.10
snpdata = sub.10@.Data
head(subject.support)
status = subject.support$cc
head(status)

# Fit a logistic regression model on a recessive (need 2 copies of minor allele 
# to confer risk) and additive scale for the 10th SNP. Make a table of the fitted 
# values versus the case/control status. Does one model fit better than the other?
  
# 3. B. fit a logistic regression model
snp10 = as.numeric(snpdata[,10])
snp10[snp10==0] = NA
glm10 = glm(status ~ snp10, family="binomial")
tidy(glm10)

# 3. C. Recessive Model
snp10_dom = (snp10 == 2)
glm10_dom = glm(status ~ snp10_dom, family="binomial")
tidy(glm10_dom)

# 3. D. Fitted values vs case/control status for logistic regression (additive model)
table(glm10$y, glm10$fitted.values)

# 3. E. Fitted values vs case/control status (Recessive model)
table(glm10_dom$y, glm10_dom$fitted.values)

# o	The recessive model fits much better since it appears that once you aggregate 
# the heterozygotes and homozygous minor alleles, there is a bigger difference in 
# the proportion of cases and controls. 

# o	The additive model fits much better since there are fewer parameters to fit 
# and the effect size is so large. 

# o	No, in all cases, the fitted values are near 0.5 and there are about an equal 
# number of cases and controls in each group. This is true regardless of whether 
# you fit a recessive or additive model. 

# o	The recessive model shows a strong effect, but the additive model shows no 
# difference so the recessive model is better.

# Question 4

# 4. A.	Load the example SNP data with the following code:
  
library(snpStats)
library(broom)
ls()
data(for.exercise)
snps.10
use <- seq(1, ncol(snps.10), 10)
head(use)
sub.10 <- snps.10[,use]
sub.10
snpdata = sub.10@.Data
head(subject.support)
status = subject.support$cc
head(status)

# Fit an additive logistic regression model to each SNP. What is the average 
# effect size? What is the max? What is the minimum?
  
# 4. B. Fit an additive logistic regression model to each SNP
results = rep(NA, dim(snpdata)[2])
for (i in 1:ncol(snpdata)){
  snpdata_i = as.numeric(snpdata[,i])
  snpdata_i[snpdata_i == 0] = NA
  glm_i = glm(status ~ snpdata_i, family = "binomial")
  results[i] = tidy(glm_i)$statistic[2]
}

# 4. C. Average effect size
mean(results)

# 4. D. Minimum effect size
min(results)

# 4. E. maximum effect size
max(results)


# o	Average effect size =  0.02, minimum = -0.88, maximum = 0.88

# o	Average effect size =  0.007, minimum = -4.25, maximum = 3.90

# o	Average effect size =  -0.02, minimum =-3.59 , maximum = 4.16

# o	Average effect size =  1.35, minimum =-6.26 , maximum = 6.26

# Question 5

# 5. A. Load the example SNP data with the following code:
  
library(snpStats)
library(broom)
ls()
data(for.exercise)
snps.10
use <- seq(1, ncol(snps.10), 10)
head(use)
sub.10 <- snps.10[,use]
sub.10
snpdata = sub.10@.Data
head(subject.support)
status = subject.support$cc
head(status)

# Fit an additive logistic regression model to each SNP and square the coefficients. 
# What is the correlation with the results from using snp.rhs.tests and chi.squared? 
# Why does this make sense?
  
# 5. B. Square the coefficients
results_coeff_squre =  results^2

# 5. C. Correlation with the results from using snp.rhs.tests and chi.squared
glm_all = snp.rhs.tests(status ~ 1, snp.data = sub.10)
cor(results_coeff_squre, chi.squared(glm_all))

# o	0.99. They are both testing for the same association using the same additive 
# regression model on the logistic scale. But it doesn't make sense since they 
# should be perfectly correlated. 

# o	0.99. They are both testing for the same association using the same additive
# regression model on the logistic scale but using slightly different tests. 

# o	0.002 It doesn't make sense since they are both testing for the same association 
# using the same additive regression model on the logistic scale. But it doesn't 
# make sense since they should be perfectly correlated. 

# o	0.81 They are both testing for the same association using the same additive 
# regression model on the logistic scale. But it doesn't make sense since they 
# should be perfectly correlated.

# Question 6

# 6. A.	Load the Montgomery and Pickrell eSet:

library(devtools)
library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
ls()
mp = montpick.eset
mp
pdata=pData(mp)
head(pdata)
edata=as.data.frame(exprs(mp))
head(edata)
fdata = fData(mp)
head(fdata)

# Do the log2(data + 1) transform and fit calculate F-statistics for the difference 
# between studies/populations using genefilter:rowFtests and using genefilter:rowttests.
# Do you get the same statistic? Do you get the same p-value?
  
# 6. B. Log2 transform
edata = log2(as.matrix(edata) + 1)

# 6. C. Perform rowttests
library(genefilter)
tstats_obj = rowttests(edata, as.factor(pdata$population))
tidyr::tibble(tstats_obj)

# 6. D. Perform rowFtests
fstats_obj = rowFtests(edata, as.factor(pdata$population))
tidyr::tibble(fstats_obj)

# 6. E. Histogram plots of T-tests and F-Tests
par(mfrow=c(1,2))
hist(tstats_obj$statistic, col=2)
hist(fstats_obj$statistic, col=2)

# o You get the same p-values and statistics. This is because the F-statistic and
# t-statistic are the exact same in this case. 

# o	You get different p-values and statistics. The F-statistic and t-statistic are 
# testing the same thing but do it totally differently. 

# o	You get the same p-value but different statistics. This is because the F-statistic 
# and t-statistic test the same thing when doing a two group test and one is a transform 
# of the other. 

# o	You get different p-values but the same statistic. This is because the F-statistic 
# and t-statistic test the same thing when doing a two group test and one is a transform 
# of the other.

# Question 7

# 7. A. Load the Montgomery and Pickrell eSet:
  
library(devtools)
library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
ls()
mp = montpick.eset
mp
pdata=pData(mp)
head(pdata)
edata=as.data.frame(exprs(mp))
head(edata, n = 1)
fdata = fData(mp)
head(fdata)

# First test for differences between the studies using the DESeq2 package using 
# the DESeq function. Then do the log2(data + 1) transform and do the test for 
# differences between studies using the limma package and the lmFit, ebayes and 
# topTable functions. What is the correlation in the statistics between the two 
# analyses? Are there more differences for the large statistics or the small 
# statistics (hint: Make an MA-plot).

# 7. B.Load Libraries
library(DESeq2)
library(limma)
library(edge)
library(genefilter)

# 7. C. Using DESeq2 test the differences between the studies
de = DESeqDataSetFromMatrix(edata, pdata, design =~study)
glm_de = DESeq(de)
result_de = results(glm_de)
result_de

# 7. D. Using limma test the differences
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ as.factor(pdata$study))
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma) 
top = topTable(ebayes_limma,number=dim(edata)[1], sort.by="none")
head(top)

# 7. E. Correlation in the statistics between two analyses
cor(result_de$stat, top$t)

# 7. F. Make an MA-plot
y = cbind(result_de$stat, top$t)
limma::plotMA(y)

# o	0.75. There are more differences for the large statistics.

# o	0.63. There are more differences for the large statistics.

# o	0.93. There are more differences for the small statistics. 

# o	0.36. There are more differences for the large statistics.

# Question 8

# 8. A. DESeq analysis
fp_bh = p.adjust(result_de$pvalue, method="BH")
sum(fp_bh < 0.05)

# 8. B. limma analysis
fp_bh = p.adjust(top$P.Value, method="BH")
sum(fp_bh < 0.05)


# o	DESeq = 1995 significant; 

# limma = 2807 significant

# o	DESeq = 12 significant; 

# limma = 3significant

# o	DESeq = 0 significant; 

# limma = 0 significant

# o	DESeq = 1119 significant; 

# limma = 2328 significant

# Question 9

# Is the number of significant differences surprising for the analysis comparing 
# studies from Question 8? Why or why not?
  
# o	Yes. This is a very large number of genes different between studies and we 
#   don't have a good explanation.

# o	Yes and no. It is surprising because there is a large fraction of the genes 
#   that are significantly different, but it isn't that surprising because we would 
#   expect that when comparing measurements from very different batches. 

# o	No. We are testing many genes so we expect a huge fraction to be different 
#   between studies. 

# o	No. There are very few genes different between studies and that is what we 
#   would expect.


