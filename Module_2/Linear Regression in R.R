# Linear Regression in R

# 1. Install Packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) 

#   install.packages("BiocManager")

# BiocManager::install(version = "3.22")

install.packages("broom")

# 2. Load Libraries
library(devtools)
library(Biobase)
library(broom)

# 3. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 4. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 5. Making the circles on the plots filled solid
par(pch=19)

# 6. Download Data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
head(bm)
pdata=pData(bm)
head(pdata)
edata=as.data.frame(exprs(bm))
head(edata)
fdata = fData(bm)
head(fdata)
ls()

# 7. Fitting a linear regression model genes vs age - Model 1

# 7. A. Convert data set to matrix (easier to deal with numeric values)
edata = as.matrix(edata)

# 7. B. Fit model (on first gene in matrix - ENSG00000000003)
lm1 = lm(edata[1,] ~ pdata$age)

# 7. C. Display model
lm1

# 7. D. Information about model
tidy(lm1)

# 7. E. Plot age vs the first gene
plot(pdata$age,edata[1,], col=1)

# 7. F. Fit linear model to this data
abline(lm1, col=2,lwd=3)

# 8. Look at Gender variable, a factor variable
pdata$gender

# 9. Count data for Gender
table(pdata$gender)

# 10. Boxplot of 1st gene vs Gender
boxplot(edata[1,] ~ pdata$gender)

# 11. Boxplot of 1st gene vs Gender overlaid with data points &
# linear model to relate categorical variable to outcome variable
points(edata[1,] ~ jitter(as.numeric(pdata$gender)), col=as.numeric(pdata$gender))

# 12. Creating male dummy variable
dummy_m = pdata$gender=="M"
dummy_m

# 13. Multiply dummy variable by 1 gives 1's and 0's
dummy_m*1

# 14. Creating female dummy variable
dummy_f = pdata$gender=="F"
dummy_f

# 15. Categorical Variables - Gender – Model 2

# 15. A. R modelling formulae create dummy variables automatically for categorical variables
lm2 = lm(edata[1,] ~ pdata$gender)
tidy(lm2)

# 15. B. Peaking under the hood to look at the variables used in the model
mod2 = model.matrix(~pdata$gender)
mod2

# 16. Categorical variables - Tissue Type - Model 3

# 16. A. Displaying the levels of the tissue type categorical variable
table(pdata$tissue.type)

# 16. B. You can define multiple dummy variables one for each tissue type
pdata$tissue.type == "adipose"

# 16. C. R modelling formulae create coefficients for all levels of a categorical variable with multiple levels
tidy(lm(edata[1,] ~ pdata$tissue.type ))

# 16. D. Adjusting for variables
lm3 = lm(edata[1,] ~ pdata$age + pdata$gender)
tidy(lm3)

# 16. E. Fitting interaction terms in the linear model
lm4 = lm(edata[1,] ~ pdata$age*pdata$gender)
tidy(lm4)

# 17. Fitting a linear regression model gene 6 vs age - Model 4

# 17. A. Constructing a linear model for gene 6 vs age
lm4 = lm(edata[6,] ~ pdata$age)

# 17. B. Plot
plot(pdata$age,edata[6,],col=2)

# 17. C. Plot Regression line
abline(lm4,col=1,lwd=3)

# 18. Fitting model of gene 6 across 19 samples - Model 5

# 18. A.  Plot of gene 6 across the 19 samples (indexes)
index = 1:19
lm5 = lm(edata[6,] ~ index)
plot(index,edata[6,],col=2)

# 18. B. Fit the model line to plot
abline(lm5,col=1,lwd=3)

# 19. Fitting model without 19th point - Model 6

# 19. A. Re-fit the model without the 19th data point
lm6 = lm(edata[6,-19] ~ index[-19])
abline(lm6,col=3,lwd=3)

# 19. B. Adding a legend to the plot
legend(5,1000,c("With outlier","Without outlier"),col=c(1,3),lwd=3)

# 19. C. Histogram of residuals of the models with & without outlier
par(mfrow=c(1,2))
hist(lm6$residuals,col=2)
hist(lm5$residuals,col=3)

# 20. Log Transform the data and look at the residuals again - Model 7
gene1 = log2(edata[1,]+1)
lm7 = lm(gene1 ~ index)
hist(lm7$residuals,col=4)

# 21. Too many covariates – Model 8

# 21. A. Fitting a model with too many covariates
lm8 = lm(gene1 ~ pdata$tissue.type + pdata$age)
tidy(lm8)

# 21. B. Dimensions of the model matrix you're fitting
dim(model.matrix(~ pdata$tissue.type + pdata$age))

# 22. Plot a linear model for age and plot the residuals - Model 9
colramp = colorRampPalette(1:4)(17)
lm9 = lm(edata[2,] ~ pdata$age)
plot(lm9$residuals,col=colramp[as.numeric(pdata$tissue.type)])

























