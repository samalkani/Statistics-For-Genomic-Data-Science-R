# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing "sva", "bladderbatch", "snpStats"
# BiocManager::install(c("sva", "bladderbatch", "snpStats"), type = "source", force = TRUE)

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load Libraries
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

# 5. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 6. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 7. Making the circles on the plots filled solid
par(pch=19)

# 8. Load data from snpStats package
data(for.exercise)

# 9. The controls dataset
controls <- rownames(subject.support)[subject.support$cc==0]

# 10. Dimensions
dim(controls)
length(controls)

# 11. Display the controls dataset
head(controls)

# 12. Take every 10th value in the dataset
use <- seq(1, ncol(snps.10), 10)

# 13. Take a subset controls / disease
ctl.10 <- snps.10[controls,use]

# 14. Principal components

# 14. A. Linear algebra calculation using xxt command
xxmat <- xxt(ctl.10, correct.for.missing=FALSE)

# 14. B. Eigen decomposition to the above calculation
evv <- eigen(xxmat, symmetric=TRUE)

# 14. C. Calculate the Eigenvectors
pcs <- evv$vectors[,1:5]

# 14. D. Dimensions of PC matrix
dim(pcs)

# 14. E. 1st row of PC matrix
pcs[1,]

# 14. F. Look at the population the PC's come from
pop <- subject.support[controls,"stratum"]

# 14. G. Dimension
dim(pop)

# 14. H. Display pop variable
head(pop)
length(pop)

# 14. I. Plot PC1 vs PC2
par(mfrow=c(1,1))
plot(pcs[,1],pcs[,2],col=as.numeric(pop),
     xlab="PC1",ylab="PC2")
legend(0,0.15,legend=levels(pop),pch=19,col=1:2)










