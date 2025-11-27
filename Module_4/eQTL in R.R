# eQTL in R

# 1. Installing Bioconductor package
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()

# 2. Installing MatrixEQTL
# BiocManager::install(c("MatrixEQTL"), force = TRUE)

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Plotting Parameters
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
par(pch=19)

# 5. Load Libraries
library(devtools)
library(Biobase)
library(MatrixEQTL)

# 6. Download the Data
base.dir = find.package("MatrixEQTL")
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="")
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="")
output_file_name = tempfile()

# 7. Understanding the filepaths to the data
base.dir
paste(base.dir, "/data/SNP.txt", sep="")

# 8. Load and display the expression data using the file paths outlined above
expr = read.table(expression_file_name,sep="\t",
                  header=T,row.names=1)
expr[1,]

# 9. Load and display the SNP data using file paths outlined above
snps = read.table(SNP_file_name,sep="\t",
                  header=T,row.names=1)
snps[1,]

# 10. Load and display the covariates data using the file paths outlined above
cvrt = read.table(covariates_file_name,sep="\t",
                  header=T,row.names=1)
cvrt[1,]

# 11. Extracting and displaying the 1st row of the gene expression and SNP matrix
e1 = as.numeric(expr[1,])
e1
s1 = as.numeric(snps[1,])
s1

# 12. Construct a linear model relating gene expression (e1) to SNPs (s1)
library(broom)
lm1 = lm(e1 ~ s1)
tidy(lm1)

# 13. Plot the relationship between gene expression level (e1) and genotype SNP (s1)
plot(e1 ~ jitter(s1),
     col=(s1+1),xaxt="n",xlab="Genotype",ylab="Expression")
axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
lines(lm1$fitted ~ s1,type="b",pch=15,col="darkgrey")

# 14. Dimensions of Expression data
dim(expr)

# 15. Setting up parameters for constructing multiple linear regression models

# 15. A. P - value threshold
pvOutputThreshold = 1e-2

# 15. B. Setting up an independence error model
errorCovariance = numeric()

# 15. C. Setting up an additive linear model vs dominant or recessive model
useModel = modelLINEAR

# 16. Setting up files for the MatrixEQTL object to make it computationally efficient (SNP data)
snps = SlicedData$new()
snps$fileDelimiter = "\t"                   # the TAB character
snps$fileOmitCharacters = "NA"              # denote missing values;
snps$fileSkipRows = 1                       # one row of column labels
snps$fileSkipColumns = 1                    # one column of row labels
snps$fileSliceSize = 2000                   # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

# 17. Setting up files for the MatrixEQTL object to make it computationally efficient (Gene Expression)
gene = SlicedData$new()
gene$fileDelimiter = "\t"                   # the TAB character
gene$fileOmitCharacters = "NA"              # denote missing values;
gene$fileSkipRows = 1                       # one row of column labels
gene$fileSkipColumns = 1                    # one column of row labels
gene$fileSliceSize = 2000                   # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name)

# 18. Create an empty object for adjustment covariates the model will have no covariates
cvrt = SlicedData$new()

# 19. Running the MatrixEQTL
me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = NULL,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE, 
  noFDRsaveMemory = FALSE)

# 20. Plot all the p-values from the tests
plot(me)

# 21. Displaying the 'me' object
names(me)

# 22. How long it took to run the object
me$time.in.sec

# 23. Description of parameters
me$param

# 24. Number of tests performed
names(me$all)
me$all$ntests

# 25. Number of EQTLs it found, the number that passed the significance test
me$all$neqtls

# 26. Displaying the EQTLs
me$all$eqtls


