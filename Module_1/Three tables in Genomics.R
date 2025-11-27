# Three tables in genomics (in R)

# 1. Installing Biobase package

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")

BiocManager::install("Biobase")
library(Biobase)


# 2. Load the data
con=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
ls()

# 3. Close the connection
close(con)

# 4. Shorten the variable name for ease of analysis
bm = bodymap.eset

# 5. Display the expression set
bm

# 6. Extract out the gene expression data
exp_data = exprs(bm)


# 7. Displaying dimensions of the matrix
dim(exp_data)

# 8. Look at the top 5 rows of matrix
head(exp_data,n=5)

# 9. Look at the other two tables, first phenotype data
pheno_data = pData(bm)
dim(pheno_data)

# 10. Display the top five rows of the phenotype data
head(pheno_data)

# 11. Look at the second table the feature table
feature_data = fData(bm)
dim(feature_data)

# 12. Display the first 15 rows of the feature table
feature_data[1:15,1]




