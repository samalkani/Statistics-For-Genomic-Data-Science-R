# Clustering

# 1. Installing Bioconductor packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(version = "3.22")

# 2. Installing dendextend
BiocManager::install(c("dendextend"), force = TRUE)

# 3. Troubleshooting installing packages
# BiocManager::valid()

# 4. Load Libraries
library(devtools)
library(Biobase)
library(dendextend)

# 5. Making the plots prettier
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")

# 6. Use the palette command to direct R to use those colors outlined above
palette(tropical)

# 7. Making the circles on the plots filled solid
par(pch=19)

# 8. Load Data from URL
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)

# 9. Body Map eset
bm = bodymap.eset

# 10. Phenotype Data
pdata=pData(bm)

# 11. Gene Expression Data
edata=as.data.frame(exprs(bm))

# 12. Feature Data
fdata = fData(bm)
ls()

# 13. Filter out lowly expressed genes
edata = edata[rowMeans(edata) > 5000,]

# 14. Dimensions
dim(edata)

# 15. Log2 transform
edata = log2(edata + 1)

# 16. By default calculates the distance between rows
dist1 = dist(t(edata))
dist1

# 17. Look at distance matrix
colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)

# 18. Apply hierarchical clustering on distance matrix
hclust1 = hclust(dist1)

# 19. Plot the dendrogram
plot(hclust1)

# 20. Re-plot dendrogram with the labels at the same level
plot(hclust1,hang=-1)

# 21. Add colour to the Dendrogram labels
dend = as.dendrogram(hclust1)	
dend = color_labels(hclust1,4,col=1:4)
plot(dend)

# 22. Add colour to the Dendrogram labels three clusters
dend = as.dendrogram(hclust1)
dend = color_labels(hclust1,3,col=1:3)
plot(dend)

# 23. Plot dendrogram with pre-specified groups
labels_colors(dend) = c(rep(1,10),rep(2,9))
plot(dend)

# 24. K-means clustering and the clustering object
kmeans1 = kmeans(edata,centers=3)
names(kmeans1)

# 25. Plotting the cluster centres
matplot(t(kmeans1$centers),col=1:3,type="l",lwd=3)

# 26. How many genes belong to each cluster
table(kmeans1$cluster)

# 27. Which cluster each gene belongs to
kmeans1$cluster[1:10]

# 28. Re-ordering the rows according to the cluster relationship
newdata = as.matrix(edata)[order(kmeans1$cluster),]
dim(newdata)

# 29. Create Heatmap
heatmap(newdata, col = colramp, Colv = NA, Rowv = NA)

# 30. Repeat clustering with same data but get different result
kmeans2 = kmeans(edata,centers=3)
table(kmeans1$cluster,kmeans2$cluster)








