#Read OTU abundance table
cols <- c("OTU", "W1", "F1", "F2", "W2", "Domain", "Kingdom", 
          "Division", "Class", "Order", "Family", "Genus_Species")
data <- read.table("OTU_table_for_biom_stu.txt", 
                   col.names = cols, skip = 2, fill = TRUE)
head(data)


#Dominance (count to ratio)
for(i in 2:5){
  data[, i] <- data[, i] / sum(data[, i])
}
#Blank to NA
for(i in 6:12){
  data[which(data[, i] == ""), i] <- NA
  data[, i] <- as.factor(data[, i])
}
summary(data)


#PCA
X <- as.matrix(t(data[, 2:5]))
pca <- prcomp(X)
summary(pca)
pca$x


#plot principal component score
plot(pca$x[, c(1, 2)], xlim=c(-0.05, 0.05), ylim=c(-0.025, 0.025), pch=16)
text(pca$x[, 1], pca$x[, 2], rownames(pca$x))
#plot principal component loading (x0.05)
ld_1 <- data.frame(pca$rotation)
rownames(ld_1) <- data$OTU
ld_1.select <- ld_1[order(abs(ld_1$PC1), decreasing = TRUE), ]
arrows(0, 0, ld_1.select$PC1[1:10]*0.05, ld_1.select$PC2[1:10]*0.05, col="black", length = 0.05)
text(ld_1.select$PC1[1:10]*0.05, ld_1.select$PC2[1:10]*0.05, rownames(ld_1.select)[1:10])

#Search OTU
data[data$OTU %in% rownames(ld_1.select)[1:10], ]


#Additional
data.select <- data[data$Class %in% c("Bacillales", "Rhizobiales"), ]
apply(data.select[data.select$Class == "Bacillales", 2:5], 2, sum)
apply(data.select[data.select$Class == "Rhizobiales", 2:5], 2, sum)
