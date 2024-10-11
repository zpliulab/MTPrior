#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")

##### Input Data #####
HCC_mRNA1<- read.delim("DATA/HCC/HCC_mRNA_genename.txt", header = TRUE, sep = "\t")
dim(HCC_mRNA1)
HCC_mRNA1[1:10,1:10]
## ## ## ##
expression_data<- HCC_mRNA1[,-1]
rownames(expression_data)<- HCC_mRNA1[,1]
dim(expression_data)

##### Normalize the gene expression data #####
norm_factors <- calcNormFactors(expression_data)
normalized_data <- t(t(expression_data) / norm_factors)

##### Log-transform the expression values #####
log_transformed_data <- log2(1 + normalized_data)
dim(log_transformed_data)


##### Select the tumor data not the normal data #####
HCC_sample_type<- read.delim("DATA/HCC/HCC_sample_type.txt", header = TRUE, sep = "\t")
dim(HCC_sample_type)
unique(HCC_sample_type$LIHC_mRNA.sample_type)
#dim(HCC_sample_type[HCC_sample_type$LIHC_mRNA.sample_type %in% "Solid Tissue Normal",])
row_numbers <- which(HCC_sample_type$LIHC_mRNA.sample_type != "Solid Tissue Normal")
expression_data_T<- log_transformed_data[,row_numbers]
dim(expression_data_T)


##Normal data part
Normal_rows<- which(HCC_sample_type$LIHC_mRNA.sample_type == "Solid Tissue Normal")
expression_data_Normal<- log_transformed_data[,Normal_rows]
dim(expression_data_Normal)




##### Filter out low-expression genes #####
# Calculate total counts for each sample
total_counts <- colSums(expression_data_T)
# Perform total count normalization
normalized_data <- t(t(expression_data_T) / total_counts)
##### Set threshold and min_samples values
threshold <- 0.5 # CPM threshold for RNA-seq data
min_samples <- 3 # Minimum number of samples with expression above threshold
expression_data_filtered <- normalized_data[rowSums(cpm(normalized_data) > threshold) >= min_samples, ]
dim(expression_data_filtered)




#####scale the expression values ###### 
##it can be  z-score normalization  or  min-max scaling
#expression_data_scaled <- scale(log_transformed_data)
## Define a function to scale data between 0 and 1
scale_to_01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# Apply the function to each column of the data matrix
expression_data_scaled <- apply(expression_data_filtered, 2, scale_to_01)
dim(expression_data_scaled)
expression_data_scaled[1:10,1:10]
## ## ## ##



####### Filter according variant in Genes #####
# Calculate variance of gene expression across samples
gene_variance <- apply(expression_data_scaled, 1, var)
# Set a threshold for variance
hist(gene_variance, breaks = 30, main = "Distribution of Gene Variances", xlab = "Variance", ylab = "Frequency")

threshold <- 0.005 # Adjust the threshold as needed

# Identify genes with variance above the threshold
genes_above_threshold <- names(gene_variance[gene_variance >= threshold])
length(genes_above_threshold)
# Filter gene expression data based on variance threshold
filtered_gene_expression_data <- expression_data_scaled[genes_above_threshold, ]
dim(filtered_gene_expression_data)




####Test if all Disease gene are in filtered data or not
HCC_Dgeans_A<- read.table("DATA/HCC/HCC_Diseasegeans_DiGeNet_GWAS.txt", header = F)
head(HCC_Dgeans_A)
dim(HCC_Dgeans_A)
filtered_genes<- data.frame(G=rownames(filtered_gene_expression_data))
length(filtered_genes[filtered_genes$G %in% HCC_Dgeans_A$V1,])
## Common Gene between HCC Disease Genes and filtered Genes
commonG<-filtered_genes[filtered_genes$G %in% HCC_Dgeans_A$V1,]
length(commonG)
## Genes in Disease Gene that they are not in filtered data
RemainDeGe<- setdiff(HCC_Dgeans_A$V1, commonG)
length(RemainDeGe)
## Rows for Remaining Genes
extraRows <- subset(expression_data_T, rownames(expression_data_T) %in% RemainDeGe)
dim(extraRows)
extraRows_scaled <- apply(extraRows, 2, scale_to_01)
dim(extraRows_scaled)
extraRows_scaled[1:10,1:10]
dim(filtered_gene_expression_data[rownames(filtered_gene_expression_data) %in% rownames(extraRows_scaled),])



write.table(filtered_gene_expression_data,"DATA/filteredGeneExpression.txt",
            quote=F,sep="\t",row.names=FALSE,col.names = FALSE)