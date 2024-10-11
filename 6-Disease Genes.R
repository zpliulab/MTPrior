setwd("/home/fatemeh/Fatemeh/0--SecoundProject")

##data from GWAS Catalog
##https://www.ebi.ac.uk/gwas/docs/fileheaders#_file_headers_for_catalog_v1_0
data <- read.delim2("DATA/HCC/gwas_catalog_v1.0-associations_e111_r2024-02-11.tsv", header = TRUE, sep = "\t")

##### Find rows contain the specific disease #####
D_data<-data.frame(data$DISEASE.TRAIT)
dim(D_data)
D_data_U<-data.frame(unique(D_data))
result <- D_data_U[grep("Hepatocellular", D_data_U$data.DISEASE.TRAIT), ]
result
data_filter<- data[data$DISEASE.TRAIT %in% result,]
dim(data_filter)

##### select related columns #####
colnames(data_filter)
columns <- c("REPORTED.GENE.S.", "MAPPED_GENE","DISEASE.TRAIT")
gene_related_D<-data_filter[,columns]
gene_related_D

##### Clear data #####
gene_related_D$MAPPED_GENE <- gsub(";", ",", gene_related_D$MAPPED_GENE)
gene_related_D$MAPPED_GENE <- gsub(" - ", ",", gene_related_D$MAPPED_GENE)

dim(gene_related_D) 
gene_related_D

##### Split values in the first and second columns ##### 
# split_column1 <- strsplit(gene_related_D$REPORTED.GENE.S., ",")
# split_column2 <- strsplit(gene_related_D$MAPPED_GENE, ",")
# 

split_rows1 <- strsplit(gene_related_D[["REPORTED.GENE.S."]], ", ")
gene_Values1 <-data.frame(genes=unlist(split_rows1)) 
gene_Values1

split_rows2 <- strsplit(gene_Values1[["genes"]], "or")
gene_Values2 <-data.frame(genes=unlist(split_rows2)) 
gene_Values2

split_rows3 <- strsplit(gene_related_D[["MAPPED_GENE"]], ",")
gene_Values3 <-data.frame(genes=unlist(split_rows3)) 
gene_Values3

result3 <- gene_Values3[grep("No mapped genes", gene_Values3$genes), ]
gene_Values4<-data.frame(genes=gene_Values3[!(gene_Values3$genes %in% result3),])
gene_Values4

gene_Values<- rbind(gene_Values2,gene_Values4)
gene_Values<-data.frame(unique(gene_Values))
gene_Values

gene_Values$genes<- trimws(gene_Values$genes)

write.table(gene_Values, file = "DATA/HCC/HCC_diseasegenes_GWAS.txt",  sep = "\t", quote = FALSE,row.names=FALSE,col.names = FALSE)


#gene_Values<- read.delim2("DATA/HCC/HCC_diseasegenes_GWAS.txt",header=F)

####data from DiGeNet (GWAS & OMIM)
## http://www.disgenet.org/
DisGeNet_Data<-read.delim2("DATA/HCC/DisGeNet/C0279607__C3273019__C0861876__C1711391_disease_gda_evidences.tsv",header = T)
dim(DisGeNet_Data)
DisGeNet_Genes<-data.frame(DisGeNet_Data$Gene)

DisGeNet_Genes<- unique(DisGeNet_Genes)

colnames(gene_Values)<- c("gene")
colnames(DisGeNet_Genes)<- c("gene")
datasetsforDisease<-rbind(DisGeNet_Genes,gene_Values)
datasetsforDisease<-unique(datasetsforDisease)
dim(datasetsforDisease)

####...
####...
####...

HCC_all_genes<- read.delim("DATA/HCC/HCC_all_genes.txt",header = F, sep = "\t")
AllGenesN<-HCC_all_genes$V1

#HCC_Dgeans_inALL1<- data.frame(genes=HCC_all_genes[HCC_all_genes$V1 %in% gene_Values$genes,])
HCC_Dgeans_inALL1<- data.frame(genes=HCC_all_genes[HCC_all_genes$V1 %in% datasetsforDisease$gene,])
dim(HCC_Dgeans_inALL1)


#HCC_Dgeans_inALL<- rbind(HCC_Dgeans_inALL1,HCC_Dgeans_inALL2,HCC_Dgeans_inALL3,HCC_Dgeans_inALL4,HCC_Dgeans_inALL5)
#dim(HCC_Dgeans_inALL1)+dim(HCC_Dgeans_inALL2)+dim(HCC_Dgeans_inALL3)+dim(HCC_Dgeans_inALL4)+dim(HCC_Dgeans_inALL5)

#HCC_Dgeans_inALL<-data.frame(unique(HCC_Dgeans_inALL))
HCC_Dgeans_inALL<-data.frame(unique(HCC_Dgeans_inALL1))
dim(HCC_Dgeans_inALL)

write.table(HCC_Dgeans_inALL, file = "DATA/HCC/HCC_Diseasegeans_DiGeNet_GWAS.txt", sep = "\t", quote = FALSE,row.names=FALSE,col.names = FALSE)

##### Label Matrix contains disease and non-disease genes
lable_G<-HCC_all_genes
lable_G$V2<-0
lable_G$V3<-0
lable_G[lable_G$V1 %in% HCC_Dgeans_inALL$genes,c("V2")]<-1
sum(lable_G$V2)

Not_Dgeans<- data.frame(genes=setdiff(lable_G$V1,HCC_Dgeans_inALL$genes))
dim(Not_Dgeans)

lable_G[lable_G$V1 %in% Not_Dgeans$genes,c("V3")]<-1
sum(lable_G$V3)

dim(lable_G)


#write.table(lable_G, file = "DATA/HCC/Output/lable_G.txt", sep = "\t", quote = FALSE,row.names=FALSE,col.names = FALSE)

lable_Ready<-lable_G[,c(2,3)]
write.table(lable_Ready, file = "DATA/HCC/Output/lable_Ready.txt", sep = "\t", quote = FALSE,row.names=FALSE,col.names = FALSE)








