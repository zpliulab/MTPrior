#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")


## Meta-paths
## lncRNA-mRNA-lncRNA
## lncRNA-miRNA-lncRNA
## lncRNA-miRNA-mRNA-miRNA-lncRNA


lncRNA_IDs <- read.table("DATA/HCC/lncRNAData/filteredLncNames.txt", header = FALSE, sep = "\t")
lnc_All<- lncRNA_IDs$V1


###### lncRNA-mRNA-lncRNA Metapath Matrix ########
lncRNA_mRNA_data <- read.table("DATA/HCC/lncRNAData/f_lncRNA_mRNA.txt", header = FALSE, sep = "\t")
head(lncRNA_mRNA_data)
colnames(lncRNA_mRNA_data)<- c("lncRNA","mRNA")
mRNA_lncRNA_data<- lncRNA_mRNA_data[,c("mRNA","lncRNA")]

lnc_m_lnc <- merge(lncRNA_mRNA_data, mRNA_lncRNA_data, by.x = "mRNA", by.y = "mRNA", all.x = FALSE,all.y=FALSE)
lnc_m_lnc<-lnc_m_lnc[,c("lncRNA.x","lncRNA.y")]
lnc_m_lnc<-unique(lnc_m_lnc)
dim(lnc_m_lnc)


lnc_m_lnc_matrix <- matrix(0, nrow = length(lnc_All), ncol = length(lnc_All), dimnames = list(lnc_All, lnc_All))
dim(lnc_m_lnc_matrix)

# Identify edges in lnc_m_lnc and set the corresponding cells in the interaction matrix to 1
edges <- as.matrix(lnc_m_lnc[,c("lncRNA.x","lncRNA.y")])
lnc_m_lnc_matrix[edges] <- 1
sum(lnc_m_lnc_matrix>0)

write.table(lnc_m_lnc_matrix,"DATA/HCC/lncRNAData/MP_HCC_lncmlnc_matrix_2.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)


###### lncRNA-miRNA-lncRNA Metapath Matrix ########
miRNA_lncRNA_data <- read.table("DATA/HCC/lncRNAData/f_lncRNA_miRNA.txt", header = FALSE, sep = "\t")
head(miRNA_lncRNA_data)
colnames(miRNA_lncRNA_data)<- c("miRNA","lncRNA")
lncRNA_miRNA_data<- miRNA_lncRNA_data[,c("lncRNA","miRNA")]
head(lncRNA_miRNA_data)


lnc_mi_lnc <- merge(lncRNA_miRNA_data, miRNA_lncRNA_data, by.x = "miRNA", by.y = "miRNA", all.x = FALSE,all.y=FALSE)
lnc_mi_lnc<-lnc_mi_lnc[,c("lncRNA.x","lncRNA.y")]
lnc_mi_lnc<-unique(lnc_mi_lnc)
dim(lnc_mi_lnc)


lnc_mi_lnc_matrix <- matrix(0, nrow = length(lnc_All), ncol = length(lnc_All), dimnames = list(lnc_All, lnc_All))
dim(lnc_mi_lnc_matrix)
# Identify edges in lnc_m_lnc and set the corresponding cells in the interaction matrix to 1
edges <- as.matrix(lnc_mi_lnc[,c("lncRNA.x","lncRNA.y")])
lnc_mi_lnc_matrix[edges] <- 1
sum(lnc_mi_lnc_matrix>0)

write.table(lnc_mi_lnc_matrix,"DATA/HCC/lncRNAData/MP_HCC_lncmilnc_matrix_2.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)



###### lncRNA-miRNA-mRNA-miRNA-lncRNA Metapath Matrix ########
miRNA_mRNA_data <- read.table("DATA/HCC/lncRNAData/f_miRNA_mRNA.txt", header = FALSE, sep = "\t")
head(miRNA_mRNA_data)
colnames(miRNA_mRNA_data)<- c("miRNA","mRNA")
mRNA_miRNA_data<- miRNA_mRNA_data[,c("mRNA","miRNA")]
head(mRNA_miRNA_data)


lnc_mi_m <- merge(lncRNA_miRNA_data, miRNA_mRNA_data, by.x = "miRNA", by.y = "miRNA", all.x = FALSE,all.y=FALSE)
lnc_mi_m<-lnc_mi_m[,c("lncRNA","mRNA")]
lnc_mi_m<-unique(lnc_mi_m)
dim(lnc_mi_m)

lnc_mi_m_mi <- merge(lnc_mi_m, mRNA_miRNA_data , by.x = "mRNA", by.y = "mRNA", all.x = FALSE,all.y=FALSE)
lnc_mi_m_mi<-lnc_mi_m_mi[,c("lncRNA","miRNA")]
lnc_mi_m_mi<-unique(lnc_mi_m_mi)
dim(lnc_mi_m_mi)


lnc_mi_m_mi_lnc <- merge(lnc_mi_m_mi, miRNA_lncRNA_data , by.x = "miRNA", by.y = "miRNA", all.x = FALSE,all.y=FALSE)
lnc_mi_m_mi_lnc<-lnc_mi_m_mi_lnc[,c("lncRNA.x","lncRNA.y")]
lnc_mi_m_mi_lnc<-unique(lnc_mi_m_mi_lnc)
dim(lnc_mi_m_mi_lnc)


lnc_mi_m_mi_lnc_matrix <- matrix(0, nrow = length(lnc_All), ncol = length(lnc_All), dimnames = list(lnc_All, lnc_All))
dim(lnc_mi_m_mi_lnc_matrix)

# Identify edges in lnc_m_lnc and set the corresponding cells in the interaction matrix to 1
edges <- as.matrix(lnc_mi_m_mi_lnc[,c("lncRNA.x","lncRNA.y")])
lnc_mi_m_mi_lnc_matrix[edges] <- 1
sum(lnc_mi_m_mi_lnc_matrix>0)

write.table(lnc_mi_m_mi_lnc_matrix,"DATA/HCC/lncRNAData/MP_HCC_lncmimmilnc_matrix_2.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)






##Counts
dim(miRNA_lncRNA_data)
dim(miRNA_mRNA_data)
dim(lncRNA_mRNA_data)

d1<-unique(data.frame(m=miRNA_mRNA_data$mRNA))
d2<- unique(data.frame(m=lncRNA_mRNA_data$mRNA))
dim(unique(rbind(d1,d2))) #mRNA count

d1<-unique(data.frame(m=miRNA_mRNA_data$miRNA))
d2<- unique(data.frame(m=miRNA_lncRNA_data$miRNA))
dim(unique(rbind(d1,d2))) #miRNA count








