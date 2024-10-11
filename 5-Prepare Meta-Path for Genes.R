
##### Gene-miRNA-Gene metapath matrix #####
miRNA_mRNA_data <- read.table("DATA/HCC/HCC_miRNA_mRNA.txt", header = FALSE, sep = "\t")
dim(miRNA_mRNA_data)
head(miRNA_mRNA_data)
colnames(miRNA_mRNA_data)<- c("miRNA","Gene")
mRNA_miRNA_data<- miRNA_mRNA_data[,c("Gene","miRNA")]

G_mi_G <- merge(mRNA_miRNA_data, miRNA_mRNA_data, by.x = "miRNA", by.y = "miRNA", all.x = FALSE,all.y=FALSE)
G_mi_G1<- G_mi_G
dim(G_mi_G)
head(G_mi_G1)

G_mi_G<-G_mi_G[,c("Gene.x","Gene.y")]
G_mi_G<-unique(G_mi_G)
dim(G_mi_G)
head(G_mi_G)

write.table(G_mi_G,"DATA/HCC/Metapaths/MP_HCC_G_mi_G.txt", quote=F,sep="\t",row.names=FALSE,col.names = TRUE)

AllGenes <- read.table("DATA/HCC/HCC_all_genes.txt", header = FALSE, sep = "\t")
AllGenesN<-AllGenes$V1
length(AllGenesN)

igi_metapath_matrix <- matrix(0, nrow = length(AllGenesN), ncol = length(AllGenesN), dimnames = list(AllGenesN, AllGenesN))

# Identify edges in G_mi_G and set the corresponding cells in the interaction matrix to 1
edges <- as.matrix(G_mi_G[, c("Gene.x", "Gene.y")])
igi_metapath_matrix[edges] <- 1
sum(igi_metapath_matrix>0)

write.table(igi_metapath_matrix,"DATA/HCC/Metapaths/MP_HCC_migenemi_matrix.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)




##### Gene-miRNA-circRNA-miRNA-Gene metapath matrix #####
circRNA_miRNA_data <- read.table("DATA/HCC/HCC_circRNA_miRNA.txt", header = FALSE, sep = "\t")
dim(circRNA_miRNA_data)
head(circRNA_miRNA_data)
colnames(circRNA_miRNA_data)<- c("circRNA","miRNA")
head(circRNA_miRNA_data)
miRNA_circRNA_data<- circRNA_miRNA_data[,c("miRNA","circRNA")]
head(miRNA_circRNA_data)

G_mi_circ <- merge(mRNA_miRNA_data, miRNA_circRNA_data, by.x = "miRNA", by.y = "miRNA", all.x = FALSE,all.y=FALSE)
G_mi_circ1<- G_mi_circ

G_circ <- G_mi_circ[,c("Gene","circRNA")]
dim(G_circ)
head(G_circ)
G_circ<-unique(G_circ)
dim(G_circ)
head(G_circ)


G_mi_circ_mi <- merge(G_circ, circRNA_miRNA_data, by.x = "circRNA", by.y = "circRNA", all.x = FALSE,all.y=FALSE)
G_mi_circ_mi1<- G_mi_circ_mi
head(G_mi_circ_mi)
dim(G_mi_circ_mi)

G_mi_circ_mi_S<- G_mi_circ_mi[,c("Gene","miRNA")]
dim(G_mi_circ_mi_S)
G_mi_circ_mi_S<-unique(G_mi_circ_mi_S)
dim(G_mi_circ_mi_S)
head(G_mi_circ_mi_S)
write.table(G_mi_circ_mi_S,"DATA/HCC/G_mi_circ_mi_S.txt", quote=F,sep="\t",row.names=FALSE,col.names = TRUE)

head(mRNA_miRNA_data)
#G_mi_circ_mi_T <- read.table("DATA/HCC/G_mi_circ_mi_S.txt", header = TRUE, sep = "\t")
head(G_mi_circ_mi_T)


##Run in another code and get the result
# G_mi_circ_mi_G <- merge(G_mi_circ_mi_T, mRNA_miRNA_data, by.x = "miRNA", by.y = "miRNA", all.x = FALSE,all.y=FALSE)
#library("dplyr")
# G_mi_circ_mi_G<-inner_join(G_mi_circ_mi_T, mRNA_miRNA_data, by =NULL)



dim(G_mi_circ_mi_G)
head(G_mi_circ_mi_G)

lj_b2<-G_mi_circ_mi_G


gmicmig_metapath_matrix <- matrix(0, nrow = length(AllGenesN), ncol = length(AllGenesN), dimnames = list(AllGenesN, AllGenesN))

# Identify edges in G_mi_G and set the corresponding cells in the interaction matrix to 1
edges <- as.matrix(lj_b2[, c("Gene", "Gene2")])
gmicmig_metapath_matrix[edges] <- 1
sum(gmicmig_metapath_matrix>0)

write.table(gmicmig_metapath_matrix,"HCC/Metapaths/MP_HCC_gmicmig_matrix.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)





##### Gene-lncRNA-Gene metapath matrix #####
lncRNA_mRNA_data <- read.table("DATA/HCC/HCC_lncRNA_mRNA.txt", header = FALSE, sep = "\t")
dim(lncRNA_mRNA_data)
head(lncRNA_mRNA_data)
colnames(lncRNA_mRNA_data)<- c("lncRNA","Gene")
mRNA_lncRNA_data<- lncRNA_mRNA_data[,c("Gene","lncRNA")]
head(mRNA_lncRNA_data)

G_lnc_G <- merge(lncRNA_mRNA_data, mRNA_lncRNA_data, by.x = "lncRNA", by.y = "lncRNA", all.x = FALSE,all.y=FALSE)
G_lnc_G1<- G_lnc_G
dim(G_lnc_G)
head(G_lnc_G)

G_lnc_G<-G_lnc_G[,c("Gene.x","Gene.y")]
G_lnc_G<-unique(G_lnc_G)
dim(G_lnc_G)
head(G_lnc_G)

write.table(G_lnc_G,"DATA/HCC/Metapaths/MP_HCC_G_lnc_G.txt", quote=F,sep="\t",row.names=FALSE,col.names = TRUE)

# AllGenes <- read.table("DATA/HCC/HCC_all_genes.txt", header = FALSE, sep = "\t")
# AllGenesN<-AllGenes$V1
# length(AllGenesN)

GlncG_metapath_matrix <- matrix(0, nrow = length(AllGenesN), ncol = length(AllGenesN), dimnames = list(AllGenesN, AllGenesN))

# Identify edges in G_mi_G and set the corresponding cells in the interaction matrix to 1
edgesL <- as.matrix(G_lnc_G[, c("Gene.x", "Gene.y")])
GlncG_metapath_matrix[edgesL] <- 1
sum(GlncG_metapath_matrix>0)

write.table(GlncG_metapath_matrix,"DATA/HCC/Metapaths/MP_HCC_glncg_matrix.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)




##### Gene-miRNA-lncRNA-miRNA-Gene metapath matrix #####
lncRNA_miRNA_data <- read.table("DATA/HCC/HCC_lncRNA_miRNA.txt", header = FALSE, sep = "\t")
dim(lncRNA_miRNA_data)
head(lncRNA_miRNA_data)
colnames(lncRNA_miRNA_data)<- c("lncRNA","miRNA")
head(lncRNA_miRNA_data)
miRNA_lncRNA_data<- lncRNA_miRNA_data[,c("miRNA","lncRNA")]
head(miRNA_lncRNA_data)

G_mi_lnc <- merge(mRNA_miRNA_data, miRNA_lncRNA_data, by.x = "miRNA", by.y = "miRNA", all.x = FALSE,all.y=FALSE)
G_mi_lnc1<- G_mi_lnc

G_lnc <- G_mi_lnc[,c("Gene","lncRNA")]
dim(G_lnc)
head(G_lnc)
G_lnc<-unique(G_lnc)
dim(G_lnc)
head(G_lnc)

dim(unique(rbind(data.frame(g=lncRNA_mRNA_data$V1),data.frame(g=lncRNA_miRNA_data$V1))))


