#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")
# install.packages("tidyr")
# install.packages("dplyr")
# library(tidyr)
# library(dplyr)

##### Input mRNA Data #####
HCC_mRNA<- read.delim("DATA/HCC/LIHC_mRNA.txt", header = TRUE, sep = "\t")
# length(unique(HCC_mRNA$X))
#### ### ## #

HCC_mRNA1<-HCC_mRNA

# Gene_Ens<- data.frame(unique(rownames(HCC_mRNA)))
# write.table(Gene_Ens,"DATA/HCC/geneENS.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

##### Input GeneEnsambl and GeneSymbol and Gene ID #####
G_ENS_ID<- read.delim("DATA/HCC/GeneEns-GID.txt", header = TRUE, sep = "\t")
G_ENS_ID_S<-G_ENS_ID[,c("query","entrezgene","symbol")]
colnames(G_ENS_ID_S)<- c("Ensemblgene","geneID","symbol")
head(G_ENS_ID_S)
###  ************** #HERE


##### Input Gene ID and Gene name Data 
gene_id_name<- read.delim("DATA/HCC/gene_id_name.txt", header = TRUE, sep = "\t")
gene_id_name<-gene_id_name[,2:3]
colnames(gene_id_name)<- c("GeneName","GeneID")


##### Replace row names in HCC_mRNA with Gene ID
HCC_mRNA1<- data.frame(GeneN=rownames(HCC_mRNA1),HCC_mRNA1)

HCC_mRNA1 <- merge(HCC_mRNA1, gene_id_name, by.x = "GeneN", by.y = "GeneName", all.x = FALSE)
HCC_mRNA1$GeneN <- HCC_mRNA1$GeneID

HCC_mRNA1<-HCC_mRNA1[,!(names(HCC_mRNA1) %in% c("rownames.HCC_mRNA1."))]
HCC_mRNA1 <- HCC_mRNA1[, !colnames(HCC_mRNA1) %in% "GeneID"]

# install.packages("dplyr")
# library(dplyr)
# HCC_mRNA2 <- HCC_mRNA1 %>% group_by(GeneN) %>% summarise_all(mean)
HCC_mRNA2 <- aggregate(. ~ GeneN, data = HCC_mRNA1, FUN = mean)
HCC_mRNA1<-HCC_mRNA2

#write.table(HCC_mRNA1,"DATA/HCC/HCC_mRNA_genename.txt", quote=F,sep="\t",row.names=FALSE,col.names = TRUE)


##### Input GGI Data 
GGI<- read.delim("DATA/9606.protein.links.v12.0.txt", header = TRUE, sep = " ")

###Contains Protein ID and Gene names
proteinInfo<- read.delim("DATA/9606.protein.info.v12.0.txt", header = TRUE, sep = "\t")
###Keep just Protein ID and Gene names
proteinInfo<-proteinInfo[,1:2]
###Save a backup of GGI Data
GGI_BackUP<-GGI

#### ### ## #
# proteinEns <- data.frame(sub("^9606\\.", "", proteinInfo$X.string_protein_id))
# proteinEns<-data.frame(P=unique(proteinEns[,1]))
# write.table(proteinEns,"DATA/HCC/proteinEns.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

P_ENS_ID<- read.delim("DATA/HCC/ProteinENS-GID.txt", header = TRUE, sep = "\t")
###  ************** #HERE


##### Replace Protein ID with Gene Names in GGI Dataframe  #####
###Merge for protein1 column
GGI <- merge(GGI, proteinInfo, by.x = "protein1", by.y = "X.string_protein_id", all.x = FALSE)
GGI$protein1 <- GGI$preferred_name
GGI <- GGI[, !(names(GGI) %in% c("X.string_protein_id", "preferred_name"))]
###Merge for protein2 column
GGI <- merge(GGI, proteinInfo, by.x = "protein2", by.y = "X.string_protein_id", all.x = FALSE)
GGI$protein2 <- GGI$preferred_name
GGI <- GGI[, !(names(GGI) %in% c("X.string_protein_id", "preferred_name"))]
GGI<-GGI[,1:2]
colnames(GGI)<-c("Gene1","Gene2")
dim(GGI)
head(GGI)

#### ### ## #


##### Extract and Save Gene node list and gene-gene edge list #####
###gene node list
###Extract gene names from HCC_mRNA
genes_HCC_mRNA <- unique(HCC_mRNA1$GeneN)
genes_HCC_mRNA_df<-data.frame(genes=genes_HCC_mRNA)
dim(genes_HCC_mRNA_df)
##### GENE ORDER #####
Gene_order<- genes_HCC_mRNA

###Extract gene names from both columns in GGI
genes_GGI <- unique(c(GGI$Gene1, GGI$Gene2))

all_genes_df<-genes_HCC_mRNA_df
all_genes <- all_genes_df$genes
dim(all_genes_df)

filtered_GGI <- GGI[GGI$Gene1 %in% all_genes,]
filtered_GGI <- filtered_GGI[filtered_GGI$Gene2 %in% all_genes,]
dim(filtered_GGI)

###Save gene node list
#write.table(genes_HCC_mRNA_df,"DATA/HCC/HCC_all_genes.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

###Save Gene-Gene Edges List
#write.table(filtered_GGI,"DATA/HCC/HCC_filtered_GGI.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

#### ### ## #



##### Create GGI matrix contains all genes #####
# Create an empty interaction matrix
interaction_matrix <- matrix(0, nrow = length(all_genes), ncol = length(all_genes), dimnames = list(all_genes, all_genes))

# Identify edges in filtered_GGI and set the corresponding cells in the interaction matrix to 1
edges <- as.matrix(filtered_GGI[, c("Gene1", "Gene2")])
interaction_matrix[edges] <- 1
sum(interaction_matrix>0)

#write.table(interaction_matrix,"DATA/HCC/HCC_interaction_matrix.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

#### ### ## #



##### Input miRNA data from TCGA #####
# HCC_miRNA<- read.delim("DATA/HCC/LIHC_miRNA.txt", header = TRUE, sep = "\t")


##### load miRNA_mRNA interactions ########
####       2 dataset    
GeneID_miRNA<- read.delim("DATA/GeneID_miRNA.txt", header = TRUE, sep = "\t")
GeneID_miRNA <- GeneID_miRNA[grepl("^hsa", GeneID_miRNA$miRNA), ]
GeneID_miRNA<-GeneID_miRNA[,c("miRNA","Gene.Symbol")]
colnames(GeneID_miRNA)<- c("miRNA", "Gene")
dim(GeneID_miRNA)
head(GeneID_miRNA)


miRNA_mRNA<- read.delim("DATA/miRTarBase_MTI.txt", header = TRUE, sep = "\t")
miRNA_mRNA <- miRNA_mRNA[grepl("^hsa", miRNA_mRNA$miRNA), ]
miRNA_mRNA<-miRNA_mRNA[,c("miRNA","Target.Gene")]
colnames(miRNA_mRNA)<- c("miRNA","Gene")
dim(miRNA_mRNA)
head(miRNA_mRNA)


##### Save edges and nodes list of miRNA-mRNA (All data) #####

###Find all edges of miRNA-mRNA
###Merge the dataframes using a full outer join
miRNA_mRNA_df<- rbind(GeneID_miRNA, miRNA_mRNA)

###Remove duplicate rows
unique_miRNA_mRNA <- unique(miRNA_mRNA_df)

miRNA_mRNA_All<-unique_miRNA_mRNA[,c("miRNA", "Gene")]
dim(miRNA_mRNA_All)
head(miRNA_mRNA_All)



all_miRNA <- unique(c(miRNA_mRNA_All$miRNA))
all_miRNA_df<-data.frame(all_miRNA)


###save edge list of miRNA_mRNA
#write.table(miRNA_mRNA_All,"DATA/HCC/miRNA_mRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###save node list of miRNA
#write.table(all_miRNA_df,"DATA/HCC/all_miRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

#### ### ## #


##### Filter miRNA-mRNA edge list according to previous gene list #####

miRNA_mRNA_Gfilter<- miRNA_mRNA_All[miRNA_mRNA_All$Gene %in% all_genes,]
dim(miRNA_mRNA_Gfilter)

#write.table(miRNA_mRNA_Gfilter,"DATA/HCC/miRNA_mRNA_Gfilter.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

miRNA_Gfilter<- unique(c(miRNA_mRNA_Gfilter$miRNA))
miRNA_Gfilter_df<- data.frame(miRNA_Gfilter)

#write.table(miRNA_Gfilter_df,"DATA/HCC/miRNA_Gfilter.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)



##DELETE LAST PART OF MIRNA
# miRNA_Gfilter_df_C<-miRNA_Gfilter_df
# miRNA_Gfilter_df_C$miRNA_Gfilter <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", miRNA_Gfilter_df_C$miRNA_Gfilter)

miRNA_mRNA_Gfilter_C<- miRNA_mRNA_Gfilter
miRNA_mRNA_Gfilter_C$miRNA<- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", miRNA_mRNA_Gfilter_C$miRNA)

## low case
# miRNA_Gfilter_df_C_l<-miRNA_Gfilter_df_C
# miRNA_Gfilter_df_C_l$miRNA_Gfilter<- tolower(miRNA_Gfilter_df_C_l$miRNA_Gfilter)

miRNA_mRNA_Gfilter_C_l<-miRNA_mRNA_Gfilter_C
miRNA_mRNA_Gfilter_C_l$miRNA<- tolower(miRNA_mRNA_Gfilter_C_l$miRNA)

miRNA_mRNA_Gfilter_C_l<-unique(miRNA_mRNA_Gfilter_C_l)

miRNA_Gfilter_df_C_0<-data.frame(unique(miRNA_mRNA_Gfilter_C_l$miRNA))

#         
#write.table(miRNA_mRNA_Gfilter_C_l,"DATA/HCC/miRNA_mRNA_Gfilter_Clean.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

#write.table(miRNA_Gfilter_df_C_0,"DATA/HCC/miRNA_Gfilter_Clean.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)




##### load miRNA-Disease Data #####

###The whole dataset of miRNA-disease association data 
D_miRNA<- read.delim("DATA/HMDD_miRNA_Disease.txt", header = TRUE, sep = "\t")
# dim(miRNA_Gfilter_df_C[miRNA_Gfilter_df_C$miRNA_Gfilter %in% D_miRNA$miRNA,])
# dim(D_miRNA[D_miRNA$miRNA %in% miRNA_Gfilter_df_C$miRNA_Gfilter,])


D_Names<-data.frame(unique(D_miRNA$disease))
specific_word <- "Hepatocellular"
# Find rows where the first column contains the specific word
rows_with_specific_word <- which(grepl(specific_word, D_Names[,1], ignore.case = TRUE))
# Extract the rows from the dataframe that shows all expression contains the specific word
result_dataframe <- D_Names[rows_with_specific_word, ]
### select just miRNA are related to the specific disease
miRNA_Specific<-D_miRNA[D_miRNA$disease %in% result_dataframe,]
dim(miRNA_Specific)

miRNA_Specific<- data.frame(unique(miRNA_Specific$miRNA))
colnames(miRNA_Specific)<-c("miRNA")
dim(miRNA_Specific)


### %%%%%%%%%%%%% mir miR Chat 

miRNA_Specific$miRNA <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", miRNA_Specific$miRNA)
dim(miRNA_Gfilter_df_C_0[miRNA_Gfilter_df_C_0$unique.miRNA_mRNA_Gfilter_C_l.miRNA. %in% miRNA_Specific$miRNA ,])
dim(miRNA_Specific[miRNA_Specific$miRNA  %in% miRNA_Gfilter_df_C_0$unique.miRNA_mRNA_Gfilter_C_l.miRNA.,])
### %%%%%%%%%%%%%

###Select miRNA-mRNA contain miRNA_Specific
# miRNA_mRNA_HCC<-miRNA_mRNA_Gfilter_C_l[miRNA_mRNA_Gfilter_C_l$miRNA %in% miRNA_Specific$miRNA,]  ###Zero!
# dim(miRNA_mRNA_HCC)

miRNA_mRNA_Gfilter_C_l[miRNA_mRNA_Gfilter_C_l$miRNA %in% miRNA_Specific$miRNA,]
miRNA_Specific[miRNA_Specific$miRNA %in% miRNA_mRNA_Gfilter_C_l$miRNA,]
## *****************
##????
# dtest<-data.frame(unique(miRNA_mRNA_Gfilter_C_l$miRNA))
# dtest[dtest$unique.miRNA_mRNA_Gfilter_C_l.miRNA. %in% miRNA_Specific$miRNA,]


HCC_miRNA_mRNA1<- miRNA_mRNA_Gfilter_C_l[miRNA_mRNA_Gfilter_C_l$miRNA %in% miRNA_Specific$miRNA,]
## ****** ******

### mir2disease data of miRNA-disease association
mir2disease_miRNA<-read.delim("DATA/MicroRNA-target1.txt", header = TRUE, sep = "\t")


D_Names1<-data.frame(unique(mir2disease_miRNA$Reference))
specific_word1 <- "hepatocellular"
# Find rows where the first column contains the specific word
rows_with_specific_word1 <- which(grepl(specific_word1, D_Names1[,1], ignore.case = TRUE))
# Extract the rows from the dataframe that shows all expression contains the specific word
result_dataframe1 <- D_Names1[rows_with_specific_word1, ]
result_dataframe1
### select just miRNA are related to the specific disease
miRNA_Specific1<-mir2disease_miRNA[mir2disease_miRNA$Reference %in% result_dataframe1,]
dim(miRNA_Specific1)


miRNA_Specific1$miRNA <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", miRNA_Specific1$miRNA)
miRNA_Specific1$miRNA<-tolower(miRNA_Specific1$miRNA)
miRNA_mRNA_Gfilter_C_l[miRNA_mRNA_Gfilter_C_l$miRNA %in% miRNA_Specific1$miRNA,]
miRNA_Specific1[miRNA_Specific1$miRNA%in% miRNA_mRNA_Gfilter_C_l$miRNA,]

dim(miRNA_mRNA_Gfilter_C_l)
dim(miRNA_mRNA_Gfilter_C_l[miRNA_mRNA_Gfilter_C_l$miRNA %in% miRNA_Specific1$miRNA,])
## *****************************


HCC_miRNA_mRNA2<- miRNA_mRNA_Gfilter_C_l[miRNA_mRNA_Gfilter_C_l$miRNA %in% miRNA_Specific1$miRNA,]
## ****** ******


#####Save edges and nodes list of miRNA-mRNA (Related to disease data)
HCC_miRNA_mRNA<- rbind(HCC_miRNA_mRNA1,HCC_miRNA_mRNA2)
HCC_miRNA_mRNA<- unique(HCC_miRNA_mRNA)

HCC_miRNA<-(data.frame(unique(HCC_miRNA_mRNA$miRNA)))

#write.table(HCC_miRNA_mRNA,"DATA/HCC/HCC_miRNA_mRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
#write.table(HCC_miRNA,"DATA/HCC/HCC_miRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)





##### load circRNA-miRNA interactions ########
circRNA_miRNA<- read.delim("DATA/circRNA_miRNA interaction.csv",sep=",")


### keep rows that their miRNA value is in all_miRNA list
circRNA_miRNA_C<-circRNA_miRNA
circRNA_miRNA_C$mir1<- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", circRNA_miRNA$mir1)
circRNA_miRNA_C$mir1<- tolower(circRNA_miRNA_C$mir1)

circRNA_miRNA_Clean<- circRNA_miRNA_C[circRNA_miRNA_C$mir1 %in% miRNA_Gfilter_df_C_0[,1],]
dim(circRNA_miRNA_Clean)
circRNA_miRNA_Clean<-unique(circRNA_miRNA_Clean)

### keep circRNA and miRNA columns
circRNA_miRNA_Clean<-circRNA_miRNA_Clean[,1:2]

### circRNA nodes
all_circRNA <- unique(c(circRNA_miRNA_Clean$CircID))
all_circRNA_df<-data.frame(all_circRNA)

###save edge list of circRNA-miRNA
#write.table(circRNA_miRNA_Clean,"DATA/HCC/circRNA_miRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###save node list of circRNA
#write.table(all_circRNA_df,"DATA/HCC/all_circRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

#### ### ## #


#####    load circRNA-Disease      ##### 
D_circRNA<- read.delim("DATA/HCC/LIHC_circRNA.csv", header = TRUE, sep = ",")
#dim(circRNA_miRNA_Clean[circRNA_miRNA_Clean$CircID %in% D_circRNA$CircID,])
circRNA_miRNA_LIHC<-circRNA_miRNA_Clean[circRNA_miRNA_Clean$CircID %in% D_circRNA$CircID,]
LIHC_circRNA<-(unique(circRNA_miRNA_LIHC$CircID))
LIHC_circRNA_df<-data.frame(LIHC_circRNA)
LIHC_circRNA

###save edge list of circRNA-miRNA for HCC
#write.table(circRNA_miRNA_LIHC,"DATA/HCC/HCC_circRNA_miRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###save node list of circRNA for HCC
#write.table(LIHC_circRNA_df,"DATA/HCC/HCC_circRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)

#### ### ## #
##HERE
circRNA_miRNA_LIHC_f<- (circRNA_miRNA_LIHC[circRNA_miRNA_LIHC$CircID %in% LIHC_circRNA,])
circRNA_miRNA_LIHC_f2<- (circRNA_miRNA_LIHC_f[circRNA_miRNA_LIHC_f$mir1 %in% miRNA_Specific$miRNA,])
dim(unique(circRNA_miRNA_LIHC_f2))
length(unique(circRNA_miRNA_LIHC_f2$CircID))
circRNA_miRNA_LIHC_f3<- circRNA_miRNA_LIHC_f2[circRNA_miRNA_LIHC_f2$mir1 %in% lncRNA_miRNA_HCC$miRNA,]
length(unique(circRNA_miRNA_LIHC_f3$CircID))


#####      load lncRNA-mRNA interactions      ########
lncRNA_mRNA<- read.delim("DATA/lncRNA_target.txt", header = TRUE, sep = "\t")
head(lncRNA_mRNA)
dim(lncRNA_mRNA)

lncRNA_mRNA_c<-lncRNA_mRNA[lncRNA_mRNA$Species %in% c("9606"),]
dim(lncRNA_mRNA_c)

##
lncRNA_name_EnsembleID<-lncRNA_mRNA_c[,c("LncRNA_official_symbol","Ensembl_ID")]
#write.table(lncRNA_name_EnsembleID,"DATA/HCC/lncRNA_name_EnsembleID.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
##

lncRNA_mRNA_c<-lncRNA_mRNA_c[,c("LncRNA_official_symbol","Target_official_symbol")]
lncRNA_mRNA_c<-unique(lncRNA_mRNA_c)
colnames(lncRNA_mRNA_c)<- c("LncRNA","genes")
head(lncRNA_mRNA_c)
dim(lncRNA_mRNA_c)


lncRNA_mRNA_c<-lncRNA_mRNA_c[lncRNA_mRNA_c$genes %in% all_genes,]
dim(lncRNA_mRNA_c)

all_lncRNA <- unique(c(lncRNA_mRNA_c$LncRNA))
all_lncRNA_df<-data.frame(all_lncRNA)
dim(all_lncRNA_df)

###save edge list of lncRNA-mRNA for HCC
#write.table(lncRNA_mRNA_c,"DATA/HCC/lncRNA_mRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###save node list of circRNA for HCC
#write.table(all_lncRNA_df,"DATA/HCC/lncRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)





###### load lncRNA-Disease  #####
D_lncRNA<- read.delim("DATA/lncDiseasev3.tsv", header = TRUE, sep = "\t")
dim(D_lncRNA)
head(D_lncRNA)
###What category of data does it have?
unique(D_lncRNA$ncRNA.Category)  ## "LncRNA"  "CircRNA" "miRNA" ## 

### Just keep Homo sapiens data 
D_lncRNA_c<- D_lncRNA[D_lncRNA$Species %in% c("Homo sapiens"),]
dim(D_lncRNA_c)



### separate categories
lnc_D_lncRNA_c<- D_lncRNA_c[D_lncRNA_c$ncRNA.Category %in% c("LncRNA"),]
dim(lnc_D_lncRNA_c)
circ_D_lncRNA_c<- D_lncRNA_c[D_lncRNA_c$ncRNA.Category %in% c("CircRNA"),]
dim(circ_D_lncRNA_c)
miRNA_D_lncRNA_c<- D_lncRNA_c[D_lncRNA_c$ncRNA.Category %in% c("miRNA"),]
dim(miRNA_D_lncRNA_c)


### find rows contain specific disease
specific_word2 <- "hepatocellular"
# Find rows where the first column contains the specific word
rows_with_specific_word2 <- which(grepl(specific_word2, lnc_D_lncRNA_c$Disease.Name, ignore.case = TRUE))
# Extract the rows from the dataframe that shows all expression contains the specific word
result_dataframe2 <- lnc_D_lncRNA_c[rows_with_specific_word2, ]
head(result_dataframe2)
dim(result_dataframe2)



### select just lncRNA are related to the specific disease
D_lncRNA_HCC <- result_dataframe2
HCC_lncRNA<-unique(D_lncRNA_HCC$ncRNA.Symbol)
HCC_lncRNA_df<-data.frame(HCC_lncRNA)


### lncRNA-mRNA interaction related to HCC
lncRNA_mRNA_HCC<- lncRNA_mRNA_c[lncRNA_mRNA_c$LncRNA %in% HCC_lncRNA_df$HCC_lncRNA,]
lncRNA_mRNA_HCC<- unique(lncRNA_mRNA_HCC)



###save edge list of lncRNA-mRNA for HCC
#write.table(lncRNA_mRNA_HCC,"DATA/HCC/HCC_lncRNA_mRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###save node list of lncRNA for LGG
#write.table(HCC_lncRNA_df,"DATA/HCC/HCC_lncRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)





#####      load lncRNA-miRNA interactions     ########
lncRNA_miRNA<- read.delim("DATA/lncRNASNP_lncRNA_miRNA_interactions.txt", header = TRUE, sep = "\t")
dim(lncRNA_miRNA)
head(lncRNA_miRNA)

lncRNA_miRNA_c<- lncRNA_miRNA[,c("lncRNA","miRNA")]

#here
lncRNA_miRNA_c$miRNA<- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", lncRNA_miRNA_c$miRNA)
lncRNA_miRNA_c$miRNA<- tolower(lncRNA_miRNA_c$miRNA)


lncRNA_miRNA_miFiltered<- lncRNA_miRNA_c[lncRNA_miRNA_c$miRNA %in% miRNA_Gfilter_df_C_0[,1],]
dim(lncRNA_miRNA_miFiltered)
lncRNA_miRNA_miFiltered<-unique(lncRNA_miRNA_miFiltered)
dim(lncRNA_miRNA_miFiltered)

all_lncRNA2 <- unique(c(lncRNA_miRNA_miFiltered$lncRNA))
all_lncRNA2_df<-data.frame(all_lncRNA2)
dim(all_lncRNA2_df)

###save edge list of lncRNA-miRNA
#write.table(lncRNA_miRNA_miFiltered,"DATA/HCC/lncRNA_miRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###save node list of lncRNA
#write.table(all_lncRNA2_df,"DATA/HCC/lncRNA2mi.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)



#####      lncRNA interaction related to HCC     ##### 
D_lncRNA2<- read.delim("DATA/lncRNASNP_lncRNA_diseases.txt", header = TRUE, sep = "\t")
D_lncRNA2[1:10,]



### find rows contain specific disease
D_lncRNA2_HCC <- D_lncRNA2[grep(specific_word, D_lncRNA2$disease,ignore.case = TRUE), ]

### find rows contain specific disease
specific_word2 <- "hepatocellular"
# Find rows where the first column contains the specific word
rows_with_specific_word3 <- which(grepl(specific_word2, D_lncRNA2_HCC$disease, ignore.case = TRUE))
# Extract the rows from the dataframe that shows all expression contains the specific word
result_dataframe3 <- D_lncRNA2_HCC[rows_with_specific_word3, ]
head(result_dataframe3)
dim(result_dataframe3)




### lncRNA-miRNA interaction related to HCC
lncRNA_miRNA_HCC<-lncRNA_miRNA_miFiltered[lncRNA_miRNA_miFiltered$lncRNA %in% result_dataframe3$transcript,]
dim(lncRNA_miRNA_HCC)



HCC_lncRNA<-unique(lncRNA_miRNA_HCC$lncRNA)
HCC_lncRNA_df<-data.frame(HCC_lncRNA)
dim(HCC_lncRNA_df)
###save edge list of lncRNA-miRNA for LGG
#write.table(lncRNA_miRNA_HCC,"DATA/HCC/HCC_lncRNA_miRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
###save node list of lncRNA2 for LGG
#write.table(HCC_lncRNA_df,"DATA/HCC/HCC_lncRNA.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)



