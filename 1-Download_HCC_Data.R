####https://portal.gdc.cancer.gov/repository
setwd("F:/Code/HCC/")

library(TCGAbiolinks)
library(SummarizedExperiment)


##### Download mRNA DATA #####

query_mRNA <- GDCquery(project = "TCGA-LIHC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts"
)

GDCdownload(query_mRNA )

LIHC_mRNA <- GDCprepare(query_mRNA , save = T,save.filename = "DATA/LIHC_mRNA.Rdata")

mRNA_dataPrep <- TCGAanalyze_Preprocessing(object = LIHC_mRNA,
                                           cor.cut = 0.6,
                                           datatype = "stranded_first")

write.table(mRNA_dataPrep,"DATA/LIHC_mRNA.txt", quote=F,sep="\t",row.names=TRUE,col.names = TRUE)

mRNA_LIHC<-read.delim("DATA/LIHC_mRNA.txt", header = TRUE, sep = "\t")


#normalized_counts <- assay(LIHC_mRNA)


###### Input the data-set contains ENSG ID and RNA Type #####
ENSG_RNA<-read.delim("DATA/mart_export.txt", header = TRUE, sep = "\t")
dim(ENSG_RNA)
head(ENSG_RNA)


ENSG_RNA_S<-ENSG_RNA[,c("Gene.stable.ID.version","Gene.type")]
dim(ENSG_RNA_S)
head(ENSG_RNA_S)
colnames(ENSG_RNA_S)<- c("ENSG","Gene.type")

###### ENSG ID of feature matrix and RNA Type
HCC_ENS<- data.frame(rownames(mRNA_LIHC))
dim(HCC_ENS)
head(HCC_ENS)
colnames(HCC_ENS)<- c("ENSG")

ENSGType<-merge(x = HCC_ENS, y = ENSG_RNA_S, by = "ENSG")
dim(ENSGType)
head(ENSGType)

unique(ENSGType$Gene.type)

lncRNA_ENSG<-ENSGType[ENSGType$Gene.type %in% c("lncRNA"),]
dim(data.frame(lncRNA_ENSG))
head(lncRNA_ENSG)


mRNA_LIHC_namesAdd<- data.frame(ENSG=rownames(mRNA_LIHC),mRNA_LIHC)
dim(mRNA_LIHC_namesAdd)


mRNA_LIHC_namesTypeAdd<- merge(x=mRNA_LIHC_namesAdd, y=lncRNA_ENSG, by="ENSG")
dim(mRNA_LIHC_namesTypeAdd)
mRNA_LIHC_namesTypeAdd[1,]
mRNA_LIHC_namesTypeAdd_C<- mRNA_LIHC_namesTypeAdd[,-ncol(mRNA_LIHC_namesTypeAdd)]

dim(mRNA_LIHC_namesTypeAdd_C)

write.table(mRNA_LIHC_namesTypeAdd_C,"HCC/DATA/HCC_lncRNAfeature_TCGA.txt", quote=F,sep="\t",row.names=T,col.names = TRUE)







