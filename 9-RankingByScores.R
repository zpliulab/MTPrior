#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")

##Import ranking result
RanksResult <- read.table("DATA/HCC/ranked_indices_All.txt", header = F, sep = "\t")
min(RanksResult)
max(RanksResult)
RanksResult1<-RanksResult+1
min(RanksResult1)
max(RanksResult1)

RanksResultID<-RanksResult1
RanksResultID$Rank<-c(1:nrow(RanksResultID))
colnames(RanksResultID)<- c("IDNumber","rank")
head(RanksResultID)

## Import Disease Genes
HCC_Dgeans_A<- read.table("DATA/HCC/HCC_Diseasegeans_DiGeNet_GWAS.txt", header = F)
head(HCC_Dgeans_A)
dim(HCC_Dgeans_A)


## Import Genes' Names
GeneNames <- read.table("DATA/HCC/filteredGenenNames.txt", header = F, sep = "\t")
GeneNamesID<-GeneNames
GeneNamesID$N<-c(1:nrow(GeneNamesID))
colnames(GeneNamesID)<- c("Genes","IDNumber")
head(GeneNamesID)


## Combine Gene names and Rank result to find ranks of genes
IDrankGenes<-full_join(RanksResultID,GeneNamesID,by="IDNumber")

Genes20<-IDrankGenes[IDrankGenes$rank<=20,]
dim(Genes20[Genes20$Genes %in% HCC_Dgeans_A$V1,])

