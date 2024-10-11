#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")

##Import ranking result
ProbabilityResult <- read.table("DATA/HCC/probabilities_All.txt", header = F, sep = "\t")
dim(ProbabilityResult)
head(ProbabilityResult)


ProbabilityResultID<-ProbabilityResult
ProbabilityResultID$IDNumber<-c(1:nrow(ProbabilityResult))
colnames(ProbabilityResultID)<- c("prob","IDNumber")
head(ProbabilityResultID)
ProbabilityResultID<-ProbabilityResultID[,c("IDNumber","prob")]


## Import Genes' Names
GeneNames <- read.table("DATA/HCC/GeneNamesF.txt", header = F, sep = "\t")
GeneNamesID<-GeneNames
# GeneNamesID$N<-c(1:nrow(GeneNamesID))
# colnames(GeneNamesID)<- c("Genes","IDNumber")
# head(GeneNamesID)
# 




## Import Disease Genes
HCC_Dgeans_A<- read.table("DATA/HCC/HCC_Diseasegeans_DiGeNet_GWAS.txt", header = F)
head(HCC_Dgeans_A)
dim(HCC_Dgeans_A)



probIDrank<- data.frame(IDNumber=ProbabilityResultID[,1],GeneName=GeneNamesID[,1],prob=ProbabilityResultID[,2])
head(probIDrank)


probIDrank$Type<- 0
probIDrank[probIDrank$GeneName %in% HCC_Dgeans_A$V1, ]$Type<- 1
#sum(probIDrank$Type)


library(dplyr)
df_sorted_desc <- probIDrank %>%
  arrange(desc(prob))

probIDrank_A <- df_sorted_desc
probIDrank_A$rank<- c(1:nrow(probIDrank_A))
head(probIDrank_A)





#####lnc prob Data
DlncProb_A <-read.table("Data/DlncProb_n.txt", header=T)
dim(DlncProb_A)
head(DlncProb_A)

DlncProb_Or<-DlncProb_A %>%
  arrange(desc(prob))

DlncProb_Or$rank<- c(1:nrow(DlncProb_Or))





library("ggplot2")


DisGeneProb<- probIDrank_A[probIDrank_A$Type==1,]
dim(DisGeneProb)

violindata1<- data.frame(prob=DisGeneProb[,3],C="Genes")


DlncProb<-DlncProb_A[DlncProb_A$Type == 1,]

violindata2<- data.frame(prob=DlncProb[,2],C="LncRNA")

violindata<-rbind(violindata1,violindata2)

custom_colors=c("#1E90FF", "#FFD700")
pviol<- ggplot(violindata, aes(x = C, y = prob, alpha = 0.5)) +
  geom_violin(aes(fill = C),scale="width", trim = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = NULL,
       x = NULL,
       y = "Value") +
  scale_fill_manual(values = custom_colors) +
  theme_bw()+
  theme(legend.position = "none")


pviol

#pdf("Code/Plots/Distribution3.pdf",width = 4,height =4)
pviol
#dev.off()








# How many in top
library(dplyr)


# probIDrank_A
# DlncProb_Or



# IDrankGenesADD<-IDrankGenes
# IDrankGenesADD$DN<-0
# IDrankGenesADD[IDrankGenesADD$Genes %in% HCC_Dgeans_A$V1,]$DN<-1

# df <- IDrankGenesADD[order(IDrankGenesADD$rank), ]
# rank_values <- unique(df$rank)
# # Initialize an empty dataframe to store the result
# result_df <- data.frame(rank = unique(df$rank), disease_gene_count = NA)


C_gene<- probIDrank_A[,c("prob","Type","rank")]
C_lnc<- DlncProb_Or[,c("prob","Type","rank")]
C_gene$C<-"Gene"
C_lnc$C<-"LncRNA"

gene_count<-C_gene[C_gene$Type==1,]
lnc_count<-C_lnc[C_lnc$Type==1,]

gene_count200<-gene_count[gene_count$rank<=200,]
lnc_count200<-lnc_count[lnc_count$rank<=200,]

gene_count200$count<- NA
lnc_count200$count<- NA

# Loop through each rank value
for (i in 1:nrow(gene_count200)) {
  # Count the disease genes with rank less than or equal to the current rank value
  gene_count200$count[i] <- sum(gene_count200$rank <= gene_count200$rank[i] )
}

# Loop through each rank value
for (i in 1:nrow(lnc_count200)) {
  # Count the disease LNC with rank less than or equal to the current rank value
  lnc_count200$count[i] <- sum(lnc_count200$rank <= lnc_count200$rank[i] )
}



gene_count200_C<-gene_count200[,c("rank","count","C")]
lnc_count200_C<-lnc_count200[,c("rank","count","C")]



addlnc<-data.frame(rank=199,count=111,C="LncRNA")

lnc_count200_C_ADD<-rbind(lnc_count200_C,addlnc)



coun200<-rbind(gene_count200_C,lnc_count200_C_ADD)
coun200$C<-factor(coun200$C)

custom_colors=c("#1E90FF", "#FFD700")

# Create the line plot
countTop_Plot<- ggplot(coun200, aes(x = rank, y =count, group = C )) +
  geom_line(aes(color = C),linewidth=1) + 
  geom_area(aes(fill = C), alpha = 0.07, position = "identity") +
  labs(x = "Rank", y = "Count", title = NULL) +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors)+
  theme_bw()+
  theme(legend.position = c(0.86, 0.15))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 201))+ #+  # Ensure the x-axis starts at 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200))

countTop_Plot

# pdf("Code/Plots/countTop1.pdf",width = 4,height =4)
# countTop_Plot
# dev.off()




#####################################################
# install.packages("VennDiagram")
# library(VennDiagram)

## Define the sets
set1 <- probIDrank_A[probIDrank_A$rank<=100,]$GeneName

#set1 <- probIDrank_A$GeneName
set2 <- HCC_Dgeans_A$V1

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(Set1 = set1, Set2 = set2),
  category.names = c("",""),  # No titles for the circles
  fill = c("#1E90FF", "#FFD700"),  # Blue and yellow
  alpha = 0.5,
  cat.pos = 0,
  cat.dist = 0.05,
  cex = 1.2,  # Smaller font size for the inside text
  cat.cex = 1.5,
  lwd = 0,  # No border lines
  filename = NULL
)

dev.off()
grid.draw(venn.plot)


# Plot the Venn diagram
# pdf("Code/Plots/venn_G_top100_2.pdf",width = 2,height =2)
# grid.draw(venn.plot)
# dev.off()


### lncRNA venn diagram

addlnc<-data.frame(rank=199,count=111,C="LncRNA")

lnc_count200_C_ADD<-rbind(lnc_count200_C,addlnc)


add_DlncProb_Or<- DlncProb_Or
add_DlncProb_Or$I<- c(1:nrow(DlncProb_Or))

extralnc<- data.frame(Type=0, rank=c(169:200),I=c(169:200))

extralnc_C<-rbind(add_DlncProb_Or[,c("Type","rank","I")],extralnc)

## Define the sets
#top ranked
set1 <- extralnc_C[extralnc_C$rank<=2100, ]$I

#known disease lncRNA
set2 <- extralnc_C[extralnc_C$Type==1,]$I

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(Set1 = set1, Set2 = set2),
  category.names = c("",""),  # No titles for the circles
  fill = c("#1E90FF", "#FFD700"),  # Blue and yellow
  alpha = 0.5,
  cat.pos = 0,
  cat.dist = 0.05,
  cex = 1.2,  # Smaller font size for the inside text
  cat.cex = 1.5,
  lwd = 0,  # No border lines
  filename = NULL
)

dev.off()
grid.draw(venn.plot)


# Plot the Venn diagram
# pdf("Code/Plots/venn_L_top200_2.pdf",width = 2,height =2)
# grid.draw(venn.plot)
# dev.off()


