#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")


######      load lncRNA-Disease      #####
#lncDiseasev3
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
dim(HCC_lncRNA_df)
head(HCC_lncRNA_df)

##################
HCC_lncRNA_df1<-HCC_lncRNA_df
#################



##lncRNASNP2
#####      lncRNA interaction related to HCC     ##### 
D_lncRNA2<- read.delim("DATA/lncRNASNP_lncRNA_diseases.txt", header = TRUE, sep = "\t")
D_lncRNA2[1:10,]

dim(D_lncRNA2)
length(unique(D_lncRNA2$transcript))
length(unique(D_lncRNA2$disease))


### find rows contain specific disease
D_lncRNA2_HCC <- D_lncRNA2[grep(specific_word2, D_lncRNA2$disease,ignore.case = TRUE), ]

### find rows contain specific disease
specific_word2 <- "hepatocellular"
# Find rows where the first column contains the specific word
rows_with_specific_word3 <- which(grepl(specific_word2, D_lncRNA2_HCC$disease, ignore.case = TRUE))
# Extract the rows from the dataframe that shows all expression contains the specific word
result_dataframe3 <- D_lncRNA2_HCC[rows_with_specific_word3, ]
head(result_dataframe3)
dim(result_dataframe3)

#######################
HCC_lncRNA_df2 <- data.frame(result_dataframe3$transcript)
#######################



###### lncRNA ENSG and ENST ID dataset ######
mart_data<-read.delim("DATA/HCC/lncRNAData/mart_export.txt", header = T, sep = "\t")
head(mart_data)
unique(mart_data$Gene.type)
dim(mart_data)
lnc_mart_data<- mart_data[mart_data$Gene.type %in% c("lncRNA"),]
dim(lnc_mart_data)
lnc_mart_data_C<- lnc_mart_data[,c("Gene.stable.ID","Transcript.stable.ID")]
dim(lnc_mart_data_C)
head(lnc_mart_data_C)

#write.table(lnc_mart_data_C,"DATA/HCC/lncRNAData/lncConvert1.txt", quote=F,sep="\t",row.names=FALSE,col.names = TRUE)




###### lncRNA name and ID dataset 1######
##     http://dibresources.jcbose.ac.in/zhumur/lncrbase2/start2.php   ##
lncRNA_ENSG_ID<-read.delim("DATA/HCC/lncRNAData/human_alias_info.txt", header = F, sep = "\t")
head(lncRNA_ENSG_ID)
dim(lncRNA_ENSG_ID)

lncRNA_ENSG_ID123<- lncRNA_ENSG_ID [,c(1,2,3,6,7)]
colnames(lncRNA_ENSG_ID123)<- c("N1","N2","N3","N6","N7")
lncRNA_ENSG_ID124<- lncRNA_ENSG_ID [,c(1,2,4,6,7)]
colnames(lncRNA_ENSG_ID124)<- c("N1","N2","N3","N6","N7")
lncRNA_ENSG_ID_12367<- rbind(lncRNA_ENSG_ID123,lncRNA_ENSG_ID124)
dim(lncRNA_ENSG_ID_12367)
lncRNA_ENSG_ID_12367<-unique(lncRNA_ENSG_ID_12367)

lncRNA_ENSG_ID_m1<- lncRNA_ENSG_ID_12367[,c(1,2,3,4)]
lncRNA_ENSG_ID_m2<- lncRNA_ENSG_ID_12367[,c(1,2,3,5)]
colnames(lncRNA_ENSG_ID_m1)<- c("N1","N2","N3","N4")
colnames(lncRNA_ENSG_ID_m2)<- c("N1","N2","N3","N4")
lncRNA_ENSG_ID_m<-rbind(lncRNA_ENSG_ID_m1,lncRNA_ENSG_ID_m2)
head(lncRNA_ENSG_ID_m)
#delete characters that are after . and :
lncRNA_ENSG_ID_m$N4 <- sub("[\\.:].*", "", lncRNA_ENSG_ID_m$N4)

#write.table(lncRNA_ENSG_ID_m,"DATA/HCC/lncRNAData/lncConvert2.txt", quote=F,sep="\t",row.names=FALSE,col.names = TRUE)





HCC_lncRNA_Symb_ID<-merge(x=HCC_lncRNA_df1, y=lncRNA_ENSG_ID_m, by.x = "HCC_lncRNA", by.y = "N4")
dim(HCC_lncRNA_Symb_ID)
head(HCC_lncRNA_Symb_ID)
HCC_lncRNA_Symb_ID<-merge(x=HCC_lncRNA_Symb_ID, y=lnc_mart_data_C, by.x = "N2", by.y = "Transcript.stable.ID")

HCC_lncID1<- data.frame(lncID=HCC_lncRNA_Symb_ID$Gene.stable.ID)
dim(HCC_lncID1)
###### lncRNA name and ID dataset 2 ######
HCC_lncRNA_Symb_ID2<-merge(x=HCC_lncRNA_df2, y=lncRNA_ENSG_ID_m, by.x = "result_dataframe3.transcript", by.y = "N3")
dim(HCC_lncRNA_Symb_ID2)
head(HCC_lncRNA_Symb_ID2)
HCC_lncRNA_Symb_ID2<-merge(x=HCC_lncRNA_Symb_ID2, y=lnc_mart_data_C, by.x = "N2", by.y = "Transcript.stable.ID")
dim(HCC_lncRNA_Symb_ID2)
head(HCC_lncRNA_Symb_ID2)

HCC_lncID2<- data.frame(lncID=HCC_lncRNA_Symb_ID2$Gene.stable.ID)
dim(HCC_lncID2)
head(HCC_lncID2)


HCC_lncID_All<-rbind(HCC_lncID1,HCC_lncID2)
HCC_lncID_All<-unique(HCC_lncID_All) #110
dim(HCC_lncID_All)
head(HCC_lncID_All)

#write.table(HCC_lncID_All,"DATA/HCC/lncRNAData/HCC_lncID_All.txt", quote=F,sep="\t",row.names=FALSE,col.names = TRUE)






