#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")

library(tidyr)
library(dplyr)
library(ggplot2)

data <- read.table("Code/Plots/Classification_genes_2.txt", header = TRUE)

data_long <- data %>%
  gather(key="Metric", value="Score", -Classification)

data_long$Classification <- factor(data_long$Classification, levels=unique(data$Classification))

data_long$Metric <- factor(data_long$Metric, levels=c("ACC", "AUROC", "AUPRC", "MCC", "F1","Sensitivity","Specificity"))


pCLAg<- ggplot(data_long, aes(x=Metric, y=Score, fill=Classification)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  labs(x=NULL,y="Scores", title=NULL) +
  theme_bw()+
  scale_fill_brewer(palette="YlGnBu",)+
  theme(text = element_text(family = "Times"),
        legend.position = "bottom",         # Position legend at the bottom
        legend.direction = "horizontal",    # Arrange legend items in a row
        legend.title = element_blank(),     # Remove the legend title
        legend.background = element_rect(color = "black", size = 0.3)
  )

#GnBu #YlGnBu
pCLAg

# pdf("Code/Plots/ClassificationGenes_2.pdf",width = 6,height =6)
# pCLAg
# dev.off()

svg("Code/Plots/ClassificationGenes_S.svg", width = 6, height = 6)
# Print the plot
print(pCLAg)
# Close the SVG device to save the file
dev.off()



#colors<- viridis::viridis(factor,option="D")




# Generate the custom color palette
# custom_palette <- colorRampPalette(c("yellow","darkgreen", "blue"))
# 
# # Assuming 'data_long' is your dataset and it has a 'Classification' column with unique levels
# unique_classifications <- unique(data_long$Classification)
# number_of_colors <- length(unique_classifications)
# 
# # Use the custom palette to generate colors
# colors <- custom_palette(number_of_colors)
# 
# pCLAg<-ggplot(data_long, aes(x=Metric, y=Score, fill=Classification)) +
#   geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
#   labs(x=NULL, y="Scores", title=NULL) +
#   theme_bw() +
#   scale_fill_manual(values = colors) +
#   theme(
#     legend.position = "bottom",         # Position legend at the bottom
#     legend.direction = "horizontal",    # Arrange legend items in a row
#     legend.title = element_blank(),     # Remove the legend title
#     legend.background = element_rect(color = "black", size = 0.3)
#   )

#pCLAg

pdf("Code/Plots/ClassificationGenes_4.pdf",width = 7,height =4)
pCLAg
dev.off()

###lncRNA

data <- read.table("Code/Plots/Classification_lncRNA_2.txt", header = TRUE)

data_long <- data %>%
  gather(key="Metric", value="Score", -Classification)

data_long$Classification <- factor(data_long$Classification, levels=unique(data$Classification))

data_long$Metric <- factor(data_long$Metric, levels=c("ACC", "AUROC", "AUPRC", "MCC", "F1","Sensitivity","Specificity"))


pCLAlnc<- ggplot(data_long, aes(x=Metric, y=Score, fill=Classification)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  labs(x=NULL,y="Scores", title=NULL) +
  theme_bw()+
  scale_fill_brewer(palette="YlGnBu",)+
  theme(text = element_text(family = "Times"),
        legend.position = "bottom",         # Position legend at the bottom
        legend.direction = "horizontal",    # Arrange legend items in a row
        legend.title = element_blank(),     # Remove the legend title
        legend.background = element_rect(color = "black", size = 0.3)
  )

# # Generate the custom color palette
# custom_palette <- colorRampPalette(c("lightblue", "blue"))
# 
# # Assuming 'data_long' is your dataset and it has a 'Classification' column with unique levels
# unique_classifications <- unique(data_long$Classification)
# number_of_colors <- length(unique_classifications)
# 
# # Use the custom palette to generate colors
# colors <- custom_palette(number_of_colors)
# 
# pCLAlnc<-ggplot(data_long, aes(x=Metric, y=Score, fill=Classification)) +
#   geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
#   labs(x=NULL, y="Scores", title=NULL) +
#   theme_bw() +
#   scale_fill_manual(values = colors) +
#   theme(
#     legend.position = "bottom",         # Position legend at the bottom
#     legend.direction = "horizontal",    # Arrange legend items in a row
#     legend.title = element_blank(),     # Remove the legend title
#     legend.background = element_rect(color = "black", size = 0.3)
#   )

pCLAlnc

pdf("Code/Plots/ClassificationlncRNA_4.pdf",width = 7,height =4)
pCLAlnc
dev.off()

svg("Code/Plots/ClassificationLNC_S.svg", width = 6, height = 6)
# Print the plot
print(pCLAlnc)
# Close the SVG device to save the file
dev.off()

