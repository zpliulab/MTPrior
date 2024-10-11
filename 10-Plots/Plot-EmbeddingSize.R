

# Install and load the necessary packages
#install.packages("tidyr")
library(tidyr)
library(ggplot2)


#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")

data <- read.table("Code/Plots/EmbeddingSizeResults_Genes.txt", header = TRUE)



# Reshape the data from wide to long format
data_long <- pivot_longer(data, cols = -Embedding_size, names_to = "Metric", values_to = "Values")

data_long$Embedding_size <- factor(data_long$Embedding_size, levels = c(32, 64, 128, 256, 512, 1024))
#data_long$Metric<-factor(data_long$Metric,level=c("MCC","F1","ACC","AUPRC","AUROC"))

# Create the plot

custom_colors <- c("MCC" = "#f1c40f", "F1" = "#2ECC71","ACC" = "#f1c40f", "AUPRC" = "#2ECC71", "AUROC" = "#1E90FF")


pEMBg<- ggplot(data_long, aes(x = as.numeric(as.factor(Embedding_size)), y = Values, color = Metric, shape = Metric, group = Metric)) +
  geom_line(aes(group = Metric), size = 0.7) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:6, labels = c("32", "64", "128", "256", "512", "1024")) +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(x = "Embedding size", y = "Values")+ 
  theme_bw() +
  theme(text = element_text(family = "Times"),
        panel.grid.major = element_line(linetype =  "dashed", size = 0.3, color = "#A9A9A9"), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.95, 0.05),  # Positioning the legend inside the plot
        legend.justification = c(1, 0),
        legend.title = element_blank(),    # Remove the legend title
        legend.background = element_rect(color = "black", size = 0.3))+
  scale_shape_manual(values = c(19, 15, 18, 6, 17, 19))+scale_color_manual(values = custom_colors)
#scale_color_brewer(palette="Set2")+


pEMBg

# pdf("Code/Plots/EmbeddingGenes_2.pdf",width = 6,height =6)
# pEMBg
# dev.off()

svg("Code/Plots/EmbeddingGeneS.svg", width = 6, height = 6)
#pEMBg
# Print the plot
print(pEMBg)
# Close the SVG device to save the file
dev.off()

#ggsave("Code/Plots/EmbeddingGeneS.svg", plot = pEMBg, device = "svg", width = 6, height = 6)
#####################################

#install.packages("pheatmap")
#library(pheatmap)
data <- read.table("Code/Plots/EmbeddingSizeResults_Genes_2.txt", header = TRUE)
data$Embedding_size <- factor(data$Embedding_size, levels = c(1024,512,256,128,64,32))


# Reshape the data from wide to long format
# data_long <- pivot_longer(data, cols = -Embedding_size, names_to = "Metric", values_to = "Values")
# 
# data_long$Embedding_size <- factor(data_long$Embedding_size, levels = c(1024,512,256,128,64,32))
#data_long$Metric<-factor(data_long$Metric,level=c("MCC","F1","ACC","AUPRC","AUROC"))

data_m<-data[,-1]
rownames(data_m)<-data[,1]

# Specify the new order of rows
new_order <- rev(rownames(data_m))
reordered_data_matrix <- data_m[new_order, ]


# Plot the heatmap
hEMBg<- pheatmap(reordered_data_matrix, 
                 main = "",
                 color = colorRampPalette(c("yellow","green", "blue"))(50),
                 cluster_rows = FALSE, 
                 cluster_cols = FALSE,
                 show_rownames = TRUE, 
                 show_colnames = TRUE)

hEMBg

# pdf("Code/Plots/EmbeddingGenes_pheat_3.pdf",width = 4,height =4)
# hEMBg
# dev.off()

###########################################


data <- read.table("Code/Plots/EmbeddingSizeResults_lncRNA.txt", header = TRUE)

# Reshape the data from wide to long format
data_long <- pivot_longer(data, cols = -Embedding_size, names_to = "Metric", values_to = "Values")

data_long$Embedding_size <- factor(data_long$Embedding_size, levels = c(32, 64, 128, 256, 512, 1024))

# Create the plot
pEMBl<- ggplot(data_long, aes(x = as.numeric(as.factor(Embedding_size)), y = Values, color = Metric, shape = Metric, group = Metric)) +
  geom_line(aes(group = Metric), size = 0.7) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:6, labels = c("32", "64", "128", "256", "512", "1024")) +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(x = "Embedding size", y = "Values") +
  theme_bw() +
  theme(text = element_text(family = "Times"),
        panel.grid.major = element_line(linetype =  "dashed", size = 0.3, color = "#A9A9A9"), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.95, 0.05),  # Positioning the legend inside the plot
        legend.justification = c(1, 0),
        legend.title = element_blank(),    # Remove the legend title
        legend.background = element_rect(color = "black", size = 0.3))+
  scale_shape_manual(values = c(19, 15, 18, 6, 17, 19))+scale_color_manual(values = custom_colors)
#scale_color_brewer(palette="Set1")


pEMBl

# pdf("Code/Plots/EmbeddinglncRNA_2.pdf",width = 6,height =6)
# pEMBl
# dev.off()

svg("Code/Plots/EmbeddingLNCS.svg", width = 6, height = 6)
#pEMBg
# Print the plot
print(pEMBl)
# Close the SVG device to save the file
dev.off()


# data <- read.table("Code/Plots/EmbeddingSizeResults_lncRNA_2.txt", header = TRUE)
# 
# 
# data_m<-data[,-1]
# rownames(data_m)<-data[,1]
# 
# new_order <- rev(rownames(data_m))
# reordered_data_matrix <- data_m[new_order, ]
# 
# 
# # Plot the heatmap
# hEMBl<- pheatmap(reordered_data_matrix, 
#                  main = "",
#                  color = colorRampPalette(c("yellow","green", "blue"))(50),
#                  cluster_rows = FALSE, 
#                  cluster_cols = FALSE,
#                  show_rownames = TRUE, 
#                  show_colnames = TRUE)
# 
# hEMBl

# pdf("Code/Plots/Embeddinglnc_pheat_3.pdf",width = 4,height =4)
# hEMBl
# dev.off()

###########################################
library(gridExtra)
combined_plots_size<- grid.arrange(pEMBg, pEMBl, ncol = 2)

pdf("Code/Plots/Embeddingsize_grid.pdf",width = 9,height =6)
grid.arrange(pEMBg, pEMBl, ncol = 2)
dev.off()
##################################
library(ggplot2)
library(gridExtra)
library(grid)

# Create text grobs for A and B
label_A <- textGrob("A", x = unit(0.02, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
label_B <- textGrob("B", x = unit(0.52, "npc"), y = unit(0.95, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))

# Arrange the plots
combined_plots_size <- grid.arrange(pEMBg, pEMBl, ncol = 2)

# Combine plots and add labels
pdf("Code/Plots/Embeddingsize_grid.pdf", width = 9, height = 6)

grid.arrange(combined_plots_size, top = grid.rect(gp = gpar(col = "white")), ncol = 1) # To avoid title overwriting
grid.draw(label_A)
grid.draw(label_B)

dev.off()
