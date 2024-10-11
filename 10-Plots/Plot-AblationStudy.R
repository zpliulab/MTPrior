#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")

library(ggplot2)

data_ab_g<- read.table("Code/Plots/Ablation_Gene.txt", header = TRUE)

# Add an index column to keep the order of the data
data_ab_g$Index <- rev(1:nrow(data_ab_g))





# Create the vertical line plot
pAblG<- ggplot(data_ab_g, aes(x = Index)) +
  annotate("rect", xmin = 1, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "#1E90FF") +
  annotate("rect", xmin = 2.5, xmax = 8.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#1E90FF") +
  annotate("rect", xmin = 8.5, xmax = 14, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#1E90FF") +
  geom_line(aes(y = AUROC, color = "AUROC")) +
  geom_point(aes(y = AUROC, color = "AUROC")) +
  geom_line(aes(y = ACC, color = "ACC")) +
  geom_point(aes(y = ACC, color = "ACC")) +
  coord_flip() + # This flips the coordinates to make the plot vertical
  theme_bw() +
  labs(title = NULL,
       x = NULL,
       y = NULL) +
  scale_color_manual(name = "Metrics", values = c("AUROC" = "#1E90FF", "ACC" =  "#2ECC71")) +
  scale_x_continuous(breaks = 1:nrow(data_ab_g))+ # Show all index values on the x-axis
  theme(text = element_text(family = "Times"),
        legend.position = "top",         # Position legend at the bottom
        legend.direction = "horizontal",    # Arrange legend items in a row
        legend.title = element_blank(),     # Remove the legend title
        legend.background = element_rect(color = "black", size = 0.3)
  )



pAblG


# 
# pdf("Code/Plots/Ablation_G_3.pdf",width = 3,height =6)
# pAblG
# dev.off()

svg("Code/Plots/Ablation_G_S.svg", width = 6, height = 6)

# Print the plot
print(pAblG)
# Close the SVG device to save the file
dev.off()


### lncRNA ######

data_ab_g<- read.table("Code/Plots/Ablation_lncRNA.txt", header = TRUE)

# Add an index column to keep the order of the data
data_ab_g$Index <- rev(1:nrow(data_ab_g))





# Create the vertical line plot
pAblL<- ggplot(data_ab_g, aes(x = Index)) +
  annotate("rect", xmin = 1, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.4, fill = "#1E90FF") +
  annotate("rect", xmin = 1.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "#1E90FF") +
  annotate("rect", xmin = 3.5, xmax = 8.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#1E90FF") +
  annotate("rect", xmin = 8.5, xmax = 13, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "#1E90FF") +
  geom_line(aes(y = AUROC, color = "AUROC")) +
  geom_point(aes(y = AUROC, color = "AUROC")) +
  geom_line(aes(y = ACC, color = "ACC")) +
  geom_point(aes(y = ACC, color = "ACC")) +
  coord_flip() + # This flips the coordinates to make the plot vertical
  theme_bw() +
  labs(title = NULL,
       x = NULL,
       y = NULL) +
  scale_color_manual(name = "Metrics", values = c("AUROC" = "#1E90FF", "ACC" =  "#2ECC71")) +
  scale_x_continuous(breaks = 1:nrow(data_ab_g))+
  theme(text = element_text(family = "Times"),
        legend.position = "top",         # Position legend at the bottom
        legend.direction = "horizontal",    # Arrange legend items in a row
        legend.title = element_blank(),     # Remove the legend title
        legend.background = element_rect(color = "black", size = 0.3)
  )



pAblL



# pdf("Code/Plots/Ablation_L_3.pdf",width = 3,height =6)
# pAblL
# dev.off()

svg("Code/Plots/Ablation_L_S.svg", width = 6, height = 6)

# Print the plot
print(pAblL)
# Close the SVG device to save the file
dev.off()


