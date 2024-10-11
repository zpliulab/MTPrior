#set directory
setwd("/home/fatemeh/Fatemeh/0--SecoundProject/")

# Sample data
set.seed(123)  # for reproducibility
#labels <- data.frame(label = sample(0:1, 100, replace = TRUE))
labels <- read.table("DATA/HCC/TRUE_label.txt", header = F, sep = "\t")
labels <- data.frame(labels[,1])
colnames(labels)<-c("label")
sum(labels>0)


# Function to split indices into train, validation, and test sets
split_indices <- function(labels, train_prop , val_prop ) {
  n <- nrow(labels)
  n_class1 <- sum(labels$label == 1)
  
  # Calculate number of samples for each set in class 1
  train_size1 <- round(train_prop * n_class1)
  val_size1 <- round(val_prop * n_class1)
  test_size1 <- n_class1 - train_size1 - val_size1
  
  # Generate indices for class 1
  class1_indices <- which(labels$label == 1)
  class1_indices <- sample(class1_indices)
  
  # Split class 1 indices
  train_idx_class1 <- class1_indices[1:train_size1]
  val_idx_class1 <- class1_indices[(train_size1 + 1):(train_size1 + val_size1)]
  
  # Calculate remaining samples for test set
  remaining_idx <- setdiff(class1_indices, c(train_idx_class1, val_idx_class1))
  
  test_idx_class1 <- remaining_idx
  # If there are remaining indices, assign them to the test set
  # if (length(remaining_idx) > test_size) {
  #   test_idx_class1 <- remaining_idx[1:test_size]
  # } else {
  #   test_idx_class1 <- remaining_idx
  # }
  
  # Generate indices for class 0
  class0_indices <- which(labels$label == 0)
  n_class0 <- sum(labels$label == 0)
  # Calculate number of samples for each set in class 1
  train_size0 <- round(train_prop * n_class0)
  val_size0 <- round(val_prop * n_class0)
  test_size0 <- n_class0 - train_size0 - val_size0
  
  # Sample equal number of class 0 indices for each set
  train_idx_class0 <- sample(class0_indices, train_size0)
  val_idx_class0 <- sample(setdiff(class0_indices, train_idx_class0), val_size0)
  
  # Calculate remaining samples for test set
  remaining_idx0 <- setdiff(class0_indices, c(train_idx_class0, val_idx_class0))
  
  test_idx_class0 <- remaining_idx0
  
  # Combine indices for each set
  train_idx <- c(train_idx_class0, train_idx_class1)
  val_idx <- c(val_idx_class0, val_idx_class1)
  test_idx <- c(test_idx_class0,test_idx_class1)
  
  # Return indices
  return(list(train_idx = train_idx, val_idx = val_idx, test_idx = test_idx))
}

train_prop = 0.2 
val_prop = 0.1
# Split indices
indices <- split_indices(labels,train_prop, val_prop)

# Create train, validation, and test matrices
train_idx <- indices$train_idx
val_idx <- indices$val_idx
test_idx <- indices$test_idx

# Display the lengths of train, validation, and test indices
cat("Train set size:", length(train_idx), "\n")
cat("Validation set size:", length(val_idx), "\n")
cat("Test set size:", length(test_idx), "\n")

train.idx<-(t(data.frame(train_idx)))
val.idx<-(t(data.frame(val_idx)))
test.idx<-(t(data.frame(test_idx)))

train.idx<- data.frame(train.idx-1)
val.idx<- data.frame(val.idx-1)
test.idx<- data.frame(test.idx-1)

write.table(train.idx,"DATA/HCC/Ready_train.idxT.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
write.table(val.idx,"DATA/HCC/Ready_val.idxT.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)
write.table(test.idx,"DATA/HCC/Ready_test.idxT.txt", quote=F,sep="\t",row.names=FALSE,col.names = FALSE)



