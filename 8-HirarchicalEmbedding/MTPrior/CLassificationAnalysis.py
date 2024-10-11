import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve
import numpy as np
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
#import xgboost as xgb

# Define the Multi-Layer Perceptron (MLP) model
class MLP(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, output_size)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        out = self.sigmoid(out)
        return out

# Custom Dataset class for loading feature and label matrices
class GeneDataset(Dataset):
    def __init__(self, features, labels):
        self.features = torch.tensor(features, dtype=torch.float32)
        self.labels = torch.tensor(labels, dtype=torch.float32)

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return self.features[idx], self.labels[idx]

# Function to train the model
def train_model(model, train_loader, criterion, optimizer, num_epochs=10):
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
        print(f"Epoch {epoch+1}/{num_epochs}, Loss: {running_loss/len(train_loader)}")

# Function to calculate accuracy
def calculate_accuracy(predictions, targets):
    binary_predictions = torch.round(predictions)
    correct = (binary_predictions == targets).sum().item()
    total = targets.size(0) * targets.size(1)
    accuracy = correct / total
    return accuracy
def calculate_accuracy_s(predictions, targets):
    binary_predictions = torch.round(predictions)
    correct = (binary_predictions == targets).sum().item()
    total = targets.size(0)  # Only use the first dimension
    accuracy = correct / total
    return accuracy

# Function to calculate AUC
def calculate_auc(predictions, targets):
    fpr, tpr, thresholds = roc_curve(targets.ravel(), predictions.ravel())
    auc = roc_auc_score(targets, predictions)
    return auc

'''
# Function to calculate AUC with threshold
def calculate_auc(predictions, targets, threshold=0.9):
    binary_predictions = torch.where(predictions >= threshold, 1, 0)
    fpr, tpr, thresholds = roc_curve(targets.ravel(), binary_predictions.ravel())
    auc = roc_auc_score(targets, binary_predictions)
    return auc
'''
def calculate_f1(predictions, targets, threshold=0.5):
    binary_predictions = torch.where(predictions >= threshold, 1, 0)
    tp = torch.sum(binary_predictions * targets).item()
    fp = torch.sum(binary_predictions * (1 - targets)).item()
    fn = torch.sum((1 - binary_predictions) * targets).item()
    precision = tp / (tp + fp + 1e-10)
    recall = tp / (tp + fn + 1e-10)
    f1 = 2 * (precision * recall) / (precision + recall + 1e-10)
    return f1

from sklearn.metrics import precision_recall_curve, auc
def calculate_auprc(predictions, targets):
    precision, recall, _ = precision_recall_curve(targets.ravel(), predictions.ravel())
    auprc_score = auc(recall, precision)  # Rename variable to avoid conflict
    return auprc_score

def calculate_mcc(predictions, targets, threshold=0.5):
    binary_predictions = torch.where(predictions >= threshold, 1, 0)
    tp = torch.sum(binary_predictions * targets).item()
    tn = torch.sum((1 - binary_predictions) * (1 - targets)).item()
    fp = torch.sum(binary_predictions * (1 - targets)).item()
    fn = torch.sum((1 - binary_predictions) * targets).item()
    mcc = (tp * tn - fp * fn) / ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    return mcc
# Function to calculate true negatives, false positives, sensitivity, specificity
def calculate_metrics(predictions, targets, threshold=0.5):
    binary_predictions = torch.where(predictions >= threshold, 1, 0)
    tp = torch.sum(binary_predictions * targets).item()
    tn = torch.sum((1 - binary_predictions) * (1 - targets)).item()
    fp = torch.sum(binary_predictions * (1 - targets)).item()
    fn = torch.sum((1 - binary_predictions) * targets).item()
    sensitivity = tp / (tp + fn + 1e-10)
    specificity = tn / (tn + fp + 1e-10)
    return tn, fp, sensitivity, specificity

def RClassification(feature_matrix, label_matrix):


    feature_matrix = feature_matrix
    label_matrix = label_matrix


    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(feature_matrix, label_matrix, test_size=0.2, random_state=42, stratify=label_matrix)

    # Calculate class weights
    class_counts = np.sum(y_train, axis=0)
    total_samples = len(y_train)
    class_weights = torch.tensor([total_samples / (len(class_counts) * count) for count in class_counts], dtype=torch.float32)

    # Define model parameters
    input_size = len(X_train[0])
    hidden_size =  64
    output_size = len(y_train[0])

    # Create model, criterion, optimizer
    model = MLP(input_size, hidden_size, output_size)
    criterion = nn.BCELoss(weight=class_weights)
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Create DataLoader for training data
    train_dataset = GeneDataset(X_train, y_train)
    train_loader = DataLoader(train_dataset, batch_size=1, shuffle=True)

    # Train the model
    train_model(model, train_loader, criterion, optimizer)

    # Evaluate the model on test data
    model.eval()
    with torch.no_grad():
        test_input = torch.tensor(X_test, dtype=torch.float32)
        test_labels = torch.tensor(y_test, dtype=torch.float32)
        predictions = model(test_input)

        accuracy = calculate_accuracy(predictions, test_labels)
        print("Accuracy:", accuracy)

        auc = calculate_auc(predictions, test_labels)
        print("AUC:", auc)


        auprc = calculate_auprc(predictions, test_labels)
        print("AUPRC:", auprc)

        mcc = calculate_mcc(predictions, test_labels)
        print("MCC:", mcc)

        f1 = calculate_f1(predictions, test_labels)
        print("F1 Score:", f1)

        # Calculate true negatives, false positives, sensitivity, specificity
        tn, fp, sensitivity, specificity = calculate_metrics(predictions, test_labels)
        print("True Negatives:", tn)
        print("False Positives:", fp)
        print("Sensitivity:", sensitivity)
        print("Specificity:", specificity)

        # Prepare data for sensitivity vs specificity plot
        fpr, tpr, thresholds = roc_curve(test_labels.ravel(), predictions.ravel())
        # Save data to CSV files
        np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr_g_5.csv", fpr, delimiter=',')
        np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr_g_5.csv", tpr, delimiter=',')
        np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds_g_5.csv", thresholds, delimiter=',')

        # Convert predictions tensor to numpy array
        predictions_np = predictions.numpy()
        # Extract probabilities of being in class 1
        # class_1_probabilities = predictions_np[:, 1]  # Assuming class 1 is at index 1
        class_1_probabilities = predictions_np[:, 0]  # Assuming class 1 is at index 0
        #print(class_1_probabilities)
        # Specify the file path
        file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_g_5.txt"
        # Write the array to the file
        # np.savetxt(file_path_export_rank, class_1_probabilities, fmt='%d')
        np.savetxt(file_path_export_rank, class_1_probabilities)

        ##For Consider similar rank for similar prediction value
        # Round the probabilities to two decimal places
        rounded_probabilities = np.round(class_1_probabilities, decimals=2)

        # Rank inputs based on probabilities
        ranked_indices = np.argsort(rounded_probabilities)[::-1]  # Descending order

        # Initialize rank counter
        current_rank = 1

        # Initialize previous probability
        previous_probability = None

        # Initialize list to store ranks
        ranks = []

        # Assign ranks
        for idx in ranked_indices:
            probability = rounded_probabilities[idx]
            # Check if the current probability is different from the previous one
            if probability != previous_probability:
                # If different, assign the current rank
                ranks.append(current_rank)
                current_rank += 1
            else:
                # If same, assign the same rank as the previous one
                ranks.append(current_rank - 1)
            # Update previous probability
            previous_probability = probability

        # Convert ranks list to numpy array
        ranks = np.array(ranks)

        # Specify the file path
        file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_g_5.txt"
        # Write the array to the file
        np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

        return auc

        # Prepare data for sensitivity vs specificity plot
        fpr, tpr, thresholds = roc_curve(test_labels.ravel(), predictions.ravel())
        #print("fpr",fpr)
        #print("tpr",tpr)
        #print("thresholds",thresholds)
        # Save data to CSV files
        #np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr7.csv", fpr, delimiter=',')
        #np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr7.csv", tpr, delimiter=',')
        #np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds7.csv", thresholds, delimiter=',')

def SVMClassification(feature_matrix, label_matrix):
    # Convert one-hot encoded labels to class labels
    label_matrix = np.argmax(label_matrix, axis=1)


    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(feature_matrix, label_matrix, test_size=0.2, random_state=42, stratify=label_matrix)

    # Initialize SVM model with class weighting
    #class_weights = {i: len(y_train) / (len(np.unique(y_train)) * np.sum(y_train == i)) for i in np.unique(y_train)}
    #model = SVC(kernel='linear', class_weight=class_weights, probability=True)
    model = SVC(kernel='linear', probability=True)

    # Train the SVM model
    model.fit(X_train, y_train)

    # Make predictions
    predictions = model.predict_proba(X_test)[:, 1]

    # Calculate metrics
    accuracy = calculate_accuracy_s(torch.tensor(predictions), torch.tensor(y_test))
    auc = calculate_auc(torch.tensor(predictions), torch.tensor(y_test))
    auprc = calculate_auprc(predictions, y_test)
    mcc = calculate_mcc(torch.tensor(predictions), torch.tensor(y_test))
    f1 = calculate_f1(torch.tensor(predictions), torch.tensor(y_test))
    tn, fp, sensitivity, specificity = calculate_metrics(torch.tensor(predictions), torch.tensor(y_test))

    print("SVM Results")
    print("Accuracy:", accuracy)
    print("AUC:", auc)
    print("AUPRC:", auprc)
    print("MCC:", mcc)
    print("F1 Score:", f1)
    print("True Negatives:", tn)
    print("False Positives:", fp)
    print("Sensitivity:", sensitivity)
    print("Specificity:", specificity)

    # Prepare data for sensitivity vs specificity plot
    fpr, tpr, thresholds = roc_curve(y_test.ravel(), predictions.ravel())

    # Save data to CSV files
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr_svm.csv", fpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr_svm.csv", tpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds_svm.csv", thresholds, delimiter=',')

    # Convert predictions to numpy array
    predictions_np = predictions

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_svm.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, predictions_np)

    # Rank inputs based on probabilities
    ranked_indices = np.argsort(predictions_np)[::-1]  # Descending order

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_svm.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

    return auc



def RFClassification(feature_matrix, label_matrix):
    # Convert one-hot encoded labels to class labels if necessary
    if len(label_matrix.shape) > 1 and label_matrix.shape[1] > 1:
        label_matrix = np.argmax(label_matrix, axis=1)

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(feature_matrix, label_matrix, test_size=0.2, random_state=42, stratify=label_matrix)

    # Initialize Random Forest model
    model = RandomForestClassifier(n_estimators=100, random_state=42)

    # Train the Random Forest model
    model.fit(X_train, y_train)

    # Make predictions
    predictions = model.predict_proba(X_test)[:, 1]

    # Calculate metrics
    accuracy = calculate_accuracy_s(torch.tensor(predictions), torch.tensor(y_test, dtype=torch.float32))
    auc_score = calculate_auc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    auprc = calculate_auprc(predictions, y_test)
    mcc = calculate_mcc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    f1 = calculate_f1(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    tn, fp, sensitivity, specificity = calculate_metrics(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))

    print("Random Forest Results")
    print("Accuracy:", accuracy)
    print("AUC:", auc_score)
    print("AUPRC:", auprc)
    print("MCC:", mcc)
    print("F1 Score:", f1)
    print("True Negatives:", tn)
    print("False Positives:", fp)
    print("Sensitivity:", sensitivity)
    print("Specificity:", specificity)

    # Prepare data for sensitivity vs specificity plot
    fpr, tpr, thresholds = roc_curve(y_test.ravel(), predictions.ravel())

    # Save data to CSV files
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr_rf.csv", fpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr_rf.csv", tpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds_rf.csv", thresholds, delimiter=',')

    # Convert predictions to numpy array
    predictions_np = predictions

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_rf.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, predictions_np)

    # Rank inputs based on probabilities
    ranked_indices = np.argsort(predictions_np)[::-1]  # Descending order

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_rf.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

    return auc_score

def LRClassification(feature_matrix, label_matrix):
    # Convert one-hot encoded labels to class labels if necessary
    if len(label_matrix.shape) > 1 and label_matrix.shape[1] > 1:
        label_matrix = np.argmax(label_matrix, axis=1)

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(feature_matrix, label_matrix, test_size=0.2, random_state=42, stratify=label_matrix)

    # Initialize Logistic Regression model
    model = LogisticRegression(random_state=42, max_iter=10000)

    # Train the model
    model.fit(X_train, y_train)

    # Make predictions
    predictions = model.predict_proba(X_test)[:, 1]

    # Calculate metrics
    accuracy = calculate_accuracy_s(torch.tensor(predictions), torch.tensor(y_test, dtype=torch.float32))
    auc_score = calculate_auc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    auprc = calculate_auprc(predictions, y_test)
    mcc = calculate_mcc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    f1 = calculate_f1(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    tn, fp, sensitivity, specificity = calculate_metrics(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))

    print("Logistic Regression Results")
    print("Accuracy:", accuracy)
    print("AUC:", auc_score)
    print("AUPRC:", auprc)
    print("MCC:", mcc)
    print("F1 Score:", f1)
    print("True Negatives:", tn)
    print("False Positives:", fp)
    print("Sensitivity:", sensitivity)
    print("Specificity:", specificity)

    # Prepare data for sensitivity vs specificity plot
    fpr, tpr, thresholds = roc_curve(y_test.ravel(), predictions.ravel())

    # Save data to CSV files
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr_lr.csv", fpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr_lr.csv", tpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds_lr.csv", thresholds, delimiter=',')

    # Convert predictions to numpy array
    predictions_np = predictions

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_lr.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, predictions_np)

    # Rank inputs based on probabilities
    ranked_indices = np.argsort(predictions_np)[::-1]  # Descending order

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_lr.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

    return auc_score


def AdaBoostClassification(feature_matrix, label_matrix):
    # Convert one-hot encoded labels to class labels if necessary
    if len(label_matrix.shape) > 1 and label_matrix.shape[1] > 1:
        label_matrix = np.argmax(label_matrix, axis=1)

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(feature_matrix, label_matrix, test_size=0.2, random_state=42, stratify=label_matrix)

    # Initialize AdaBoost model
    model = AdaBoostClassifier(n_estimators=50, random_state=42)

    # Train the model
    model.fit(X_train, y_train)

    # Make predictions
    predictions = model.predict_proba(X_test)[:, 1]

    # Calculate metrics
    accuracy = calculate_accuracy_s(torch.tensor(predictions), torch.tensor(y_test, dtype=torch.float32))
    auc_score = calculate_auc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    auprc = calculate_auprc(predictions, y_test)
    mcc = calculate_mcc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    f1 = calculate_f1(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    tn, fp, sensitivity, specificity = calculate_metrics(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))

    print("AdaBoost Classifier Results")
    print("Accuracy:", accuracy)
    print("AUC:", auc_score)
    print("AUPRC:", auprc)
    print("MCC:", mcc)
    print("F1 Score:", f1)
    print("True Negatives:", tn)
    print("False Positives:", fp)
    print("Sensitivity:", sensitivity)
    print("Specificity:", specificity)

    # Prepare data for sensitivity vs specificity plot
    fpr, tpr, thresholds = roc_curve(y_test.ravel(), predictions.ravel())

    # Save data to CSV files
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr_ab.csv", fpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr_ab.csv", tpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds_ab.csv", thresholds, delimiter=',')

    # Convert predictions to numpy array
    predictions_np = predictions

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_ab.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, predictions_np)

    # Rank inputs based on probabilities
    ranked_indices = np.argsort(predictions_np)[::-1]  # Descending order

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_ab.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

    return auc_score

from sklearn.naive_bayes import GaussianNB

def NaiveBayesClassification(feature_matrix, label_matrix):
    # Convert one-hot encoded labels to class labels if necessary
    if len(label_matrix.shape) > 1 and label_matrix.shape[1] > 1:
        label_matrix = np.argmax(label_matrix, axis=1)

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(feature_matrix, label_matrix, test_size=0.2, random_state=42, stratify=label_matrix)

    # Initialize Naive Bayes model
    model = GaussianNB()

    # Train the model
    model.fit(X_train, y_train)

    # Make predictions
    predictions = model.predict_proba(X_test)[:, 1]

    # Calculate metrics
    accuracy = calculate_accuracy_s(torch.tensor(predictions), torch.tensor(y_test, dtype=torch.float32))
    auc_score = calculate_auc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    auprc = calculate_auprc(predictions, y_test)
    mcc = calculate_mcc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    f1 = calculate_f1(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    tn, fp, sensitivity, specificity = calculate_metrics(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))

    print("Naive Bayes Classifier Results")
    print("Accuracy:", accuracy)
    print("AUC:", auc_score)
    print("AUPRC:", auprc)
    print("MCC:", mcc)
    print("F1 Score:", f1)
    print("True Negatives:", tn)
    print("False Positives:", fp)
    print("Sensitivity:", sensitivity)
    print("Specificity:", specificity)

    # Prepare data for sensitivity vs specificity plot
    fpr, tpr, thresholds = roc_curve(y_test.ravel(), predictions.ravel())

    # Save data to CSV files
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr_nb.csv", fpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr_nb.csv", tpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds_nb.csv", thresholds, delimiter=',')

    # Convert predictions to numpy array
    predictions_np = predictions

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_nb.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, predictions_np)

    # Rank inputs based on probabilities
    ranked_indices = np.argsort(predictions_np)[::-1]  # Descending order

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_nb.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

    return auc_score


def XGBoostClassification(feature_matrix, label_matrix):
    # Convert one-hot encoded labels to class labels if necessary
    if len(label_matrix.shape) > 1 and label_matrix.shape[1] > 1:
        label_matrix = np.argmax(label_matrix, axis=1)

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(feature_matrix, label_matrix, test_size=0.2, random_state=42, stratify=label_matrix)

    # Initialize XGBoost model
    model = xgb.XGBClassifier(n_estimators=100, random_state=42)

    # Train the model
    model.fit(X_train, y_train)

    # Make predictions
    predictions = model.predict_proba(X_test)[:, 1]

    # Calculate metrics
    accuracy = calculate_accuracy_s(torch.tensor(predictions), torch.tensor(y_test, dtype=torch.float32))
    auc_score = calculate_auc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    auprc = calculate_auprc(predictions, y_test)
    mcc = calculate_mcc(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    f1 = calculate_f1(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))
    tn, fp, sensitivity, specificity = calculate_metrics(torch.tensor(predictions, dtype=torch.float32), torch.tensor(y_test, dtype=torch.float32))

    print("XGBoost Classifier Results")
    print("Accuracy:", accuracy)
    print("AUC:", auc_score)
    print("AUPRC:", auprc)
    print("MCC:", mcc)
    print("F1 Score:", f1)
    print("True Negatives:", tn)
    print("False Positives:", fp)
    print("Sensitivity:", sensitivity)
    print("Specificity:", specificity)

    # Prepare data for sensitivity vs specificity plot
    fpr, tpr, thresholds = roc_curve(y_test.ravel(), predictions.ravel())

    # Save data to CSV files
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr_xg.csv", fpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr_xg.csv", tpr, delimiter=',')
    np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds_xg.csv", thresholds, delimiter=',')

    # Convert predictions to numpy array
    predictions_np = predictions

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_xg.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, predictions_np)

    # Rank inputs based on probabilities
    ranked_indices = np.argsort(predictions_np)[::-1]  # Descending order

    # Specify the file path
    file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_xg.txt"
    # Write the array to the file
    np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

    return auc_score

'''
if __name__ == "__main__":

    #file_path_embedding = "/home/fatemeh/MTPrior/data/HCC/Export/Emb_g_test_256.txt"
    file_path_embedding = "/home/fatemeh/MTPrior/data/HCC/Export/Emb_lnc_512.txt"

    #file_path_label = "/home/fatemeh/MTPrior/data/HCC/Import_filterAgain2/Ready_label.txt"
    file_path_label = "/home/fatemeh/MTPrior/data/HCC/lncRNA-import/PN2/PN_label_2.txt"

    feature_matrix = np.loadtxt(file_path_embedding)
    label_matrix = np.loadtxt(file_path_label)

    RClassification(feature_matrix, label_matrix)
    #SVMClassification(feature_matrix, label_matrix)
    #RFClassification(feature_matrix, label_matrix)
    #LRClassification(feature_matrix, label_matrix)
    #AdaBoostClassification(feature_matrix, label_matrix)
    #XGBoostClassification(feature_matrix, label_matrix)
    #NaiveBayesClassification(feature_matrix, label_matrix)

'''

'''
# prediction for whole dataset according to trained model
X_model = feature_matrix
y_model = label_matrix

if __name__ == "__main__":
    model.eval()
    with torch.no_grad():
        All_input = torch.tensor(X_model, dtype=torch.float32)
        # Convert y_test to tensor
        All_labels = torch.tensor(y_model, dtype=torch.float32)
        predictions = model(All_input)

        All_accuracy = calculate_accuracy(predictions, All_labels)
        print("Accuracy All:", All_accuracy)

        #threshold = 0.9  # Define your threshold value
        #auc = calculate_auc(predictions, All_labels, threshold=threshold)
        auc = calculate_auc(predictions, All_labels)
        print("AUC:", auc)
        #print("predictions", predictions[5:])

        f1 = calculate_f1(predictions, All_labels)
        print("F1 Score:", f1)

        #auprc = calculate_auprc(predictions, test_labels)
        #print("AUPRC:", auprc)

        mcc = calculate_mcc(predictions, All_labels)
        print("MCC:", mcc)

        # Calculate true negatives, false positives, sensitivity, specificity
        tn, fp, sensitivity, specificity = calculate_metrics(predictions, All_labels)
        print("True Negatives:", tn)
        print("False Positives:", fp)
        print("Sensitivity:", sensitivity)
        print("Specificity:", specificity)

        # Prepare data for sensitivity vs specificity plot
        fpr, tpr, thresholds = roc_curve(All_labels.ravel(), predictions.ravel())
        #print("fpr", fpr)
        #print("tpr", tpr)
        #print("thresholds", thresholds)
        # Save data to CSV files
        #np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/fpr3.csv", fpr, delimiter=',')
        #np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/tpr3.csv", tpr, delimiter=',')
        #np.savetxt("/home/fatemeh/MTPrior/data/HCC/Export/thresholds3.csv", thresholds, delimiter=',')

# Convert predictions tensor to numpy array
predictions_np = predictions.numpy()

# Extract probabilities of being in class 1
#class_1_probabilities = predictions_np[:, 1]  # Assuming class 1 is at index 1
class_1_probabilities = predictions_np[:, 0]  # Assuming class 1 is at index 0
print(class_1_probabilities)
# Specify the file path
file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/probabilities_All_ClassWeighting_Again7.txt"
# Write the array to the file
#np.savetxt(file_path_export_rank, class_1_probabilities, fmt='%d')
np.savetxt(file_path_export_rank, class_1_probabilities)


# Rank inputs based on probabilities
ranked_indices = np.argsort(class_1_probabilities)[::-1]  # Descending order

# Now, ranked_indices contains the indices of inputs sorted based on the probability of being in class 1
# You can use these indices to access corresponding inputs in feature_matrix
ranked_inputs = X_model[ranked_indices]

#print(ranked_indices)

# Specify the file path
##file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_indices_All_ClassWeighting_Again7.txt"
# Write the array to the file
##np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

##For Consider similar rank for similar prediction value
# Round the probabilities to two decimal places
rounded_probabilities = np.round(class_1_probabilities, decimals=2)

# Rank inputs based on probabilities
ranked_indices = np.argsort(rounded_probabilities)[::-1]  # Descending order

# Initialize rank counter
current_rank = 1

# Initialize previous probability
previous_probability = None

# Initialize list to store ranks
ranks = []

# Assign ranks
for idx in ranked_indices:
    probability = rounded_probabilities[idx]
    # Check if the current probability is different from the previous one
    if probability != previous_probability:
        # If different, assign the current rank
        ranks.append(current_rank)
        current_rank += 1
    else:
        # If same, assign the same rank as the previous one
        ranks.append(current_rank - 1)
    # Update previous probability
    previous_probability = probability

# Convert ranks list to numpy array
ranks = np.array(ranks)

# Specify the file path
##file_path_export_rank = "/home/fatemeh/MTPrior/data/HCC/Export/ranked_indices_All_ClassWeighting_SimilarRank_Again7.txt"
# Write the array to the file
##np.savetxt(file_path_export_rank, ranked_indices, fmt='%d')

'''