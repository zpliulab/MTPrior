# MTPrior: A multi-task hierarchical graph embedding framework for prioritizing Hepatocellular Carcinoma associated genes and long non-coding RNA

### Highlights:
A novel multi-task graph embedding model to prioritize disease-associated coding and noncoding RNAs, including lncRNAs and genes, by accounting for their interactions.  
The model adapts to various tasks and refines network structures to enhance prediction. 
Extensive validation shows its superior performance in identifying HCC-related genes and lncRNAs compared to state-of-the-art methods.

If you have any questions about MTPrior please directly contact the corresponding author Prof. Zhi-Ping Liu with the E-mail: zpliu@sdu.edu.cn

### Data:
In the Data directory, we provide only the essential files required to run the code. 
The datasets, which were obtained from publicly available sources, are introduced and referenced in the paper, along with the respective website links where they can be accessed.

### R and Python Code to MTPrior 
The serial number (1) (2) ... (10) represents the order in which the program runs in our work.

##### Prepare data:
*	1-Download Expression data from TCGA
*	2-Prepare Gene feature matrix
*	3-Prepare lncRNA feature matrix (Input the data-set contains ENSG ID and RNA Type (download: mart_export.txt))

##### Prepare Network:
*	4-Download and clean Interaction data (Gene-Gene, mRNA-miRNA, circRNA_miRNA, lncRNA-miRNA, ..., miRNA-HCC, circRNA-HCC,...)
	
##### Prepare inputs to embedding method:
*	5-Preapre Meta-path matrices
*	6-Lable nodes by HCC related Genes and HCC related lncRNAs
*	7-Prepare Train, Validation, Test index

##### Hirarchical embedding method:
*	8- Import meta-path matrices, feature matrix, label vector and index vector and 
		Run Hirarchical embedding code 

##### Rank by scores:
*	9-Rank genes and lncRNAs by scores gets as a result of hirarchical embedding method

##### Draw plots:
*	10 -Plot using different embedding size results, 
		Plot using different classification method, 
		Plot Ablation study results, 
		Plot distributation of results
