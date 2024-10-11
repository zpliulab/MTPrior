R and Python Code to MTPrior 
The serial number (1) (2) ... (10) represents the order in which the program runs in our work.
(1)Prepare data:
	1-Download Expression data from TCGA
	2-Prepare Gene feature matrix
	3-Prepare lncRNA feature matrix (Input the data-set contains ENSG ID and RNA Type (download: mart_export.txt))
(2)Prepare Network:
	4-Download and clean Interaction data (Gene-Gene, mRNA-miRNA, circRNA_miRNA, lncRNA-miRNA, ..., miRNA-HCC, circRNA-HCC,...)
(3):Prepare inputs to embedding method:
	5-Preapre Meta-path matrices
	6-Lable nodes by HCC related Genes and HCC related lncRNAs
	7-Prepare Train, Validation, Test index

(5)Hirarchical embedding method:
	8- Import meta-path matrices, feature matrix, label vector and index vector and 
		Run Hirarchical embedding code 

(6)Rank by scores:
	9-Rank genes and lncRNAs by scores gets as a result of hirarchical embedding method

(7)Draw plots:
	10-Plot using different embedding size results, 
		Plot using different classification method, 
		Plot Ablation study results, 
		Plot distributation of results