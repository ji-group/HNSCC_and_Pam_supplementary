# HNSCC\_and\_Pam\_supplementary
Supplementary materials for manuscript "Identification of a novel palmitoylation-related prognostic signature in head and neck squamous cell carcinoma"

To reproduce all analyses, please follow these steps in a command line environment:

1. Git clone this repository

	```
	git clone https://github.com/ji-group/HNSCC_and_Pam_supplementary.git
	cd HNSCC_and_Pam_supplementary
	```
2. Pull large files

	```
	git lfs pull
	```
3.	Unzip large files

	```
	unzip DataFiles.zip
	unzip HNSC.zip
	```
4. Run R script from command line

	```
	Rscript code.R
	```
	
If you run into any errors (especially for missing any packages), please try to install them either through CRAN or Bioconductor.
There are commented commands in the code.R script for installing most of these packages if needed.

