# This repository contains the different scripts that were used for the analysis of the scRNA-seq data of Wt and Dare 12 weeks mice.

## Starting from raw data
+ Scripts included in Step0 and Step1 folders are useful for the readers who wish to start from the raw data (*.fastq files).
+ If you proceed with the re-analysis of the raw data, please perform the necessary changes to the relative file paths described inside the scripts. 

## Starting from processed files uploaded in GEO repository of the study
1. In the current repository go to _**"<>Code"**_ button and select _**"Download ZIP"**_. 
2. Unzip the downloaded file and navigate inside _**"Iliopoulou_et_Al_Wt_Dare12wk_sc_analysis-main"**_ directory.
3. Create the sub-directories _**"GEO_files/PROCESSED/WT12/"**_ and _**"GEO_files/PROCESSED/DARE12/"**_.
4. Download the files _**wt12_matrix.mtx.gz, wt12_features.tsv.gz, wt12_barcodes.tsv.gz**_ from the GEO repository of the study and place them in the first sub-directory _**"GEO_files/PROCESSED/WT12/"**_.
5. Download the files _**dare12_matrix.mtx.gz, dare12_features.tsv.gz, dare12_barcodes.tsv.gz**_ from the GEO repository of the study and place them in the second sub-directory _**"GEO_files/PROCESSED/DARE12/"**_.
6. Rename the files that were placed in the sub-directories described above as: _**matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz**_.
7. Now you can execute sequentially the scripts in the folders _**Step2_InitialSeuratObjects, ..., Step7_Annotation_of_clusters**_.
8. Please note that you need to set the working directory to the folder that contains the script that is currently executed. 
