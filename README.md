# Smutans-MUC5B-RNASeq
Instructions on how to replicate analyses from Werlang et. al. Nature Microbiology (2021) 
https://dx.doi.org/10.1038/s41564-021-00876-1

## Step 0: Map reads (BAM files) onto the UA159 genome using usegalaxy.org
### You can download pre-mapped read files from https://doi.org/10.6084/m9.figshare.14214131 (to skip substeps 0-2)
0. The raw bam files are available from the Gene Expression Omnibus under accession number GSE163258 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163258
2. Load the 'data/external/GCF_000007465.2_ASM746v2_genomic.fna' file and raw bam files into Galaxy
3. Map with BWA (Burrows-Wheeler Aligner). Options: unpaired single end short reads
4. Download the .bam files into the "raw_mapped_bam_from_Galaxy" folder
5. Rename the files to 'xxxx.mapped.bam'
6. Fill in the 'data/raw_mapped_bam_from_Galaxy/file_list_mapped_bam.csv' file with the new filenames

## Step 1: Count your mapped reads per gene using R 
0. The R script pulls from the following files: 
	- 'data/external/file_list_external_data.csv'
	- 'data/raw_mapped_bam_from_Galaxy/file_list_mapped_bam.csv'
1. Run "STEP1_ReadCounting_withComS.R"
2. The script outputs "interim/raw_read_counts.csv" 

## Step 2: Analyze relative expression between conditions using MATLAB
0. The MATLAB script pulls from the following files:
	- 'experiment_design_comparison_matrix.csv'
	Note: This file determines which comparisons will run (1=True, so a comparison will be run). 
	Note: The columns are the "treatment" and rows are "control" in the output. 
	- 'experiment_design_replicate_matrix.csv'
	Note: The file order should match that in 'data/raw_mapped_bam_from_Galaxy/file_list_mapped_bam.csv'
	Note: For the replicate matrix, 1s represent matched biological replicates. 
	- 'interim/raw_read_counts.csv'
2. Run "STEP2_RelativeExpressionAnalysis.m"
3. The script outputs several "RelExp_XvsY.csv" files into the "interim" folder, with a list of their names in "file_list_relative_expression.csv"

## Step 3: Organize data & do Pathway Enrichment Analysis using Jupyter Notebooks
0. The notebook pulls from the following files:
	- 'data/external/GCF_000007465.2_ASM746v2_genomic.gff3' gene index file
	- 'data/external/smu00001.json' KEGG pathway file 
	- The relative expression files output by MATLAB in the interim folder
1. Execute "STEP3_PathwayEnrichmentAnalysis.ipynb"
2. Files are outputted into the "processed" folder
