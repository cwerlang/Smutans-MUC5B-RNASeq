## Read countng for data from Werlang et. al. Nat. Micro. (2021)
## DOI: 10.1038/s41564-021-00876-1
## Written by C. Werlang and K. Wheeler

## install and update packages
install.packages("BiocManager")
BiocManager::install() #install core packages
BiocManager::install(c("ShortRead","easyRNASeq","Rsamtools",
                         "biomaRt","GenomicFeatures","GenomicAlignments","BiocParallel",
                         "DESeq2","magrittr","stringi","GenomicRanges"))
install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock')) 

n## load libraries
library("BiocParallel")
library("Rsamtools")
library(rstudioapi)   
library("GenomicAlignments")
library("GenomicFeatures")
library("GenomicRanges")

## ----- Set Directory --------
## If you save this R script as the same folder as your mapped ...
## bam files, then you don't have to change anything here.
mainwd <-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mainwd)

# set GFF file name
organisminfo <- read.csv(file.path('data','external','file_list_external_data.csv'),row.names=1, stringsAsFactors = FALSE)
GFF_FILE <- file.path('data','external',organisminfo$filename[[1]]) #this pulls the gff3 actually
file.exists(GFF_FILE)

## ----- Set Filenames & Load Samples-------
## the Galaxy output (mapped bam files) should have been renamed
## and match the sample names in the csv file 
FILE_LIST <- file.path('data','raw_mapped_bam_from_Galaxy', 'file_list_mapped_bam.csv') 
OUTPUT_FILE <- file.path('data','interim','raw_read_counts.csv')
BAMFILE_LOCATION <- file.path('data','raw_mapped_bam_from_Galaxy')

# change csv file name here if needed, also make sure that Run column has correct names
samplefile <- file.path(FILE_LIST)
sampleTable <- read.csv(samplefile, row.names = 1)
sampleTable

# check that our files were loaded correctly
filenames <- file.path(BAMFILE_LOCATION, paste0(sampleTable$sample, ".mapped.bam"))
file.exists(filenames)

# load the bam files as bam files (Rsamtools)
bamfiles <- BamFileList(filenames, yieldSize=100000)
# show sequence info of bam files
seqinfo(bamfiles[1])

## ----txdb----------------------------------------------------------------
# upload your gff, gtf, or gff3 file
# these have your genome indices annotated with gene name/locus
# make sure to change the "format" to match your file type

txdb <- makeTxDbFromGFF(GFF_FILE, format = "gff3", circ_seqs = character())
txdb

## ------------------------------------------------------------------------
cbt <- cdsBy(txdb, "tx") #converts txdb to CompressedGRangesList 
cbt

ids<-id2name(txdb,feature.type="cds")
head(ids)

cbt@partitioning@NAMES <-ids
cbt@elementMetadata@rownames <- ids

## --- This is special code to add in ComS to the UA159 genome #ignore in other cases

gr <- GRanges(
  seqnames = Rle(c("NC_004350.2"),c(1)),
  ranges = IRanges(62583,62663),
  strand = Rle(strand(c("+")), c(1)),
  cds_id = 3000,
  cds_name = "comS",
  exon_rank = 1
)

gr@elementMetadata@rownames <- "comS"

## Map sequencing reads to genome------------------------------------------------------------------------

se <- summarizeOverlaps(features=c(cbt,gr), reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=TRUE,
                        fragments=FALSE, param=ScanBamParam())

se
dim(se)
assayNames(se)
head(assay(se), 12)
colSums(assay(se))
###

## write out mapped unnormalized sequencing reads
write.csv(assay(se),OUTPUT_FILE)


