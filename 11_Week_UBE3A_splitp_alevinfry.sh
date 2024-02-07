rf# Build splici index for USA mode quantification
# Include extra mito sequences as these are absent in the gencode reference
# Downloaded from: https://zenodo.org/record/5799568 on 2/24/22
# In R:

library(roe)
make_splici_txome(
  gtf_path = "/Documents/genome_annotation/gencode.v39.annotation.gtf",
  genome_path = "/Documents/genome_annotation/hg38.fa",
  read_length = 100,
  flank_trim_length = 5,
  output_dir = "/Documents/genome_annotation/",
  file_name_prefix = "gencode_v39_transcriptome_splici",
  extra_spliced = "/Documents/genome_annotation/homo_sapiens_mito_seqs_for_splici.fasta",
  dedup_seqs = FALSE
)
grep("gencode_v39_transcriptome_splici", dir("/Documents/genome_annotation/"), value = TRUE)
# Created a splici FASTA and a 3-column tx2gene table:
#    "gencode_v39_transcriptome_splici_fl95_t2g_3col.tsv"
#    "gencode_v39_transcriptome_splici_fl95.fa"

# Build salmon index
#!/bin/tcsh

#BSUB -n 12
#BSUB -W 5760
#BSUB -R "rusage[mem=128GB]"

#Activate salmon conda environment
conda activate /usr/local/usrapps/salmon/conda-for-salmon
salmon index -t /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95.fa -p 12 -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx
#Deactivate salmon conda environment
conda deactivate

#create genome annotation lookup file from https://ftp.ebi.ac.uk/pub/databases/gencode/ gtf file
#change to directory with gtf file
#do following line in terminal
grep -v '#' gencode.v39.annotation.gtf | cut -f 9 | grep transcript_id | cut -d ';' -f 1,2,4 | cut -d '"' -f 2,4,6 | tr '"' '\t' | sort | uniq > gencode.v39.annotation.lookup.tsv

#In R
#read in gencode.v39.annotation.lookup.tsv and extract info
annotation_file <- read_tsv("/Documents/genome_annotation/gencode.v39.annotation.lookup.tsv",show_col_types = FALSE,col_names = FALSE)
annotation_file <- distinct(annotation_file[,c(1,3)])
#add mitochondrial genes
mito_genes <- data.frame(c("ENSG00000210082","ENSG00000211459"),c("MT-RNR2","MT-RNR1"))
colnames(mito_genes)<-colnames(annotation_file)
annotation_file2 <- rbind(annotation_file,mito_genes)
write.table(annotation_file2,"gencode.v39.annotation.lookup.withExtraMito.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


# Run splitp on FASTQ to preprocess oligo-dT/random hexamers

#!/bin/tcsh

#BSUB -n 1
#BSUB -W 5760
#BSUB -R "rusage[mem=128GB]"
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_1_S66_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_1_S66_L003_R2_001_bcSharing_splitp.fastq.gz
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_2_S67_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_2_S67_L003_R2_001_bcSharing_splitp.fastq.gz
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_3_S68_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_3_S68_L003_R2_001_bcSharing_splitp.fastq.gz
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_4_S69_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_4_S69_L003_R2_001_bcSharing_splitp.fastq.gz
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_5_S70_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_5_S70_L003_R2_001_bcSharing_splitp.fastq.gz
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_6_S71_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_6_S71_L003_R2_001_bcSharing_splitp.fastq.gz
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_7_S72_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_7_S72_L003_R2_001_bcSharing_splitp.fastq.gz
/usr/local/usrapps/Rust/bin/splitp -r /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_8_S73_L003_R2_001.fastq.gz -b /Documents/11Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_8_S73_L003_R2_001_bcSharing_splitp.fastq.gz
# Run Alevin/alevin-fry
# Use --sketch as this will enable use of all exonic and intronic reads
# Running in U/S/A mode

#!/bin/tcsh

#BSUB -n 12
#BSUB -W 5760
#BSUB -R "rusage[mem=128GB]"

#Activate salmon conda environment
conda activate /usr/local/usrapps/salmon/conda-for-salmon
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_1_S66_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_1_S66_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici --rad --sketch
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_2_S67_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_2_S67_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub2_splitp_splici --rad --sketch
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_3_S68_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_3_S68_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub3_splitp_splici --rad --sketch
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_4_S69_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_4_S69_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub4_splitp_splici --rad --sketch
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_5_S70_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_5_S70_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub5_splitp_splici --rad --sketch
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_6_S71_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_6_S71_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub6_splitp_splici --rad --sketch
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_7_S72_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_7_S72_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub7_splitp_splici --rad --sketch
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R1_8_S73_L003_R1_001.fastq.gz -1 /Documents/11Weeks/FastqFiles_NVS079A_Keung_L1L2_R2_8_S73_L003_R2_001_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/11Weeks/11Week_UBE3A_alevinfry_sub8_splitp_splici --rad --sketch

#Deactivate salmon conda environment
conda deactivate
#Activate alevin-fry conda environment
conda activate /usr/local/usrapps/alevin-fry/conda-for-alevin-fry
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub1_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub1_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub1_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub2_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub2_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub2_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub2_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub2_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub2_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub3_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub3_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub3_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub3_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub3_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub3_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub4_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub4_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub4_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub4_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub4_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub4_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub5_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub5_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub5_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub5_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub5_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub5_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub6_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub6_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub6_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub6_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub6_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub6_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub7_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub7_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub7_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub7_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub7_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub7_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
alevin-fry generate-permit-list -d both -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub8_splitp_splici --output-dir /Documents/11Weeks/11Week_UBE3A_alevinfry_sub8_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/11Weeks/11Week_UBE3A_alevinfry_sub8_splitp_splici -t 16 -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub8_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/11Weeks/11Week_UBE3A_alevinfry_sub8_splitp_splici_permit_knee -o 11Week_UBE3A_alevinfry_sub8_splitp_splici_counts -t 16 -r cr-like-em --use-mtx

#Deactivate alevin-fry conda environment
conda deactivate

# Move to R for analysis

library(data.table)
library(presto)
library(fishpond)
library(Seurat)
library(tidyverse)
library(Matrix)
library(glmGamPoi)
library(miQC)
library(SeuratWrappers)
library(flexmix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(sctransform)
library(stringr)
library(stringdist)
library(DESeq2)
library(SeuratData)


setwd("/Documents/11Weeks/")

# Read in entire output directory from alevin-fry using fishpond, will create SingleCellExperiment object
# Fishpond requires --use-mtx, otherwise mtx won't be included in output directory
sce1 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub1_splitp_splici_counts", outputFormat = "snRNA")
sce2 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub2_splitp_splici_counts", outputFormat = "snRNA")
sce3 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub3_splitp_splici_counts", outputFormat = "snRNA")
sce4 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub4_splitp_splici_counts", outputFormat = "snRNA")
sce5 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub5_splitp_splici_counts", outputFormat = "snRNA")
sce6 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub6_splitp_splici_counts", outputFormat = "snRNA")
sce7 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub7_splitp_splici_counts", outputFormat = "snRNA")
sce8 <- loadFry(fryDir="11Week_UBE3A_alevinfry_sub8_splitp_splici_counts", outputFormat = "snRNA")

# Collapse transcripts (gencode.v39) to gene symbols
tx2gene <- read.table("/Documents/genome_annotation/gencode.v39.annotation.lookup.withExtraMito.txt",header=F,sep="\t",col.names=c("tx","gene"))
exp.txId1 <- rownames(counts(sce1))
exp.geneId1 <- as.vector(tx2gene$gene[match(exp.txId1, tx2gene$tx)])
exp.tx.grp1 <- t(sparse.model.matrix(~ 0 + exp.geneId1))
exp.summarized1 <- exp.tx.grp1 %*% counts(sce1)
rownames(exp.summarized1) <- rownames(exp.summarized1) %>% str_replace_all(".+.geneId1","")

exp.txId2 <- rownames(counts(sce2))
exp.geneId2 <- as.vector(tx2gene$gene[match(exp.txId2, tx2gene$tx)])
exp.tx.grp2 <- t(sparse.model.matrix(~ 0 + exp.geneId2))
exp.summarized2 <- exp.tx.grp2 %*% counts(sce2)
rownames(exp.summarized2) <- rownames(exp.summarized2) %>% str_replace_all(".+.geneId2","")

exp.txId3 <- rownames(counts(sce3))
exp.geneId3 <- as.vector(tx2gene$gene[match(exp.txId3, tx2gene$tx)])
exp.tx.grp3 <- t(sparse.model.matrix(~ 0 + exp.geneId3))
exp.summarized3 <- exp.tx.grp3 %*% counts(sce3)
rownames(exp.summarized3) <- rownames(exp.summarized3) %>% str_replace_all(".+.geneId3","")

exp.txId4 <- rownames(counts(sce4))
exp.geneId4 <- as.vector(tx2gene$gene[match(exp.txId4, tx2gene$tx)])
exp.tx.grp4 <- t(sparse.model.matrix(~ 0 + exp.geneId4))
exp.summarized4 <- exp.tx.grp4 %*% counts(sce4)
rownames(exp.summarized4) <- rownames(exp.summarized4) %>% str_replace_all(".+.geneId4","")

exp.txId5 <- rownames(counts(sce5))
exp.geneId5 <- as.vector(tx2gene$gene[match(exp.txId5, tx2gene$tx)])
exp.tx.grp5 <- t(sparse.model.matrix(~ 0 + exp.geneId5))
exp.summarized5 <- exp.tx.grp5 %*% counts(sce5)
rownames(exp.summarized5) <- rownames(exp.summarized5) %>% str_replace_all(".+.geneId5","")

exp.txId6 <- rownames(counts(sce6))
exp.geneId6 <- as.vector(tx2gene$gene[match(exp.txId6, tx2gene$tx)])
exp.tx.grp6 <- t(sparse.model.matrix(~ 0 + exp.geneId6))
exp.summarized6 <- exp.tx.grp6 %*% counts(sce6)
rownames(exp.summarized6) <- rownames(exp.summarized6) %>% str_replace_all(".+.geneId6","")

exp.txId7 <- rownames(counts(sce7))
exp.geneId7 <- as.vector(tx2gene$gene[match(exp.txId7, tx2gene$tx)])
exp.tx.grp7 <- t(sparse.model.matrix(~ 0 + exp.geneId7))
exp.summarized7 <- exp.tx.grp7 %*% counts(sce7)
rownames(exp.summarized7) <- rownames(exp.summarized7) %>% str_replace_all(".+.geneId7","")

exp.txId8 <- rownames(counts(sce8))
exp.geneId8 <- as.vector(tx2gene$gene[match(exp.txId8, tx2gene$tx)])
exp.tx.grp8 <- t(sparse.model.matrix(~ 0 + exp.geneId8))
exp.summarized8 <- exp.tx.grp8 %*% counts(sce8)
rownames(exp.summarized8) <- rownames(exp.summarized8) %>% str_replace_all(".+.geneId8","")

# Rename cells to include sample number, and divide barcodes
exp.cnames1 = gsub(colnames(exp.summarized1),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized1) = exp.cnames1

exp.cnames2 = gsub(colnames(exp.summarized2),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized2) = exp.cnames2

exp.cnames3 = gsub(colnames(exp.summarized3),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized3) = exp.cnames3

exp.cnames4 = gsub(colnames(exp.summarized4),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized4) = exp.cnames4

exp.cnames5 = gsub(colnames(exp.summarized5),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized5) = exp.cnames5

exp.cnames6 = gsub(colnames(exp.summarized6),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized6) = exp.cnames6

exp.cnames7 = gsub(colnames(exp.summarized7),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized7) = exp.cnames7

exp.cnames8 = gsub(colnames(exp.summarized8),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="11Wks_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized8) = exp.cnames8

# Convert gene-summarized matrix to Seurat object
UBE3A_11_Week_sub1 <- CreateSeuratObject(counts = exp.summarized1)
UBE3A_11_Week_sub2 <- CreateSeuratObject(counts = exp.summarized2)
UBE3A_11_Week_sub3 <- CreateSeuratObject(counts = exp.summarized3)
UBE3A_11_Week_sub4 <- CreateSeuratObject(counts = exp.summarized4)
UBE3A_11_Week_sub5 <- CreateSeuratObject(counts = exp.summarized5)
UBE3A_11_Week_sub6 <- CreateSeuratObject(counts = exp.summarized6)
UBE3A_11_Week_sub7 <- CreateSeuratObject(counts = exp.summarized7)
UBE3A_11_Week_sub8 <- CreateSeuratObject(counts = exp.summarized8)

#combine 8 sublibraries into single Seurat object
#merge all 8 sublibraries, adding sublibrary ID to enure cell names are unique
UBE3A_11_Week <- merge(x = UBE3A_11_Week_sub1,y = list(UBE3A_11_Week_sub2, UBE3A_11_Week_sub3, UBE3A_11_Week_sub4, UBE3A_11_Week_sub5, UBE3A_11_Week_sub6, UBE3A_11_Week_sub7, UBE3A_11_Week_sub8), add.cell.ids = c("Sub1", "Sub2", "Sub3", "Sub4", "Sub5", "Sub6", "Sub7", "Sub8"), project = "11_Week_UBE3A")
# Look up barcodes (last 8, this orientation) to get sample name
# Searches for a hamming distance within 1
# If no match, cell name will be prepended with "TRASH" for subsequent filtering
# Matches will have cell name prepended with GTrep-ampBC

bclookup=read.table("bcLookup.txt",header=F,sep="\t")

ids = colnames(UBE3A_11_Week@assays$RNA)
ids_l8 = str_sub(ids, start = -8)

newNames = ids
for(i in 1:length(ids_l8)) {
  ham = stringdist(ids_l8[i],bclookup$V1,method="hamming")
  ct = length(ham[ham<2])
  if(ct==0){
    newNames[i] = paste0("TRASH_",newNames[i])
  }
  else if(ct > 1){
    print(paste0("Barcode ",ids_l8[i]," had ",ct," matches within a Hamming distance of 1"))
    newNames[i] = paste0("TRASH_",newNames[i])
  }
  else{
    bestMatch = which.min(ham)
    matchingBC = as.character(bclookup[bestMatch,1])
    newNames[i] = paste0(as.character(bclookup[bestMatch,2]),"_",newNames[i])
  }
}
UBE3A_11_Week = RenameCells(UBE3A_11_Week,new.names = newNames)

# Parse out sample metadata, drop cells without mathching barcodes
cnames = colnames(UBE3A_11_Week)
names = as.character(gsub("(.+)_(11Wks_UBE3A)_(afUSA)_(.+)","\\1",cnames,perl=T))
names = as.character(gsub("TRASH","TRASH-TRASH",names,perl=T))
sampleNames = as.character(gsub("(.+)-(.+)","\\1",names,perl=T))
reps = as.character(gsub(".{2}([[:digit:]])","Rep\\1",sampleNames,perl=T))
gt = as.character(gsub("(.{2})[[:digit:]]","\\1",sampleNames,perl=T))
ampBC = as.character(gsub("(.+)-(.+)","\\2",names,perl=T))
sublib = as.character(gsub("(.+)-(.+)_(.+)","\\3",names,perl=T))

UBE3A_11_Week = AddMetaData(UBE3A_11_Week,metadata=reps,col.name="rep")
UBE3A_11_Week = AddMetaData(UBE3A_11_Week,metadata=gt,col.name="gt")
UBE3A_11_Week = AddMetaData(UBE3A_11_Week,metadata=ampBC,col.name="ampBC")
UBE3A_11_Week = AddMetaData(UBE3A_11_Week,metadata=sampleNames,col.name="gtrep")
UBE3A_11_Week = AddMetaData(UBE3A_11_Week,metadata=sublib,col.name="sublibrary")

UBE3A_11_Week = subset(UBE3A_11_Week, subset = gt!="TRASH")


# Compute mitochondrial contamination and filter out low quality cells
UBE3A_11_Week <- PercentageFeatureSet(UBE3A_11_Week, pattern = "^MT-", col.name = "percent.mt")
UBE3A_11_Week.filtered <- subset(UBE3A_11_Week, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & percent.mt <10)

# Set up normalization via scTransform
all.list <- SplitObject(UBE3A_11_Week.filtered, split.by = "gtrep")
all.list <- lapply(X = all.list, FUN = SCTransform, vst.flavor = "v2", return.only.var.genes = F, vars.to.regress = "percent.mt")

# Set up integration, then integrate data
options(future.globals.maxSize=5242880000)
all.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 5000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = all.features)
all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", anchor.features = all.features)
UBE3A_11_Week.integrated <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT")

# Run PCA, then determine how many PCs are informative using elbow plot
UBE3A_11_Week.integrated <- RunPCA(UBE3A_11_Week.integrated, verbose = FALSE, npcs = 300)
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_elbow_plot.pdf")
ElbowPlot(UBE3A_11_Week.integrated, ndims = 300, reduction = "pca")
dev.off()

# Run UMAP, and identify cell clusters
UBE3A_11_Week.integrated <- RunUMAP(UBE3A_11_Week.integrated, dims = 1:50)
UBE3A_11_Week.integrated <- FindNeighbors(UBE3A_11_Week.integrated, dims = 1:50, verbose = FALSE)
UBE3A_11_Week.integrated <- FindClusters(UBE3A_11_Week.integrated, verbose = FALSE, resolution = 1.5, algorithm=2)

# Plot UMAP labeled by clusters
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_umap.pdf")
DimPlot(UBE3A_11_Week.integrated,reduction = "umap", label = TRUE)
DimPlot(UBE3A_11_Week.integrated,reduction = "umap", group.by = "gt")
DimPlot(UBE3A_11_Week.integrated,reduction = "umap", group.by = "gtrep")
dev.off()

#Comparison plots
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_umap_comparison.pdf",width=20)
DimPlot(UBE3A_11_Week.integrated,reduction = "umap", split.by = "gt")
DimPlot(UBE3A_11_Week.integrated,reduction = "umap", split.by = "gtrep")
dev.off()

#Median UMI count
median(UBE3A_11_Week.integrated$nCount_RNA)
12756.79
#Median gene count
median(UBE3A_11_Week.integrated$nFeature_RNA)
4558

#Count number of cells in each genotype
metadatatable <- UBE3A_11_Week.integrated@meta.data %>% as.data.table
metadatatable[, .N, by = "gt"]
   gt N
1: WT 4811
2: KO 1948

#Find differentially expressed genes for each cluster
UBE3A_11_Week.integrated <- PrepSCTFindMarkers(UBE3A_11_Week.integrated, assay = "SCT", verbose = TRUE)

for (i in as.numeric(levels(UBE3A_11_Week.integrated))) {
  markers.pre <- FindMarkers(UBE3A_11_Week.integrated, assay="SCT", slot="scale.data", ident.1=i)
  markers <- markers.pre[markers.pre[,2]>0,]
  write.table(markers,paste0("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_res1.5_Markers_c",i,".txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
  dev.off()
}

#Marker gene dot plot
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_goi_cell_type_markers_dotplot.pdf")
celltype.markers = unique(c("PTPRZ1","SLC1A3","GLI3","GLI2","VIM","NOTCH1","NES","PROM1","CLU", #RG
                     "MKI67","TOP2A","CENPK", #Cell cycle
                     "ST18","NHLH1","SOX4","SOX11", #IP-containing
                     "MYT1L","NRXN1","MAPT","DCX", "STMN2","TUBB3","MAP2","RBFOX3", #Neuron
                     "SLC17A6","GRIA2","SLC17A7", #EN
                     "BCL11B","FRMD4B","SATB2","PLXND1","TBR1","FEZF2", #EN-CL
                     "GAD1","GAD2", #IN
                     "TRPM3","HTR2C","TTR", #ChP
                     "RSPO1","RSPO2","RSPO3","MAF", #RSPO+
                     "COL3A1","COL1A2","COL5A1", #Mes
                     "HSP90AA1","HSP90AB1",  #HSP
                     "CDH1","GRHL2","KRT8", #Epithelial cells
                     "PECAM1","FLT1")) #Endothelial cells
DotPlot(UBE3A_11_Week.integrated,features=celltype.markers,assay="SCT",cols = c("blue","red"))+ RotatedAxis()
dev.off()

#label cluster cell annotations
new.cluster.ids <- c("c0_EN", "c1_Neuron/IP", "c2_ChP", "c3_Neuron","c4_EN-CL", "c5_Neuron",
                     "c6_RG", "c7_Proliferating RG", "c8_Mes","c9_RG","c10_RG",
                     "c11_IN","c12_Neuron","c13_Unknown-HSP","c14_Proliferating RG","c15_EN-CL",
                     "c16_EN","c17_RSPO+","c18_RG","c19_ChP","c20_IN","c21_RG","c22_Neuron-HSP",
                     "c23_EN-CL","c24_RG/IP", "c25_Neuron","c26_EN","c27_Mes","c28_IN",
                     "c29_EN","c30_EN","c31_Neuron","c32_ChP","c33_RG-HSP","c34_EN",
                     "c35_Mes","c36_Proliferating RG/Mes","c37_RG/Mes","c38_Mes","c39_RG",
                     "c40_RG/ChP","c41_IN","c42_Epithelial","c43_IN","c44_Endothelial")
names(new.cluster.ids) <- levels(UBE3A_11_Week.integrated)
UBE3A_11_Week.integrated <- RenameIdents(UBE3A_11_Week.integrated, new.cluster.ids)
UBE3A_11_Week.integrated$ClusterIDs <- Idents(UBE3A_11_Week.integrated)

# Plot UMAP labeled by cluster
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_umap_clusterIDs.pdf")
DimPlot(UBE3A_11_Week.integrated,reduction = "umap", label = TRUE)
dev.off()

#label cell type
new.cell.ids <- c("EN", "Neuron/IP", "ChP", "Neuron","EN-CL", "Neuron",
                     "RG", "Proliferating RG", "Mes","RG","RG",
                     "IN","Neuron","Unknown-HSP","Proliferating RG","EN-CL",
                     "EN","RSPO+","RG","ChP","IN","RG","Neuron-HSP",
                     "EN-CL","RG/IP", "Neuron","EN","Mes","IN",
                     "EN","EN","Neuron","ChP","RG-HSP","EN",
                     "Mes","Proliferating RG/Mes","RG/Mes","Mes","RG",
                     "RG/ChP","IN","Epithelial","IN","Endothelial")
names(new.cell.ids) <- levels(UBE3A_11_Week.integrated)
UBE3A_11_Week.integrated <- RenameIdents(UBE3A_11_Week.integrated, new.cell.ids)
UBE3A_11_Week.integrated$CellType <- Idents(UBE3A_11_Week.integrated)

pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_umap_CellType.pdf")
DimPlot(UBE3A_11_Week.integrated,reduction = "umap", label = TRUE)
dev.off()

#Calculate cell type proportions by replicate
proptable = as.data.frame(table(UBE3A_11_Week.integrated$CellType, UBE3A_11_Week.integrated$gtrep)) %>%
  as_tibble() %>%
  dplyr::rename("Sample" = Var2) %>%
  left_join(as_tibble(as.data.frame(table(UBE3A_11_Week.integrated$gtrep))) %>% dplyr::rename("Sample" = Var1), by="Sample") %>%
  mutate(Proportion = Freq.x / Freq.y) %>%
  dplyr::select(-Freq.x,-Freq.y) %>%
  pivot_wider(names_from=Sample,values_from=Proportion) %>%
  dplyr::rename("CellType" = Var1)
write.table(proptable,"UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_proportions_replicates.txt",row.names=F,quote=FALSE,sep="\t")

#Calculate cell type proportions by genotype
proptablegt = as.data.frame(table(UBE3A_11_Week.integrated$CellType, UBE3A_11_Week.integrated$gt)) %>%
  as_tibble() %>%
  dplyr::rename("Sample" = Var2) %>%
  left_join(as_tibble(as.data.frame(table(UBE3A_11_Week.integrated$gt))) %>% dplyr::rename("Sample" = Var1), by="Sample") %>%
  mutate(Proportion = Freq.x / Freq.y) %>%
  dplyr::select(-Freq.x,-Freq.y) %>%
  pivot_wider(names_from=Sample,values_from=Proportion) %>%
  dplyr::rename("CellType" = Var1)

# Cell type proportions stacked bar plot by genotypes
proptablegt = as.data.frame(proptablegt)
rownames(proptablegt)<-proptablegt$CellType
proptablegt$CellType <- NULL
Genotype <- c()
for(i in 1:length(colnames(proptablegt))) {
  Genotype[[i]] <- rep(colnames(proptablegt[i]),length(rownames(proptablegt)))
}
Genotype <- unlist(Genotype)
CellTypes <- rep(c(rownames(proptablegt)),length(colnames(proptablegt)))
Proportion <- as.vector(as.matrix(proptablegt))
dataforbarplot <- data.frame(Genotype,CellTypes,Proportion)
dataforbarplot$CellTypes <- factor(dataforbarplot$CellTypes, levels=unique(CellTypes))
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_proportions_barplots.pdf")
ggplot(dataforbarplot, aes(fill=CellTypes, y=Proportion, x=Genotype)) +
  geom_bar(position='stack', stat='identity')
dev.off()

# Test cell type proportions with propeller using arcsin squareroot transformation
propeller(clusters=UBE3A_11_Week.integrated$CellType, sample = UBE3A_11_Week.integrated$gtrep, group = UBE3A_11_Week.integrated$gt, transform = "asin")

#Highlight significant cell clusters on UMAP
ENcells <- WhichCells(UBE3A_11_Week.integrated, idents = c("EN"))
NeuronIP <- WhichCells(UBE3A_11_Week.integrated, idents = c("Neuron/IP"))
ChPcells <- WhichCells(UBE3A_11_Week.integrated, idents = c("ChP"))
Neuroncells <- WhichCells(UBE3A_11_Week.integrated, idents = c("Neuron"))
ENCL <- WhichCells(UBE3A_11_Week.integrated, idents = c("EN-CL"))
ProlifRG <- WhichCells(UBE3A_11_Week.integrated, idents = c("Proliferating RG"))
Mes <- WhichCells(UBE3A_11_Week.integrated, idents = c("Mes"))
RGCHP <- WhichCells(UBE3A_11_Week.integrated, idents = c("RG/ChP"))
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_UMAP_significant_cell_types_final.pdf")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(ENcells), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(NeuronIP), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(ChPcells), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(Neuroncells), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(ENCL), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(ProlifRG), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(Mes), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_11_Week.integrated, label=T, group.by="gt", cells.highlight = c(RGCHP), cols.highlight = c("red"), cols= "grey")
dev.off()

# Cell type proportions boxplot
pdf("UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5_proportions_replicates_boxplots.pdf")
proptable %>%
  pivot_longer(cols=!CellType,names_to="Sample",values_to="Proportion") %>%
  tidyr::extract(Sample,c("Genotype","Replicate"),"(^.{2})(.{1})") %>%
  mutate(Genotype = fct_relevel(Genotype,c("WT","KO"))) %>%
  ggplot(aes(fill = Genotype, x = Genotype, y=Proportion)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.2)) +
  facet_wrap(~CellType,scales="free") +
  scale_fill_manual(values=c("#f5940c","#0fb7fa")) +
  xlab("")
dev.off()

#All cells differential gene expression
UBE3A_11_Week.integrated <- SetIdent(UBE3A_11_Week.integrated, value = UBE3A_11_Week.integrated$gt)
markers.DGE.all <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="KO", ident.2="WT", logfc.threshold = 0)
write.table(markers.DGE.all,paste0("UBE3A_11_Week_WT_dKO_Differential_Gene_Expression_allcells.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")

#Cell type specfic differential gene expression
cnames = colnames(UBE3A_11_Week.integrated)
names = as.character(gsub("(.+)_(11Wks_UBE3A)_(afUSA)_(.+)","\\1",cnames,perl=T))
names = as.character(gsub("TRASH","TRASH-TRASH",names,perl=T))
sampleNames = as.character(gsub("(.+)-(.+)","\\1",names,perl=T))
gt = as.character(gsub("(.{2})[[:digit:]]","\\1",sampleNames,perl=T))
CellType.gt = paste0(as.character(UBE3A_11_Week.integrated$CellType),sep = "-",gt)
UBE3A_11_Week.integrated = AddMetaData(UBE3A_11_Week.integrated,metadata=CellType.gt,col.name="CellType.gt")
Idents(UBE3A_11_Week.integrated) <- UBE3A_11_Week.integrated$CellType.gt

markers.DGE0 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="RG-KO", ident.2="RG-WT")
write.table(markers.DGE0,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_RG.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE1 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="EN-KO", ident.2="EN-WT")
write.table(markers.DGE1,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_EN.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE2 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Proliferating RG-KO", ident.2="Proliferating RG-WT")
write.table(markers.DGE2,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_ProliferatingRG.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE3 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Neuron-KO", ident.2="Neuron-WT")
write.table(markers.DGE3,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Neuron.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE4 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Mes-KO", ident.2="Mes-WT")
write.table(markers.DGE4,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Mes.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE5 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="ChP-KO", ident.2="ChP-WT")
write.table(markers.DGE5,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_ChP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE6 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Endothelial-KO", ident.2="Endothelial-WT") #Cell group 1 has fewer than 3 cells
write.table(markers.DGE6,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Endothelial.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE7 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Proliferating RG/Mes-KO", ident.2="Proliferating RG/Mes-WT")
write.table(markers.DGE7,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Proliferating RG_Mes.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE8 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="IN-KO", ident.2="IN-WT")
write.table(markers.DGE8,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_IN.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE9 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="RG/IP-KO", ident.2="RG/IP-WT")
write.table(markers.DGE9,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_RG_IP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE10 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Epithelial-KO", ident.2="Epithelial-WT") #Cell group 2 has fewer than 3 cells
write.table(markers.DGE10,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Epithelial.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE11 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="RSPO+-KO", ident.2="RSPO+-WT") #Cell group 2 has fewer than 3 cells
write.table(markers.DGE11,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_RSPO.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE12 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="RG/ChP-KO", ident.2="RG/ChP-WT")  #Cannot find the following identities in the object: RG/ChP-WT
write.table(markers.DGE12,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_RG_ChP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE13 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="RG/Mes-KO", ident.2="RG/Mes-WT")  #Cannot find the following identities in the object: RG/Mes-WT
write.table(markers.DGE13,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_RG_Mes.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE14 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Neuron/IP-KO", ident.2="Neuron/IP-WT")
write.table(markers.DGE14,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Neuron_IP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE15 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="EN-CL-KO", ident.2="EN-CL-WT") #Cell group 1 has fewer than 3 cells
write.table(markers.DGE15,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_EN_CL.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE16 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Unknown-HSP-KO", ident.2="Unknown-HSP-WT")
write.table(markers.DGE16,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Unknown_HSP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE17 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="Neuron-HSP-KO", ident.2="Neuron-HSP-WT") #Cell group 1 has fewer than 3 cells
write.table(markers.DGE17,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_Neuron_HSP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE18 <- FindMarkers(UBE3A_11_Week.integrated, assay="RNA", slot="data", ident.1="RG-HSP-KO", ident.2="RG-HSP-WT") #Cell group 1 has fewer than 3 cells
write.table(markers.DGE18,paste0("Keung_SPLITseq_WT_dKO_Differential_Gene_Expression_RNA_RG_HSP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")

#Create volcano plot of DEGs
#Convert the rownames to a column
markers.DGE.RNA <- cbind(gene=rownames(markers.DGE.RNA), markers.DGE.RNA)
data=markers.DGE.RNA$delabel <- NA
markers.DGE.RNA$delabel[markers.DGE.RNA$diffexpressed != "Not Significant"] <- markers.DGE.RNA$gene[markers.DGE.RNA$diffexpressed != "Not Significant"]
# Remove rows that contain NA values
markers.DGE.RNA<- markers.DGE.all[complete.cases(markers.DGE.all), ]
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
markers.DGE.RNA$diffexpressed <- "Not Significant"
markers.DGE.RNA$diffexpressed[markers.DGE.RNA$avg_log2FC > 1 & markers.DGE.RNA$p_val_adj < 0.05] <- "Up"
markers.DGE.RNA$diffexpressed[markers.DGE.RNA$avg_log2FC < -1 & markers.DGE.RNA$p_val_adj < 0.05] <- "Down"
pdf("UBE3A_11_Week_WT_dKO_Volcano_Plot_Allcells_allgenes.pdf")
ggplot(data=markers.DGE.RNA, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed)) + #, label=top15upanddownregulated)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), col="black",linetype = "dotdash") +
  geom_hline(yintercept=-log10(0.05), col="black",linetype = "dotdash") +
  xlim(-300,300) +
  ggtitle("6 Weeks Differential Expression")+
  labs(x="log2FoldChange",y="-log10(adjusted p-value)",col="Differential Expression")
dev.off()

#Number of upregulated genes
>length(which(markers.DGE.RNA$diffexpressed=="Up"))
[1] 1035
#Number of downregulated genes
> length(which(markers.DGE.RNA$diffexpressed=="Down"))
[1] 5503
#Number of not significant genes
> length(which(markers.DGE.RNA$diffexpressed=="Not Significant"))
[1] 6607

#Save clustered object
saveRDS(UBE3A_11_Week.integrated,"UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5.rds")

#Heatmap of differentially expressed markers between EN-CL and Non-EN-CL Neurons
#Subset neurons
UBE3A_11_Week.integrated.neuron <- UBE3A_11_Week.integrated[,UBE3A_11_Week.integrated$CellType == "Neuron"|UBE3A_11_Week.integrated$CellType =="EN"|UBE3A_11_Week.integrated$CellType =="IN"|UBE3A_11_Week.integrated$CellType =="EN-CL"|UBE3A_11_Week.integrated$CellType =="Neuron/IP"|UBE3A_11_Week.integrated$CellType =="Neuron-HSP"]
NeuronType = as.character(UBE3A_11_Week.integrated.neuron$CellType)
NeuronType = str_replace_all(NeuronType, c("EN-CL"="Cortical","EN"="Non-EN-CL","Neuron/IP"="Non-EN-CL","Neuron-HSP"="Non-EN-CL","Neuron"="Non-EN-CL","IN"="Non-EN-CL","Cortical"="EN-CL"))
UBE3A_11_Week.integrated.neuron = AddMetaData(UBE3A_11_Week.integrated.neuron,metadata=NeuronType,col.name="EN.CL")
UBE3A_11_Week.integrated.neuron <- SetIdent(UBE3A_11_Week.integrated.neuron, value = UBE3A_11_Week.integrated.neuron$EN.CL)
ENCL.DGE <- FindMarkers(UBE3A_11_Week.integrated.neuron, assay="SCT", slot="scale.data", ident.1="EN-CL", ident.2="Non-EN-CL",recorrect_umi = FALSE)
write.table(ENCL.DGE,paste0("UBE3A_11_Week_ECNL_Neuron_Differential_Marker_Expression.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
levels(UBE3A_11_Week.integrated.neuron) <- c("EN-CL","Non-EN-CL")
ENCLNeuron.Markers= unique(c("SLC24A2","LRRC7","ANKS1B","PLXND1","SLC4A10","LRRTM4","KCNQ5","GRIA3","GRIA1","DLGAP1","CACNA1E","ABHD6","RIMS2", #synapse-related
                             "LIN7A","EPHA5","SDK1","CAP2","CTTNBP2","FLRT2", #synapse and axon/dendrite development
                             "NELL2","NFIB","MYO5B","ARPP21","EPHA7", #axon/dendrite development
                             "LIN7A","NEO1", "SHTN1","SRGAP1","SOX5","DAB1",#axon/dendrite development and cell migration/cortical layering
                             "AFF3","BCL11A","ZBTB18",# cell migration/cortical layering
                             #neuronal excitability/neurotransmission
                             "ST3GAL1","NBEA","COPG2", #golgi-related
                             "PDE1A","TG","GAREM1", #signaling-related
                             "CHD7","ZFHX4", #cell differentiation
                             "PDGFC","ADAMTS3","MLLT3","FRMD4B","DACH1","PBX3")) #other
pdf("UBE3A_11_Week_ECNL_Neuron_Differential_Marker_Expression_SCT_res1.5_ENCLDGE_heatmap.pdf",width=10,height=10)
DoHeatmap(UBE3A_11_Week.integrated.neuron, features = ENCLNeuron.Markers, group.by = "EN.CL")
dev.off()

#Save clustered object
saveRDS(UBE3A_11_Week.integrated,"UBE3A_11_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res1.5.rds")

