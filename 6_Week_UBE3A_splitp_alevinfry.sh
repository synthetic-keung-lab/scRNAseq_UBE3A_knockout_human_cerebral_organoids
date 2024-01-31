# Build splici index for USA mode quantification
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
/usr/local/usrapps/Rust/bin/splitp -r /Documents/6Weeks/2-LCC1896-LCD1244_R2.fastq.gz -b /Documents/6Weeks/bcSharing.txt -s 79 -e 86 -o | gzip  > /Documents/6Weeks/2-LCC1896-LCD1244_R2_bcSharing_splitp.fastq.gz

# Run Alevin/alevin-fry
# Use --sketch as this will enable use of all exonic and intronic reads
# Running in U/S/A mode

#!/bin/tcsh

#BSUB -n 12
#BSUB -W 5760
#BSUB -R "rusage[mem=128GB]"

#Activate salmon conda environment
conda activate /usr/local/usrapps/salmon/conda-for-salmon
salmon alevin -i /Documents/genome_annotation/gencode_v39_transcriptome_splici_fl95_idx -l A -2 /Documents/6Weeks/2-LCC1896-LCD1244_R1.fastq.gz -1 /Documents/6Weeks/2-LCC1896-LCD1244_R2_bcSharing_splitp.fastq.gz --read-geometry "2[1-end]" --umi-geometry "1[1-10]" --bc-geometry "1[11-18,49-56,79-86]" -p 32 -o /Documents/6Weeks/6Week_UBE3A_alevinfry_splitp_splici --rad --sketch
#Deactivate salmon conda environment
conda deactivate
#Activate alevin-fry conda environment
conda activate /usr/local/usrapps/alevin-fry/conda-for-alevin-fry
alevin-fry generate-permit-list -d both -i /Documents/6Weeks/6Week_UBE3A_alevinfry_splitp_splici --output-dir /Documents/6Weeks/6Week_UBE3A_alevinfry_splitp_splici_permit_knee -k
alevin-fry collate -r /Documents/6Weeks/6Week_UBE3A_alevinfry_splitp_splici -t 16 -i /Documents/6Weeks/6Week_UBE3A_alevinfry_splitp_splici_permit_knee
alevin-fry quant -m /Documents/genome_annotation/encode_v39_transcriptome_splici_fl95_t2g_3col.tsv -i /Documents/6Weeks/6Week_UBE3A_alevinfry_splitp_splici_permit_knee -o 6Week_UBE3A_alevinfry_splitp_splici_counts -t 16 -r cr-like-em --use-mtx
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

setwd("/Documents/6Weeks/")

# Read in entire output directory from alevin-fry using fishpond, will create SingleCellExperiment object
# Fishpond requires --use-mtx, otherwise mtx won't be included in output directory
sce <- loadFry(fryDir="6Week_UBE3A_alevinfry_splitp_splici_counts", outputFormat = "snRNA")

# Collapse transcripts (gencode.v39) to gene symbols
tx2gene <- read.table("/Documents/genome_annotation/gencode.v39.annotation.lookup.withExtraMito.txt",header=F,sep="\t",col.names=c("tx","gene"))
exp.txId <- rownames(counts(sce))
exp.geneId <- as.vector(tx2gene$gene[match(exp.txId, tx2gene$tx)])
exp.tx.grp <- t(sparse.model.matrix(~ 0 + exp.geneId))
exp.summarized <- exp.tx.grp %*% counts(sce)
rownames(exp.summarized) <- rownames(exp.summarized) %>% str_replace_all(".+.geneId","")

# Rename cells to include sample number, and divide barcodes
exp.cnames = gsub(colnames(exp.summarized),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="6Week_UBE3A_afUSA_\\1-\\2-\\3",perl=T)
colnames(exp.summarized) = exp.cnames

# Convert gene-summarized matrix to Seurat object
UBE3A_6_Week <- CreateSeuratObject(counts = exp.summarized)

# Look up barcodes (last 8, this orientation) to get sample name
# Searches for a hamming distance within 1
# If no match, cell name will be prepended with "TRASH" for subsequent filtering
# Matches will have cell name prepended with GTrep-ampBC

bclookup=read.table("bcLookup.txt",header=F,sep="\t")

ids = colnames(UBE3A_6_Week@assays$RNA)
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
UBE3A_6_Week = RenameCells(UBE3A_6_Week,new.names = newNames)

# Parse out sample metadata, drop cells without mathching barcodes
cnames = colnames(UBE3A_6_Week)
names = as.character(gsub("(.+)_(6Week_UBE3A)_(afUSA)_(.+)","\\1",cnames,perl=T))
names = as.character(gsub("TRASH","TRASH-TRASH",names,perl=T))
sampleNames = as.character(gsub("(.+)-(.+)","\\1",names,perl=T))
reps = as.character(gsub(".{2}([[:digit:]])","Rep\\1",sampleNames,perl=T))
gt = as.character(gsub("(.{2})[[:digit:]]","\\1",sampleNames,perl=T))
ampBC = as.character(gsub("(.+)-(.+)","\\2",names,perl=T))

UBE3A_6_Week = AddMetaData(UBE3A_6_Week,metadata=reps,col.name="rep")
UBE3A_6_Week = AddMetaData(UBE3A_6_Week,metadata=gt,col.name="gt")
UBE3A_6_Week = AddMetaData(UBE3A_6_Week,metadata=ampBC,col.name="ampBC")
UBE3A_6_Week = AddMetaData(UBE3A_6_Week,metadata=sampleNames,col.name="gtrep")

UBE3A_6_Week = subset(UBE3A_6_Week, subset = gt!="TRASH")


# Compute mitochondrial contamination and filter out low quality cells
UBE3A_6_Week <- PercentageFeatureSet(UBE3A_6_Week, pattern = "^MT-", col.name = "percent.mt")
UBE3A_6_Week.filtered <- subset(UBE3A_6_Week, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & percent.mt <10)

# Set up normalization via scTransform
all.list <- SplitObject(UBE3A_6_Week.filtered, split.by = "gtrep")
all.list <- lapply(X = all.list, FUN = SCTransform, vst.flavor = "v2", return.only.var.genes = F, vars.to.regress = "percent.mt")

# Set up integration, then integrate data
options(future.globals.maxSize=5242880000)
all.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 5000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = all.features)
all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", anchor.features = all.features)
UBE3A_6_Week.integrated <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT")

# Run PCA, then determine how many PCs are informative using elbow plot
UBE3A_6_Week.integrated <- RunPCA(UBE3A_6_Week.integrated, verbose = FALSE, npcs = 300)
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_elbow_plot.pdf")
ElbowPlot(UBE3A_6_Week.integrated, ndims = 300, reduction = "pca")
dev.off()

# Run UMAP, and identify cell clusters
UBE3A_6_Week.integrated <- RunUMAP(UBE3A_6_Week.integrated, dims = 1:50)
UBE3A_6_Week.integrated <- FindNeighbors(UBE3A_6_Week.integrated, dims = 1:50, verbose = FALSE)
UBE3A_6_Week.integrated <- FindClusters(UBE3A_6_Week.integrated, verbose = FALSE, resolution = 2, algorithm=2)

# Plot UMAP labeled by clusters
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_umap.pdf")
DimPlot(UBE3A_6_Week.integrated,reduction = "umap", label = TRUE)
DimPlot(UBE3A_6_Week.integrated,reduction = "umap", group.by = "gt")
DimPlot(UBE3A_6_Week.integrated,reduction = "umap", group.by = "gtrep")
dev.off()

#Comparison plots
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_umap_comparison.pdf",width=20)
DimPlot(UBE3A_6_Week.integrated,reduction = "umap", split.by = "gt")
DimPlot(UBE3A_6_Week.integrated,reduction = "umap", split.by = "gtrep")
dev.off()

#Median UMI count
median(UBE3A_6_Week.integrated$nCount_RNA)
8278.414
#Median gene count
median(UBE3A_6_Week.integrated$nFeature_RNA)
3558

#Count number of cells in each genotype
metadatatable <- UBE3A_6_Week.integrated@meta.data %>% as.data.table
metadatatable[, .N, by = "gt"]
   gt N
1: KO 4122
2: WT 3376

#Find differentially expressed genes for each cluster
UBE3A_6_Week.integrated <- PrepSCTFindMarkers(UBE3A_6_Week.integrated, assay = "SCT", verbose = TRUE)

for (i in as.numeric(levels(UBE3A_6_Week.integrated))) {
  markers.pre <- FindMarkers(UBE3A_6_Week.integrated, assay="SCT", slot="scale.data", ident.1=i)
  markers <- markers.pre[markers.pre[,2]>0,]
  write.table(markers,paste0("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_res2_Markers_c",i,".txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
  dev.off()
}

#Marker gene dot plot
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_goi_cell_type_markers_dotplot.pdf")
celltype.markers = unique(c("PTPRZ1","SLC1A3","GLI3","VIM","NES", #RG
                     "MKI67","TOP2A","CENPK", #Cell cycle
                     "ASCL1","DLL1","NHLH1", #IP
                     "ST18","ELAVL4", #IP/Neuron
                     "MYT1L","MAPT","RBFOX3", #Neuron
                     "SLC17A6","GRIA2", #EN
                     "GAD1","GAD2", "DLX6-AS1", #IN
                     "WLS","RSPO1","RSPO2","MAF", #RSPO+
                     "TRPM3","HTR2C","TTR", #ChP
                     "COL3A1","COL1A2","COL5A1")) #Mes
DotPlot(UBE3A_6_Week.integrated,features=celltype.markers,assay="SCT",cols = c("blue","red"))+ RotatedAxis()
dev.off()

#label cluster cell annotations
new.cluster.ids <- c("c0_RG", "c1_IN", "c2_Proliferating RG", "c3_EN","c4_Mes", "c5_ChP",
                     "c6_ChP", "c7_RG", "c8_EN/IP","c9_EN","c10_ChP",
                     "c11_RG","c12_EN","c13_EN","c14_Neuron","c15_Proliferating RG","c16_Proliferating RG",
                     "c17_RSPO+","c18_EN","c19_ChP/RG","c20_Neuron","c21_IN/IP","c22_RG","c23_RSPO+/ChP","c24_Proliferating IP",
                     "c25_RG/IP","c26_Neuron/IP","c27_IN","c28_Proliferating RG","c29_EN","c30_IN","c31_EN",
                     "c32_RSPO+","c33_IN")
names(new.cluster.ids) <- levels(UBE3A_6_Week.integrated)
UBE3A_6_Week.integrated <- RenameIdents(UBE3A_6_Week.integrated, new.cluster.ids)
UBE3A_6_Week.integrated$ClusterIDs <- Idents(UBE3A_6_Week.integrated)

# Plot UMAP labeled by cluster
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_umap_clusterIDs.pdf")
DimPlot(UBE3A_6_Week.integrated,reduction = "umap", label = TRUE)
dev.off()

#label cell type
new.cell.ids <- c("RG", "IN", "Proliferating RG", "EN","Mesenchymal", "ChP",
                     "ChP", "RG", "EN/IP","EN","ChP",
                     "RG","IN","EN","Neuron","Proliferating RG","Proliferating RG",
                     "RSPO+","EN","ChP/RG","Neuron","IN/IP","RG","RSPO+/ChP","Proliferating IP",
                     "RG/IP","Neuron/IP","IN","Proliferating RG","EN","IN","EN",
                     "RSPO+","IN")
names(new.cell.ids) <- levels(UBE3A_6_Week.integrated)
UBE3A_6_Week.integrated <- RenameIdents(UBE3A_6_Week.integrated, new.cell.ids)
UBE3A_6_Week.integrated$CellType <- Idents(UBE3A_6_Week.integrated)

pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_umap_CellType.pdf")
DimPlot(UBE3A_6_Week.integrated,reduction = "umap", label = TRUE)
dev.off()

#Calculate cell type proportions by replicate
proptable = as.data.frame(table(UBE3A_6_Week.integrated$CellType, UBE3A_6_Week.integrated$gtrep)) %>%
  as_tibble() %>%
  dplyr::rename("Sample" = Var2) %>%
  left_join(as_tibble(as.data.frame(table(UBE3A_6_Week.integrated$gtrep))) %>% dplyr::rename("Sample" = Var1), by="Sample") %>%
  mutate(Proportion = Freq.x / Freq.y) %>%
  dplyr::select(-Freq.x,-Freq.y) %>%
  pivot_wider(names_from=Sample,values_from=Proportion) %>%
  dplyr::rename("CellType" = Var1)
write.table(proptable,"UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_proportions_replicates.txt",row.names=F,quote=FALSE,sep="\t")


#Calculate cell type proportions by genotype
proptablegt = as.data.frame(table(UBE3A_6_Week.integrated$CellType, UBE3A_6_Week.integrated$gt)) %>%
  as_tibble() %>%
  dplyr::rename("Sample" = Var2) %>%
  left_join(as_tibble(as.data.frame(table(UBE3A_6_Week.integrated$gt))) %>% dplyr::rename("Sample" = Var1), by="Sample") %>%
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
dataforbarplot <- data.frame(genotypes,CellTypes,Proportion)
dataforbarplot$CellTypes <- factor(dataforbarplot$CellTypes, levels=unique(CellTypes))
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_110823_integrated_clustered_res2_Proportions_barplots.pdf")
ggplot(dataforbarplot, aes(fill=CellTypes, y=Proportion, x=Genotype)) +
  geom_bar(position='stack', stat='identity')
dev.off()

# Test cell type proportions with propeller using arcsin squareroot transformation
propeller(clusters=UBE3A_6_Week.integrated$CellType, sample = UBE3A_6_Week.integrated$gtrep, group = UBE3A_6_Week.integrated$gt, transform = "asin")

#Highlight significant cell clusters on UMAP
ProlifRG <- WhichCells(UBE3A_6_Week.integrated, idents = c("Proliferating RG"))
ProlifIP <- WhichCells(UBE3A_6_Week.integrated, idents = c("Proliferating IP"))
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_UMAP_significant_cell_types.pdf")
DimPlot(UBE3A_6_Week.integrated, label=T, group.by="gt", cells.highlight = c(ProlifRG), cols.highlight = c("red"), cols= "grey")
DimPlot(UBE3A_6_Week.integrated, label=T, group.by="gt", cells.highlight = c(ProlifIP), cols.highlight = c("red"), cols= "grey")
dev.off()

# Cell type proportions boxplot
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_proportions_replicates_boxplots.pdf")
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
UBE3A_6_Week.integrated <- SetIdent(UBE3A_6_Week.integrated, value = UBE3A_6_Week.integrated$gt)
markers.DGE.all <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="KO", ident.2="WT", logfc.threshold = 0)
write.table(markers.DGE.all,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_allcells.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")

#Cell type specfic differential gene expression
cnames = colnames(UBE3A_6_Week.integrated)
names = as.character(gsub("(.+)_(6Week_UBE3A)_(afUSA)_(.+)","\\1",cnames,perl=T))
names = as.character(gsub("TRASH","TRASH-TRASH",names,perl=T))
sampleNames = as.character(gsub("(.+)-(.+)","\\1",names,perl=T))
gt = as.character(gsub("(.{2})[[:digit:]]","\\1",sampleNames,perl=T))
CellType.gt = paste0(as.character(UBE3A_6_Week.integrated$CellType),sep = "-",gt)
UBE3A_6_Week.integrated = AddMetaData(UBE3A_6_Week.integrated,metadata=CellType.gt,col.name="CellType.gt")
Idents(UBE3A_6_Week.integrated) <- UBE3A_6_Week.integrated$CellType.gt

markers.DGE0 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="RG-KO", ident.2="RG-WT")
write.table(markers.DGE0,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_RG.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE1 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="EN-KO", ident.2="EN-WT")
write.table(markers.DGE1,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_EN.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE2 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="Proliferating RG-KO", ident.2="Proliferating RG-WT")
write.table(markers.DGE2,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_ProliferatingRG.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE3 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="Neuron-KO", ident.2="Neuron-WT")
write.table(markers.DGE3,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_Neuron.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE4 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="Mesenchymal-KO", ident.2="Mesenchymal-WT")
write.table(markers.DGE4,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_Mesenchymal.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE5 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="ChP-KO", ident.2="ChP-WT")
write.table(markers.DGE5,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_ChP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE6 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="EN/IP-KO", ident.2="EN/IP-WT")
write.table(markers.DGE6,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_EN_IP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE7 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="Proliferating IP-KO", ident.2="Proliferating IP-WT")
write.table(markers.DGE7,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_Proliferating IP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE8 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="IN-KO", ident.2="IN-WT")
write.table(markers.DGE8,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_IN.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE9 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="ChP/RG-KO", ident.2="ChP/RG-WT")
write.table(markers.DGE9,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_ChP_RG.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE10 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="IN/IP-KO", ident.2="IN/IP-WT")
write.table(markers.DGE10,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_IN_IP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE11 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="RSPO+-KO", ident.2="RSPO+-WT")
write.table(markers.DGE11,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_RSPO.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE12 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="RSPO+/ChP-KO", ident.2="RSPO+/ChP-WT")
write.table(markers.DGE12,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_RSPO_ChP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE13 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="RG/IP-KO", ident.2="RG/IP-WT")
write.table(markers.DGE13,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_RG_IP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")
markers.DGE14 <- FindMarkers(UBE3A_6_Week.integrated, assay="RNA", slot="data", ident.1="Neuron/IP-KO", ident.2="Neuron/IP-WT")
write.table(markers.DGE14,paste0("UBE3A_6_Week_WT_dKO_Differential_Gene_Expression_Neuron_IP.txt"),col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")

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
pdf("UBE3A_6_Week_WT_dKO_Volcano_Plot_Allcells_allgenes.pdf")
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
[1] 2233
#Number of downregulated genes
> length(which(markers.DGE.RNA$diffexpressed=="Down"))
[1] 1973
#Number of not significant genes
> length(which(markers.DGE.RNA$diffexpressed=="Not Significant"))
[1] 7562

#Dotplot of EMX1 expression in RG cells
#Subset RG Cells
UBE3A_6_Week.integrated.RG <- subset(x = UBE3A_6_Week.integrated, idents = c("RG","Proliferating RG"))
Idents(UBE3A_6_Week.integrated.RG) <- UBE3A_6_Week.integrated.RG$CellType.gt
levels(UBE3A_6_Week.integrated.RG) <- rev(c("RG-WT","RG-KO","Proliferating RG-WT","Proliferating RG-KO"))

pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_goi_DotPlot_EMX1.pdf")
DotPlot(UBE3A_6_Week.integrated.RG,features="EMX1",assay="SCT",cols = c("blue","red"), scale.min = 0, scale.max = 5)+ RotatedAxis()
dev.off()

#Violin plot for TERT and TP53
UBE3A_6_Week.integrated <- SetIdent(UBE3A_6_Week.integrated, value = UBE3A_6_Week.integrated$gt)
levels(UBE3A_6_Week.integrated) <- c("WT","KO")
features = c("TERT","TP53")
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_033122_integrated_clustered_res2_vlnplot_TERT_TP53.pdf")
for (i in features) {
  print(VlnPlot(UBE3A_6_Week.integrated, features = i, split.by = "gt", assay="SCT", slot="data"))
}
dev.off()

#Violin plot for GLI3
UBE3A_6_Week.integrated <- SetIdent(UBE3A_6_Week.integrated, value = UBE3A_6_Week.integrated$CellType)
levels(UBE3A_6_Week.integrated) <- c("RG","Proliferating RG", "RG/IP","Proliferating IP", "Neuron/IP", "Neuron", "EN","EN/IP","IN","IN/IP",
                                     "RSPO+","RSPO+/ChP","ChP","ChP/RG","Mes")
pdf("UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2_vlnplot_GLI3.pdf")
VlnPlot(UBE3A_6_Week.integrated, features = "GLI3", split.by = "gt", assay="SCT", slot="data")
dev.off()

#Save clustered object
saveRDS(UBE3A_6_Week.integrated,"UBE3A_6_Week_alevinfry_splitp_splici_counts_SCT_integrated_clustered_res2.rds")

