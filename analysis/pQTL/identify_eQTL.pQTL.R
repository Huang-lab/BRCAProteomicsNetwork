##### find_pQTL.R #####
# Kuan-lin Huang
# run pQTL analysis for 3 cancer types

##### dependencies #####
baseD = "/Users/khuang/Box\ Sync/Huang_lab/Huang_lab_data/CPTAC2retrospective/"
setwd("/Users/khuang/Box\ Sync/Huang_lab/manuscripts/proteomicsNetwork/analysis/pQTL")
source("pQTL.R")

library(methods)
library(reshape2)
system("mkdir out")

driver_gene = read.table(header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/Driver_BaileyCell2018/driver_gene_by_cancer.txt",fill=T)
driver_gene_BRCA = driver_gene$Gene[driver_gene$Cancer=="BRCA"]
driver_gene_OV = driver_gene$Gene[driver_gene$Cancer=="OV"]
driver_gene_CRC = driver_gene$Gene[driver_gene$Cancer=="COADREAD"]

# germline
fn = "/Users/khuang/Box\ Sync/Huang_lab/manuscripts/germlineEthnicPower/data/PCA_pathVar_integrated_filtered_adjusted_ancestry.tsv"
pathVar = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)
pathVarP = pathVar[pathVar$Overall_Classification %in% c("Pathogenic","Likely Pathogenic"),]
pathVar_matrix = dcast(pathVarP,bcr_patient_barcode ~ HUGO_Symbol, value.var = "Variant_Classification")
row.names(pathVar_matrix) = paste(gsub("-",".",gsub("TCGA-","",pathVar_matrix$bcr_patient_barcode)),"01A",sep=".")
pathVar_matrix = pathVar_matrix[,-c(which(colnames(pathVar_matrix)=="bcr_patient_barcode")),]  
pathVar_matrix[pathVar_matrix==0] = "wt"

##### BRCA #####
### Mutation matrix ###
BRCA_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"BRCA/BRCA_MC3_SOMATIC_formatted.txt",sep=""))
BRCA_mut_genes = names(rowSums(BRCA_mut!="wt" & BRCA_mut!="intronic" & BRCA_mut!="silent")[rowSums(BRCA_mut!="wt" & BRCA_mut!="intronic" & BRCA_mut!="silent") > 2])
BRCA_mut_g = BRCA_mut[(row.names(BRCA_mut) %in% BRCA_mut_genes),]

### Transcriptome ###
BRCA_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"BRCA/BRCA_mRNA_formatted_normalized.txt",sep=""))
BRCA_RNA_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_RNA,name="BRCA_mRNA_Somatic")
system("cat out/*BRCA_mRNA_Somatic*txt > out/BRCA_mRNA_Somatic_pQTL.txt")
system("rm -f out/*_BRCA_mRNA_Somatic*txt")
system("gzip out/BRCA_mRNA_Somatic_pQTL.txt")

### Proteome ###
BRCA_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"BRCA/BRCA_PRO_formatted_normalized.txt",sep=""))
BRCA_Pro_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Pro,name="BRCA_Proteome_Somatic")
system("cat out/*BRCA_Proteome_Somatic*txt > out/BRCA_Proteome_Somatic_pQTL.txt")
system("rm -f out/*_BRCA_Proteome_Somatic*txt")
system("gzip out/BRCA_Proteome_Somatic_pQTL.txt")

# germline analyses
BRCA_germ = t(pathVar_matrix[row.names(pathVar_matrix) %in% colnames(BRCA_Pro),])
BRCA_empty_germ = data.frame(matrix(ncol = sum(!colnames(BRCA_Pro) %in% row.names(pathVar_matrix)), nrow = nrow(BRCA_germ)))
row.names(BRCA_empty_germ) = row.names(BRCA_germ)
colnames(BRCA_empty_germ) = colnames(BRCA_Pro)[!colnames(BRCA_Pro) %in% row.names(pathVar_matrix)]
BRCA_empty_germ[is.na(BRCA_empty_germ)] = "wt"
BRCA_germ = cbind(BRCA_germ,BRCA_empty_germ)# need to fill the non-observed sample with NA
BRCA_germ_genes = names(rowSums(BRCA_germ!="wt")[rowSums(BRCA_germ!="wt") > 2])
BRCA_germ_g = BRCA_germ[(row.names(BRCA_germ) %in% BRCA_germ_genes),,drop=F]
BRCA_Pro_germ_diff_exp = find_diff_exp(BRCA_germ_g,BRCA_Pro,name="BRCA_Proteome_Germline")
system("cat out/*BRCA_Proteome_Germline*txt > out/BRCA_Proteome_Germline_pQTL.txt")
system("rm -f out/*_BRCA_Proteome_Germline*txt")
system("gzip out/BRCA_Proteome_Germline_pQTL.txt")

### Phosphoproteome ###
BRCA_Phosphopro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"BRCA/BRCA_PHO_formatted_normalized.txt",sep=""))
BRCA_Phosphopro_diff_exp = find_diff_exp(BRCA_mut_g,BRCA_Phosphopro,name="BRCA_Phosphoproteome_Somatic")
system("cat out/*BRCA_Phosphoproteome_Somatic*txt > out/BRCA_Phosphoproteome_Somatic_pQTL.txt")
system("rm -f out/*_BRCA_Phosphoproteome_Somatic*txt")
system("gzip out/BRCA_Phosphoproteome_Somatic_pQTL.txt")

##### OV #####
### Mutation matrix ###
OV_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"OV/OV_MC3_SOMATIC_formatted.txt",sep=""))
OV_mut_genes = names(rowSums(OV_mut!="wt" & OV_mut!="intronic" & OV_mut!="silent")[rowSums(OV_mut!="wt" & OV_mut!="intronic" & OV_mut!="silent") > 2])
OV_mut_g = OV_mut[(row.names(OV_mut) %in% OV_mut_genes),]

### Transcriptome ###
OV_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"OV/OV_mRNA_formatted_normalized.txt",sep=""))
OV_RNA_diff_exp = find_diff_exp(OV_mut_g,OV_RNA,name="OV_mRNA_Somatic")
system("cat out/*OV_mRNA_Somatic*txt > out/OV_mRNA_Somatic_pQTL.txt")
system("rm -f out/*_OV_mRNA_Somatic*txt")
system("gzip out/OV_mRNA_Somatic_pQTL.txt")

### Proteome ###
OV_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"OV/OV_PNNL_PRO_formatted_normalized.txt",sep=""))
OV_Pro_diff_exp = find_diff_exp(OV_mut_g,OV_Pro,name="OV_Proteome_Somatic")
system("cat out/*OV_Proteome_Somatic*txt > out/OV_Proteome_Somatic_pQTL.txt")
system("rm -f out/*_OV_Proteome_Somatic*txt")
system("gzip out/OV_Proteome_Somatic_pQTL.txt")

# germline analyses
OV_germ = t(pathVar_matrix[row.names(pathVar_matrix) %in% colnames(OV_Pro),])
OV_empty_germ = data.frame(matrix(ncol = sum(!colnames(OV_Pro) %in% row.names(pathVar_matrix)), nrow = nrow(OV_germ)))
row.names(OV_empty_germ) = row.names(OV_germ)
colnames(OV_empty_germ) = colnames(OV_Pro)[!colnames(OV_Pro) %in% row.names(pathVar_matrix)]
OV_empty_germ[is.na(OV_empty_germ)] = "wt"
OV_germ = cbind(OV_germ,OV_empty_germ)# need to fill the non-observed sample with NA
OV_germ_genes = names(rowSums(OV_germ!="wt")[rowSums(OV_germ!="wt") > 2])
OV_germ_g = OV_germ[(row.names(OV_germ) %in% OV_germ_genes),,drop=F]
OV_Pro_germ_diff_exp = find_diff_exp(OV_germ_g,OV_Pro,name="OV_Proteome_Germline")
system("cat out/*OV_Proteome_Germline*txt > out/OV_Proteome_Germline_pQTL.txt")
system("rm -f out/*_OV_Proteome_Germline*txt")
system("gzip out/OV_Proteome_Germline_pQTL.txt")

### Phosphoproteome ###
OV_Phosphopro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"OV/OV_PNNL_PHO_formatted_normalized.txt",sep=""))
OV_Phosphopro_diff_exp = find_diff_exp(OV_mut_g,OV_Phosphopro,name="OV_Phosphoproteome_Somatic")
system("cat out/*OV_Phosphoproteome_Somatic*txt > out/OV_Phosphoproteome_Somatic_pQTL.txt")
system("rm -f out/*_OV_Phosphoproteome_Somatic*txt")
system("gzip out/OV_Phosphoproteome_Somatic_pQTL.txt")

##### CRC #####
### Mutation matrix ###
CRC_mut = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"CRC/CRC_MC3_SOMATIC_formatted.txt",sep=""))
CRC_mut_genes = names(rowSums(CRC_mut!="wt" & CRC_mut!="intronic" & CRC_mut!="silent")[rowSums(CRC_mut!="wt" & CRC_mut!="intronic" & CRC_mut!="silent") > 2])
CRC_mut_g = CRC_mut[(row.names(CRC_mut) %in% CRC_mut_genes),]

### Transcriptome ###
CRC_RNA = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"CRC/CRC_mRNA_formatted_normalized.txt",sep=""))
CRC_RNA_diff_exp = find_diff_exp(CRC_mut_g,CRC_RNA,name="CRC_mRNA_Somatic")
system("cat out/*CRC_mRNA_Somatic*txt > out/CRC_mRNA_Somatic_pQTL.txt")
system("rm -f out/*_CRC_mRNA_Somatic*txt")
system("gzip out/CRC_mRNA_Somatic_pQTL.txt")

### Proteome ###
CRC_Pro = read.table(row.names=1, header=TRUE, sep="\t", file=paste(baseD,"CRC/CRC_PRO_formatted_normalized.txt",sep=""))
CRC_Pro_diff_exp = find_diff_exp(CRC_mut_g,CRC_Pro,name="CRC_Proteome_Somatic")
system("cat out/*CRC_Proteome_Somatic*txt > out/CRC_Proteome_Somatic_pQTL.txt")
system("rm -f out/*_CRC_Proteome_Somatic*txt")
system("gzip out/CRC_Proteome_Somatic_pQTL.txt")

# germline analyses
CRC_germ = t(pathVar_matrix[row.names(pathVar_matrix) %in% colnames(CRC_Pro),])
CRC_empty_germ = data.frame(matrix(ncol = sum(!colnames(CRC_Pro) %in% row.names(pathVar_matrix)), nrow = nrow(CRC_germ)))
row.names(CRC_empty_germ) = row.names(CRC_germ)
colnames(CRC_empty_germ) = colnames(CRC_Pro)[!colnames(CRC_Pro) %in% row.names(pathVar_matrix)]
CRC_empty_germ[is.na(CRC_empty_germ)] = "wt"
CRC_germ = cbind(CRC_germ,CRC_empty_germ)# need to fill the non-observed sample with NA
CRC_germ_genes = names(rowSums(CRC_germ!="wt")[rowSums(CRC_germ!="wt") > 2])
CRC_germ_g = CRC_germ[(row.names(CRC_germ) %in% CRC_germ_genes),,drop=F]
CRC_Pro_germ_diff_exp = find_diff_exp(CRC_germ_g,CRC_Pro,name="CRC_Proteome_Germline")
system("cat out/*CRC_Proteome_Germline*txt > out/CRC_Proteome_Germline_pQTL.txt")
system("rm -f out/*_CRC_Proteome_Germline*txt")
system("gzip out/CRC_Proteome_Germline_pQTL.txt")
