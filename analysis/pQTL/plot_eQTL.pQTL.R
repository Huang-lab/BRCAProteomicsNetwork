##### plot_pQTL.R #####
# Kuan-lin Huang
# plot pQTL vs eQTL analysis results

##### dependencies #####
setwd("/Users/khuang/Box\ Sync/Huang_lab/manuscripts/proteomicsNetwork/analysis/pQTL")
source("pQTL.R")
source("../global_aes_out.R")

library(methods)
library(reshape2)
system("mkdir out")

# Read in and pre-process pQTL results
### mRNA ###
BRCA_mRNA = read.table(header=TRUE, sep="\t", stringsAsFactors = F, file= gzfile("out/BRCA_mRNA_Somatic_pQTL.txt.gz"))
colnames(BRCA_mRNA)[1] = "Gene"
BRCA_mRNA$data = "RNA"
BRCA_mRNA$cancer = "BRCA"
BRCA_mRNA$MutGene = gsub("_.*","",BRCA_mRNA$name)

OV_mRNA = read.table(header=TRUE, sep="\t", stringsAsFactors = F, file= gzfile("out/OV_mRNA_Somatic_pQTL.txt.gz"))
colnames(OV_mRNA)[1] = "Gene"
OV_mRNA$data = "RNA"
OV_mRNA$cancer = "OV"
OV_mRNA$MutGene = gsub("_.*","",OV_mRNA$name)

CRC_mRNA = read.table(header=TRUE, sep="\t", stringsAsFactors = F, file= gzfile("out/CRC_mRNA_Somatic_pQTL.txt.gz"))
colnames(CRC_mRNA)[1] = "Gene"
CRC_mRNA$data = "RNA"
CRC_mRNA$cancer = "CRC"
CRC_mRNA$MutGene = gsub("_.*","",CRC_mRNA$name)

All_RNA = rbind(BRCA_mRNA,OV_mRNA, CRC_mRNA)
All_RNA = All_RNA[All_RNA$MutGene!="name",]
remove(BRCA_mRNA)
remove(OV_mRNA)
remove(CRC_mRNA)
All_RNA_sele = All_RNA[,c("Gene","meanChange","w_test_fdr","cancer","MutGene")]
colnames(All_RNA_sele)[2:3] = paste("RNA",colnames(All_RNA_sele)[2:3],sep="")

### proteome ###
BRCA_Pro = read.table(header=TRUE, sep="\t", stringsAsFactors = F, file= gzfile("out/BRCA_Proteome_Somatic_pQTL.txt.gz"))
colnames(BRCA_Pro)[1] = "Gene"
BRCA_Pro$data = "Protein"
BRCA_Pro$cancer = "BRCA"
BRCA_Pro$MutGene = gsub("_.*","",BRCA_Pro$name)

OV_Pro = read.table(header=TRUE, sep="\t", stringsAsFactors = F, file= gzfile("out/OV_Proteome_Somatic_pQTL.txt.gz"))
colnames(OV_Pro)[1] = "Gene"
OV_Pro$data = "Protein"
OV_Pro$cancer = "OV"
OV_Pro$MutGene = gsub("_.*","",OV_Pro$name)

CRC_Pro = read.table(header=TRUE, sep="\t", stringsAsFactors = F, file= gzfile("out/CRC_Proteome_Somatic_pQTL.txt.gz"))
colnames(CRC_Pro)[1] = "Gene"
CRC_Pro$data = "Protein"
CRC_Pro$cancer = "CRC"
CRC_Pro$MutGene = gsub("_.*","",CRC_Pro$name)

All_Pro = rbind(BRCA_Pro,OV_Pro, CRC_Pro)
All_Pro = All_Pro[All_Pro$MutGene!="name",]
remove(BRCA_Pro)
remove(OV_Pro)
remove(CRC_Pro)
All_Pro_sele = All_Pro[,c("Gene","meanChange","w_test_fdr","cancer","MutGene")]
colnames(All_Pro_sele)[2:3] = paste("PRO",colnames(All_Pro_sele)[2:3],sep="")

All_QTL_bind = rbind(All_RNA, All_Pro) 
All_QTL_merge = merge(All_RNA_sele, All_Pro_sele, by=c("Gene","cancer","MutGene"),all=F)
for (i in 4:7){
  All_QTL_merge[,i] = as.numeric(All_QTL_merge[,i])
}
for (i in 3:11){
  All_QTL_bind[,i] = as.numeric(All_QTL_bind[,i])
}

p = ggplot(data=All_QTL_merge[All_QTL_merge$RNAw_test_fdr < 0.1 | All_QTL_merge$PROw_test_fdr < 0.1,], aes(y = RNAmeanChange, x = PROmeanChange))
p = p + facet_grid(.~cancer)
p = p + geom_point(alpha=0.2,stroke=0)
p = p + labs(y="RNA Mean Change", x="Protein Mean Change") + theme_bw()
#p = p + geom_text_repel(size =2, aes(label=ifelse(PROw_test_fdr < 0.05 & MutGene != "TP53", paste(MutGene, Gene, sep=":"),NA)))
p = p + coord_fixed()
p
fn = 'out/Mut_RNA_vs_PRO_change_3cancer.pdf'
ggsave(file=fn, useDingbats=FALSE)

fdr.colors=c("NA", "#000000")
markers = All_QTL_merge$Gene[!is.na(All_QTL_merge$PROw_test_fdr) & All_QTL_merge$PROw_test_fdr< 0.05 & All_QTL_merge$MutGene != "TP53"]
p = ggplot(data=All_QTL_bind[All_QTL_bind$Gene %in% markers,])
p = p + facet_grid(data~cancer, drop=T, scales = "free", space = "free")
p = p + geom_tile(aes(x=MutGene, y=Gene, fill=as.numeric(meanChange)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-5,5))
p = p + geom_tile(aes(x=MutGene, y=Gene, color=w_test_fdr<0.05), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR < 0.05"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Marker") + theme_bw() +
  theme(axis.title = element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=6),
        axis.text.y = element_text(colour="black", size=6), axis.ticks = element_blank(),
        strip.text.x = element_text(size = 6))#element_text(colour="black", size=16))
#p = p + theme(legend.position="none")
p
fn = paste('out/Mut_RNA_vs_PRO_heatmap_p53_exclude_pan3can.pdf',sep ="_")
ggsave(file=fn, w=9, h = 4,useDingbats=FALSE)

### Focus on breast cancer ###

# 5% for mRNA and protein data
BRCA_RNA_fdr_cutoff = quantile(as.numeric(BRCA_mRNA$w_test_fdr), na.rm = T, probs = c(0.05))
cat("Top 5% w test FDR for eQTL:",BRCA_RNA_fdr_cutoff,"\n")
BRCA_PRO_fdr_cutoff = quantile(as.numeric(BRCA_Pro$w_test_fdr), na.rm = T, probs = c(0.05))
cat("Top 5% w test FDR for pQTL:",BRCA_PRO_fdr_cutoff,"\n")

All_QTL_bind_BRCA = All_QTL_bind[All_QTL_bind$cancer=="BRCA",]
All_QTL_merge_BRCA = All_QTL_merge[All_QTL_merge$cancer=="BRCA",]
cat("Mean change between RNA and Protein in breast cancer:\n")
cor.test(All_QTL_merge_BRCA$RNAmeanChange,All_QTL_merge_BRCA$PROmeanChange)
All_QTL_merge_BRCA$PROw_test_fdr_logDic = log10(All_QTL_merge_BRCA$PROw_test_fdr)
All_QTL_merge_BRCA$PROw_test_fdr_logDic[!is.na(All_QTL_merge_BRCA$PROmeanChange) & All_QTL_merge_BRCA$PROmeanChange<0] = -All_QTL_merge_BRCA$PROw_test_fdr_logDic[!is.na(All_QTL_merge_BRCA$PROmeanChange) & All_QTL_merge_BRCA$PROmeanChange<0]

All_QTL_merge_BRCA$RNAw_test_fdr_logDic = log10(All_QTL_merge_BRCA$RNAw_test_fdr)
All_QTL_merge_BRCA$RNAw_test_fdr_logDic[All_QTL_merge_BRCA$RNAmeanChange<0] = -All_QTL_merge_BRCA$RNAw_test_fdr_logDic[All_QTL_merge_BRCA$RNAmeanChange<0]
All_QTL_merge_BRCA$category = "None"
All_QTL_merge_BRCA$category[All_QTL_merge_BRCA$RNAw_test_fdr <= BRCA_RNA_fdr_cutoff] = "eQTL"
All_QTL_merge_BRCA$category[All_QTL_merge_BRCA$PROw_test_fdr <= BRCA_PRO_fdr_cutoff] = "pQTL"
All_QTL_merge_BRCA$category[All_QTL_merge_BRCA$RNAw_test_fdr <= BRCA_RNA_fdr_cutoff & All_QTL_merge_BRCA$PROw_test_fdr < BRCA_PRO_fdr_cutoff] = "eQTL_pQTL"

p = ggplot(data=All_QTL_merge_BRCA, aes(y = RNAw_test_fdr_logDic, x = PROw_test_fdr_logDic))
p = p + facet_grid(.~cancer)
p = p + geom_point(aes(color = category),alpha=0.2,stroke=0)
p = p + labs(y="-log10(eQTL FDR)", x="-log10(pQTL FDR)") + theme_bw()
p = p + geom_text_repel(size =2, aes(label=ifelse(category == "pQTL" & PROw_test_fdr < 0.03, paste(MutGene, Gene, sep=":"),NA)))
p = p + coord_fixed()
p
fn = 'out/Mut_RNA_vs_PRO_qq_plot_BRCA.pdf'
ggsave(file=fn, w=5,h =4, useDingbats=FALSE)

fdr.colors=c("NA", "#000000")
markers = All_QTL_merge_BRCA$Gene[!is.na(All_QTL_merge_BRCA$PROw_test_fdr) & All_QTL_merge_BRCA$PROw_test_fdr < 0.3 & All_QTL_merge_BRCA$MutGene != "TP53"]
p = ggplot(data=All_QTL_bind_BRCA[All_QTL_bind_BRCA$Gene %in% markers,])
p = p + facet_grid(Gene~cancer, drop=T, scales = "free", space = "free")
p = p + geom_tile(aes(y=data, x=MutGene, fill=as.numeric(meanChange)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-3,3))
p = p + geom_tile(aes(y=data, x=MutGene, color=w_test_fdr<BRCA_PRO_fdr_cutoff), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR < 0.3"),values = fdr.colors)
p = p + labs(x="Mutated Gene", y="Marker") + theme_bw() #+
  # theme(axis.title = element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=6),
  #       axis.text.y = element_text(colour="black", size=6), axis.ticks = element_blank(),
  #       strip.text.x = element_text(size = 6))#element_text(colour="black", size=16))
#p = p + theme(legend.position="none")
p
fn = paste('out/Mut_RNA_vs_PRO_heatmap_p53_exclude_BRCA.pdf',sep ="_")
ggsave(file=fn,h=4.5, w=7.5,useDingbats=FALSE)

# ##### merge everything and plot #####
# 
# ### fold change
# BRCA_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$fold_change)))
# colnames(BRCA_Pro_diff_exp_fc_m)[3] = "BRCA_PRO"
# BRCA_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$fold_change)))
# colnames(BRCA_Pho_diff_exp_fc_m)[3] = "BRCA_PHO"
# # OV_JHU_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_JHU_Pro_diff_exp$fold_change)))
# # colnames(OV_JHU_Pro_diff_exp_fc_m)[3] = "OV_JHU_PRO"
# OV_PNNL_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pro_diff_exp$fold_change)))
# colnames(OV_PNNL_Pro_diff_exp_fc_m)[3] = "OV_PRO"
# # OV_JHU_Gly_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_JHU_Gly_diff_exp$fold_change)))
# # colnames(OV_JHU_Gly_diff_exp_fc_m)[3] = "OV_JHU_GLY"
# OV_PNNL_Pho_diff_exp_fc_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pho_diff_exp$fold_change)))
# colnames(OV_PNNL_Pho_diff_exp_fc_m)[3] = "OV_PHO"
# # CRC_Pro_diff_exp_fc_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$fold_change)))
# # colnames(CRC_Pro_diff_exp_fc_m)[3] = "CRC_PRO"
# 
# fc = merge(BRCA_Pro_diff_exp_fc_m, BRCA_Pho_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# # fc = merge(fc, OV_JHU_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# fc = merge(fc, OV_PNNL_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# #fc = merge(fc, OV_JHU_Gly_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# fc = merge(fc, OV_PNNL_Pho_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# # fc = merge(fc, CRC_Pro_diff_exp_fc_m, by = c("Var1","Var2"), all=T)
# 
# fc_m = melt(fc, id.var=c("Var1","Var2"))
# colnames(fc_m)[4] = "FC"
# ## fdr
# BRCA_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pro_diff_exp$t_fdr)))
# colnames(BRCA_Pro_diff_exp_fdr_m)[3] = "BRCA_PRO"
# BRCA_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,BRCA_Pho_diff_exp$t_fdr)))
# colnames(BRCA_Pho_diff_exp_fdr_m)[3] = "BRCA_PHO"
# # OV_JHU_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_JHU_Pro_diff_exp$t_fdr)))
# # colnames(OV_JHU_Pro_diff_exp_fdr_m)[3] = "OV_JHU_PRO"
# OV_PNNL_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pro_diff_exp$t_fdr)))
# colnames(OV_PNNL_Pro_diff_exp_fdr_m)[3] = "OV_PRO"
# # OV_JHU_Gly_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_JHU_Gly_diff_exp$t_fdr)))
# # colnames(OV_JHU_Gly_diff_exp_fdr_m)[3] = "OV_JHU_GLY"
# OV_PNNL_Pho_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,OV_PNNL_Pho_diff_exp$t_fdr)))
# colnames(OV_PNNL_Pho_diff_exp_fdr_m)[3] = "OV_PHO"
# # CRC_Pro_diff_exp_fdr_m = melt(as.matrix(do.call(cbind,CRC_Pro_diff_exp$t_fdr)))
# # colnames(CRC_Pro_diff_exp_fdr_m)[3] = "CRC_PRO"
# 
# fdr = merge(BRCA_Pro_diff_exp_fdr_m, BRCA_Pho_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# # fdr = merge(fdr, OV_JHU_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# fdr = merge(fdr, OV_PNNL_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# #fdr = merge(fdr, OV_JHU_Gly_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# fdr = merge(fdr, OV_PNNL_Pho_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# # fdr = merge(fdr, CRC_Pro_diff_exp_fdr_m, by = c("Var1","Var2"), all=T)
# 
# fdr_m = melt(fdr, id.var=c("Var1","Var2"))
# colnames(fdr_m)[4] = "FDR"
# 
# fc_fdr = merge(fc_m, fdr_m, by=c("Var1","Var2","variable"))
# 
# # look for the kinase-substrate relations
# fc_fdr$KSpair = paste(fc_fdr$Var2,fc_fdr$Var1) %in% paste(k_s_table$GENE,k_s_table$SUB_GENE)
# 
# ## plot
# fdr.colors=c("NA", "#000000")
# min_d = min(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
# max_d = max(as.numeric(as.character(fc_fdr$FC)), na.rm=T)
# bound = max(c(max_d, -min_d))
# fc_fdr$sig = as.numeric(as.character(fc_fdr$FDR)) < 0.05
# fc_fdr$FC_2 = as.numeric(fc_fdr$FC)
# fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2>=2,]$FC_2=2
# fc_fdr[!is.na(fc_fdr$FC_2) & fc_fdr$FC_2<=-2,]$FC_2=-2
# fc_fdr.v = fc_fdr[rowSums(is.na(fc_fdr))<4,]
# ## dropping lvls, doesn't work yet
# # fc_fdr.v = droplevels(fc_fdr.v)
# # fc_fdr.v$Var1 = factor(fc_fdr.v$Var1)
# # fc_fdr.v$Var2 = factor(fc_fdr.v$Var2)
# # fc_fdr.v[fc_fdr.v$variable == "BRCA_PRO",]$Var2 = droplevels(fc_fdr.v[fc_fdr.v$variable == "BRCA_PRO",]$Var2) 
# 
# fn = 'out/merged_diff_exp_cgenes_kinase_substrate.pdf'
# p = ggplot(data=fc_fdr.v)
# #p = p + facet_grid(.~variable, scales = "fixed", space = "free", drop=TRUE)
# p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
# p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
# p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR < 0.05"),values = fdr.colors)
# p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=8),axis.ticks = element_blank())#element_text(colour="black", size=16))
# p = p + coord_equal()
# p
# ggsave(file=fn, width=10.5, height=18, useDingbats=FALSE)
# 
# 
# ### plot selected markers ###
# for (i in 3:6){
#   fdr[,i] = as.numeric(as.character(fdr[,i]))
# }
# markers = unique(fdr[rowSums(fdr[,c(3:6)]<=0.05, na.rm=T) >=1,]$Var1)
# fc_fdr_s = fc_fdr.v[as.character(fc_fdr.v$Var1) %in% as.character(markers),]
# #fc_fdr_s2 = fc_fdr_s[fc_fdr_s$FDR<=0.2,]
# 
# fn = 'out/merged_diff_exp_cgenes_kinase_substrate_fdr0.05in1.pdf'
# p = ggplot(data=fc_fdr_s)
# #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
# p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
# p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
# p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR < 0.05"),values = fdr.colors)
# p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
# p = p + coord_fixed()
# p
# ggsave(file=fn, width=10.5, height=6, useDingbats=FALSE)
# 
# fn = 'out/merged_diff_exp_cgenes_kinase_substrate_fdr0.05in1_p.pdf'
# p = ggplot(data=fc_fdr_s)
# #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
# p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
# p = p + geom_point(aes(x=Var2, y=Var1, fill=as.numeric(FC_2), size=-log10(as.numeric(FDR)), color=ifelse(sig, "black",NA)),pch=21) 
# p = p + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
# p = p + scale_colour_manual(values=c("black",NA))
# p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
# p = p + coord_fixed()
# p
# ggsave(file=fn, width=10.5, height=6, useDingbats=FALSE)
# 
# # ### for grant ### 
# # if (FALSE){
# #   markers = c("AR","BRAF","CDH1","CHEK2","EGFR","ERBB2","ESR1","GATA3","PIK3CA","PDGFRA","TP53")
# #   fc_fdr_s = fc_fdr.v[fc_fdr.v$Var1 %in% as.character(markers),]
# #   fc_fdr_s = fc_fdr_s[fc_fdr_s$variable %in% c("BRCA_PRO","BRCA_PHO","OV_PNNL_PRO","OV_PNNL_PHO","CRC_PRO"),]
# #   fc_fdr_s$variable = fc_fdr_s$variable[drop=T]
# #   levels(fc_fdr_s$variable) = c("BRCA_PRO","BRCA_PHO","OV_PRO","OV_PHO","CRC_PRO")
# #   
# #   fn = paste('out/merged_diff_exp_cgenes_druggable_U24_label.pdf',sep ="_")
# #   p = ggplot(data=fc_fdr_s)
# #   #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
# #   p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
# #   p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
# #   p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
# #   p = p + labs(x="Mutated Gene", y="Expression") + theme_bw() + 
# #     theme(axis.title = element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=6), 
# #           axis.text.y = element_text(colour="black", size=6), axis.ticks = element_blank(),
# #           strip.text.x = element_text(size = 6))#element_text(colour="black", size=16))
# #   p = p + coord_equal()
# #   #p = p + theme(legend.position="none")
# #   p
# #   ggsave(file=fn, width=10, height=5, useDingbats=FALSE)
# # }
# 
# # ### for PTRC grant ###
# # if (FALSE){
# #   # exclude TSG
# #   oncogenes_f = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/Volgestin_Oncogene_list.txt', header=FALSE, stringsAsFactors = F)
# #   oncogenes = as.vector(t(oncogenes_f))
# #   fc_fdr_s_onco = fc_fdr_s[fc_fdr_s$Var1 %in% oncogenes & fc_fdr_s$variable != "OV_JHU_PRO",]
# #   fc_fdr_s_onco$variable = gsub("_.*_","_", fc_fdr_s_onco$variable)
# #   fc_fdr_s_onco$variable = factor(fc_fdr_s_onco$variable, levels = c("BRCA_PRO","BRCA_PHO","OV_PRO","OV_PHO","CRC_PRO"))
# #   
# #   fn = paste('out/merged_diff_exp_cgenes_druggable_fdr0.1in1_oncogenes.pdf',sep ="_")
# #   p = ggplot(data=fc_fdr_s_onco)
# #   #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
# #   p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
# #   p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
# #   p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
# #   p = p + labs(x="Mutated Gene", y="Expression") + theme_nogrid() + 
# #     theme(axis.title = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=10), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())#element_text(colour="black", size=16))
# #   p = p + coord_fixed()
# #   p
# #   ggsave(file=fn, width=8.5, height=4, useDingbats=FALSE)
# #   
# #   
# #   fc_fdr_s = fc_fdr.v[fc_fdr.v$Var1 %in% as.character(oncogenes),]
# #   fc_fdr_s = fc_fdr_s[fc_fdr_s$variable %in% c("BRCA_PRO","BRCA_PHO","OV_PNNL_PRO","OV_PNNL_PHO","CRC_PRO"),]
# #   fc_fdr_s$variable = fc_fdr_s$variable[drop=T]
# #   levels(fc_fdr_s$variable) = c("BRCA_PRO","BRCA_PHO","OV_PRO","OV_PHO","CRC_PRO")
# #   
# #   fn = paste('out/merged_diff_exp_cgenes_druggable_oncogenes.pdf',sep ="_")
# #   p = ggplot(data=fc_fdr_s)
# #   #p = p + facet_grid(.~variable, drop=T, scales = "free", space = "free")
# #   p = p + facet_wrap(~variable, nrow=1, scales = "free_x", drop=TRUE)
# #   p = p + geom_tile(aes(x=Var2, y=Var1, fill=as.numeric(FC_2)), linetype="blank") + scale_fill_gradientn(name= "Average expression difference", na.value=NA, colours=RdBu1024, limit=c(-2,2))
# #   p = p + geom_tile(aes(x=Var2, y=Var1, color=sig), fill=NA, size=0.5) + scale_colour_manual(name=paste("FDR <= 0.1"),values = fdr.colors)
# #   p = p + labs(x="Mutated Gene", y="Expression") + theme_nogrid() + 
# #     theme(axis.title = element_text(size=8), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=6), 
# #           axis.text.y = element_text(colour="black", size=6), axis.ticks = element_blank(),
# #           strip.text.x = element_text(size = 6))#element_text(colour="black", size=16))
# #   p = p + coord_equal()
# #   #p = p + theme(legend.position="none")
# #   p
# #   ggsave(file=fn, width=7, height=5, useDingbats=FALSE)
# # }
# # ##### plot collective violin plot #####
# # ### make violin plots for selected mutated gene-marker combo###
# # 
# # plot_diff_exp_violin_c = function (Mgene, gene){
# #   BRCA_mut_gene = data.frame(t(BRCA_mut[Mgene,]))
# #   BRCA_mut_gene$carrier = as.character(BRCA_mut_gene[,1]) != "wt" & as.character(BRCA_mut_gene[,1]) != "silent"
# #   BRCA_Pro_gene = t(BRCA_Pro[gene,])
# #   BRCA_Pho_gene = t(BRCA_Pho[gene,])
# #   # merge by sample
# #   BRCA_Pro_merge = merge(BRCA_mut_gene, BRCA_Pro_gene, by = "row.names")
# #   BRCA_Pro_merge$data = "BRCA_PRO"
# #   BRCA_Pho_merge = merge(BRCA_mut_gene, BRCA_Pho_gene, by = "row.names")
# #   BRCA_Pho_merge$data = "BRCA_PHO"
# #   
# #   CRC_mut_gene = data.frame(t(CRC_mut[Mgene,]))
# #   CRC_mut_gene$carrier = as.character(CRC_mut_gene[,1]) != "wt" & as.character(CRC_mut_gene[,1]) != "silent"
# #   CRC_Pro_gene = t(CRC_Pro[gene,])
# #   CRC_Pro_merge = merge(CRC_mut_gene, CRC_Pro_gene, by = "row.names")
# #   CRC_Pro_merge$data = "CRC_PRO"
# #   
# #   OV_mut_gene = data.frame(t(OV_mut[Mgene,]))
# #   OV_mut_gene$carrier = as.character(OV_mut_gene[,1]) != "wt" & as.character(OV_mut_gene[,1]) != "silent"
# #   OV_JHU_Pro_gene = t(OV_JHU_Pro[gene,])
# #   OV_JHU_Pro_merge = merge(OV_mut_gene, OV_JHU_Pro_gene, by = "row.names")
# #   OV_JHU_Pro_merge$data = "OV_JHU_PRO"
# #   OV_PNNL_Pro_gene = t(OV_PNNL_Pro[gene,])
# #   OV_PNNL_Pro_merge = merge(OV_mut_gene, OV_PNNL_Pro_gene, by = "row.names")
# #   OV_PNNL_Pro_merge$data = "OV_PNNL_PRO"
# #   OV_PNNL_Pho_gene = t(OV_PNNL_Pho[gene,])
# #   OV_PNNL_Pho_merge = merge(OV_mut_gene, OV_PNNL_Pho_gene, by = "row.names")
# #   OV_PNNL_Pho_merge$data = "OV_PNNL_PHO"
# #   
# #   colnames(BRCA_Pro_merge)[4]=gene
# #   colnames(BRCA_Pho_merge)[4]=gene
# #   colnames(CRC_Pro_merge)[4]=gene
# #   colnames(OV_JHU_Pro_merge)[4]=gene
# #   colnames(OV_PNNL_Pro_merge)[4]=gene
# #   colnames(OV_PNNL_Pho_merge)[4]=gene
# #   
# #   gene_all_lvl = rbind(BRCA_Pro_merge,BRCA_Pho_merge,CRC_Pro_merge,OV_JHU_Pro_merge,OV_PNNL_Pro_merge,OV_PNNL_Pho_merge)
# #   colnames(gene_all_lvl) = c("Sample","Mutation_Type","Mutation_Status","Expression","Dataset")
# #   
# #   # plot violin plots faceted by marker genes
# #   fn = paste(pd, Mgene, gene, "mutational_impact_violin.pdf", sep="_")
# #   p = ggplot(data=gene_all_lvl)
# #   p = p + facet_grid(.~Dataset)
# #   p = p + geom_violin(aes(x=Mutation_Status, y=Expression, fill=Mutation_Status),alpha=0.5) + guides(fill=FALSE) 
# #   p = p + geom_jitter(aes(x=Mutation_Status, y=Expression, color=Mutation_Type)) #+ geom_point(aes(x=Status, y=value)) 
# #   p = p + labs(x = paste(Mgene,"Mutation Status"), y = paste(gene, "Expression")) + theme_bw()
# #   p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(colour="black", size=14), 
# #                 axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
# #   p
# #   ggsave(file=fn, width=14.5, limitsize=FALSE, useDingbats=FALSE)
# # }
# # 
# # # plot any pair of interest
# # if (FALSE){
# #   plot_diff_exp_violin_c("TP53","TP53")
# #   plot_diff_exp_violin_c("TP53","ESR1")
# #   plot_diff_exp_violin_c("TP53","IGF1R")
# #   plot_diff_exp_violin_c("TP53","GATA3")
# #   plot_diff_exp_violin_c("TP53","EGFR")
# #   plot_diff_exp_violin_c("TP53","CDH1")
# #   plot_diff_exp_violin_c("GATA3","EGFR")
# #   plot_diff_exp_violin_c("TP53","CHEK2")
# #   plot_diff_exp_violin_c("NF1","TP53")
# #   plot_diff_exp_violin_c("CDH1","CDH1")
# #   plot_diff_exp_violin_c("CDH1","PDGFRB")
# #   plot_diff_exp_violin_c("KRAS","MAP2K1")
# # }
# # 
# # gene="MDM2"
