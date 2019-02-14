##### analyze network results.R #####
# Kuan-lin Huang

##### dependencies #####
setwd("/Users/khuang/Box Sync/Huang_lab/manuscripts/BRCAProteomicsNetwork/analysis/analyze_network")
source("../global_aes_out.R")

library(methods)
library(reshape2)
system("mkdir out")

# Read in and pre-process pQTL results
### mRNA ###
MEGENA = read.table(header=TRUE, sep="\t", stringsAsFactors = F, file= "../../data/multiscale_module_summary.ranked.txt")
MEGENA_QTL_p_cutoff = quantile(MEGENA[,grep("FET.P_.*BRCA",colnames(MEGENA))], na.rm = T, probs = c(0.001)) # top 1%
MEGENA$module_QTL_assoc = rowSums(MEGENA[,grep("FET.P_.*BRCA",colnames(MEGENA))] < MEGENA_QTL_p_cutoff) # arbitrary for now; top 1%
MEGENA$name = paste(MEGENA$id,gsub(",.*","",MEGENA$module.hub))

MEGENA_QTL_minusLog10 = -log10(MEGENA[MEGENA$module_QTL_assoc > 0,grep("FET.P_.*BRCA",colnames(MEGENA))])
MEGENA_QTL_minusLog10$name = MEGENA$name[MEGENA$module_QTL_assoc > 0]
MEGENA_QTL_minusLog10_m = melt(MEGENA_QTL_minusLog10,by="name")
#MEGENA_QTL_minusLog10_m$mutation = paste(sep=":",gsub(".*_(.*)_[A-Z][A-Z]","\\1",MEGENA_QTL_minusLog10_m$variable),gsub("FET.P_(.*)_mut.*","\\1",MEGENA_QTL_minusLog10_m$variable))
MEGENA_QTL_minusLog10_m$mutation = gsub("FET.P_(.*)_mut.*","\\1",MEGENA_QTL_minusLog10_m$variable)
MEGENA_QTL_minusLog10_m$data = gsub(".*_BRCA_(.*)_.*_[A-Z][A-Z]","\\1",MEGENA_QTL_minusLog10_m$variable)
MEGENA_QTL_minusLog10_m$Direction = gsub(".*_BRCA_.*_.*_([A-Z][A-Z])","\\1",MEGENA_QTL_minusLog10_m$variable)
table(MEGENA_QTL_minusLog10_m$data)
table(MEGENA_QTL_minusLog10_m$mutation)

MEGENA_QTL_minusLog10_m$data = factor(MEGENA_QTL_minusLog10_m$data,levels=c("Phosphoproteome","Proteome","mRNA"))

fdr.colors=c("NA", "#000000")
getPalette = colorRampPalette(c("#FFFFFF","#ef3b2c","#cb181d","#a50f15","#67000d"))
p = ggplot(data=MEGENA_QTL_minusLog10_m)
p = p + facet_grid(mutation~Direction, drop=T, scales = "free", space = "free")
p = p + geom_tile(aes(y=data, x=name, fill=value), linetype="blank") + scale_fill_gradientn(name= "-log10(P)", na.value=NA, colours=getPalette(100))
p = p + geom_tile(aes(y=data, x=name, color=value>1.30103), fill=NA, size=0.5) + scale_colour_manual(name=paste("P < 0.05"),values = fdr.colors)
p = p + labs(y="Data level", x="Network Module") + theme_nogrid() #+ coord_flip()
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=8))
p
fn = paste('out/MEGENA_network_QTL_BRCA.pdf',sep ="_")
ggsave(file=fn,useDingbats=FALSE)
