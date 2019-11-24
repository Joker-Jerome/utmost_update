options(stringsAsFactors = F)
library(ggplot2)
library(gridExtra)
library(scales)
library(gplots)
library(grid)
library(viridis)  
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/Downstream/fqtl_twas_meta_v3_sig_count3.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/Downstream/num_sig_genes.RData")
relevant_tissue["fqtl"] = sigV1_fqtl
irrelevant_tissue["fqtl"] = sigV2_fqtl
wilcox.test(relevant_tissue$MetaXcan, relevant_tissue$MTIE,alternative = c("less"))$p.value #0.07626001
wilcox.test(irrelevant_tissue$MetaXcan, irrelevant_tissue$MTIE)$p.value #0.4193724

wilcox.test(irrelevant_tissue$TWAS, relevant_tissue$TWAS, paired=T, alternative = c("less"))$p.value #0.2025069
wilcox.test(irrelevant_tissue$MetaXcan, relevant_tissue$MetaXcan, paired=T, alternative = c("less"))$p.value #0.1187447
wilcox.test(irrelevant_tissue$fqtl, relevant_tissue$fqtl, paired=T, alternative = c("less"))$p.value #0.6264268
wilcox.test(irrelevant_tissue$MTIE, relevant_tissue$MTIE, paired=T, alternative = c("less"))$p.value #0.01950006





res_df1_1 = data.frame("Sig.Enriched","MetaXcan",relevant_tissue$MetaXcan)
res_df1_2 = data.frame("Sig.Enriched","MTIE",relevant_tissue$MTIE)
res_df1_3 = data.frame("Sig.Enriched","TWAS",relevant_tissue$TWAS)
res_df1_4 = data.frame("Sig.Enriched","fqtl",relevant_tissue$fqtl)
res_df2_1 = data.frame("NonSig.Enriched","MetaXcan",irrelevant_tissue$MetaXcan)
res_df2_2 = data.frame("NonSig.Enriched","MTIE",irrelevant_tissue$MTIE)
res_df2_3 = data.frame("NonSig.Enriched","TWAS",irrelevant_tissue$TWAS)
res_df2_4 = data.frame("NonSig.Enriched","fqtl",irrelevant_tissue$fqtl)
names(res_df1_1) = c("Group","Method","sig_cnt")
names(res_df1_2) = c("Group","Method","sig_cnt")
names(res_df1_3) = c("Group","Method","sig_cnt")
names(res_df1_4) = c("Group","Method","sig_cnt")
names(res_df2_1) = c("Group","Method","sig_cnt")
names(res_df2_2) = c("Group","Method","sig_cnt")
names(res_df2_3) = c("Group","Method","sig_cnt")
names(res_df2_4) = c("Group","Method","sig_cnt")
res_df = rbind(res_df1_3,res_df1_1,res_df1_4,res_df1_2,res_df2_3,res_df2_1,res_df2_4,res_df2_2)
res_df[,1] = as.factor(res_df[,1])
res_df[,1] = factor(res_df[,1], labels=c("Least Significant","Most Significant"))
res_df[,2] = as.factor(res_df[,2])
res_df[,2] = ordered(res_df[,2], levels = c("TWAS", "MetaXcan", "fqtl", "MTIE"))
save.image("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/Downstream/tmp.RData")
gg = ggplot(res_df, aes(x=Group, y=sig_cnt, fill=Method))
gg = gg + geom_boxplot()
gg = gg + theme_minimal(base_family="Helvetica")
gg = gg + labs(x=NULL, y="#SigGenes")
#gg = gg + scale_fill_viridis(option="magma", discrete=TRUE, begin=0.2, end=0.9)
#gg = gg + annotate("text", x=1, y=100, label="P = 0.315", size=7)
#gg = gg + annotate("text", x=2, y=100, label="p = 1.2e-4", size=7)
gg = gg + theme(axis.text=element_text(size=15), axis.title=element_text(size=20))
gg = gg + theme(legend.text=element_text(size=17), legend.title=element_blank())
gg = gg + ylim(0,300)
gg
pdf('/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/Downstream/twas_meta_fqtl_v3_sig_count.pdf', width = 8, height = 8)
gg
dev.off()
