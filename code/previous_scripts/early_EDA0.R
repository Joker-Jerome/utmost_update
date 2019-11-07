### AD study ###
t = 2
tissue <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
ssize = c(298,185,126,197,118,285, 72,100, 89,103, 96, 92, 81, 81, 93, 82,183,114,272,124,169,127,241,218,159,190, 97,278,361,256, 85,149, 87, 87,196,302, 77, 89,170,157,278, 70, 79,338)
N <- length(tissue)
trait = c("AD","SCZ","BC")
res_path = "/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/metaxcan_output/"
AD = list()
for(i in 1:44){
  tmp = read.csv(paste0(res_path,trait[t],'/', tissue[i], '.csv'), header=T)
  tmp = tmp[!duplicated(tmp),]
  tmp = tmp[tmp[,2]%in%g_sig2[[i]],]
  tmp = tmp[!is.na(tmp$pvalue),]
  AD[[i]] = data.frame(tmp$gene_name,tmp$pvalue)
  remove(tmp)
}

pp = c()
for(i in 1:44){
  pp = c(pp, AD[[i]][,2][AD[[i]][,1]=='APOE'])
}
head(AD[[i]])


for(i in 1:44){
  sig_g = AD[[i]][AD[[i]][,2]<=1.28e-7,1]
  writeLines(sig_g, paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_separate/v3/", tissue[i], '.txt'))
}

results_meta <- read.table(paste("/Users/huyiming/Box Sync/joint test/metaxcan/", trait[t], ".txt", sep=""))
for(i in 1:44){
  meta_gene = results_meta[!is.na(results_meta[,i+4]),c(2,i+4)]
  tmp1 = meta_gene[meta_gene[,2]<=2.4e-7,1]
  writeLines(tmp1, paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_separate/meta/", tissue[i], '.txt'))
}


for(i in 1:44){
  meta_gene = results_meta[!is.na(results_meta[,i+4]),c(2,i+4)]
  tmp1 = meta_gene[meta_gene[,2]<=2.4e-7,1]
  tmp2 = AD[[i]][AD[[i]][,2]<=1.28e-7,1]
  output = c("metaxcan_only","overlap","v3_only")
  if(length(tmp1)|length(tmp2)){
    ovp = intersect(tmp1,tmp2)
    pcg_only = setdiff(tmp1, ovp)
    meta_only = setdiff(tmp2, ovp)
    output[1] = paste(output[1], paste(pcg_only, collapse='\t'), sep='\t')
    output[2] = paste(output[2], paste(ovp, collapse='\t'), sep='\t')
    output[3] = paste(output[3], paste(meta_only, collapse='\t'), sep='\t')
    writeLines(output, paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/SCZ_sig/", tissue[i],".txt"))
  }
}






cgene = rep(-1,44)
for(i in 1:44){
  tmp = read.csv(paste0(res_path,trait[t],'/', tissue[i], '.csv'), header=T)
  if(length(tmp[tmp[,2]=='APOE',]$pvalue)){
    cgene[i] = tmp[tmp[,2]=='APOE',]$pvalue
  }
  
}
cgene
results_meta <- read.table(paste("/Users/huyiming/Box Sync/joint test/metaxcan/", trait[t], ".txt", sep=""))
results_meta[results_meta[,2]=='APOE',][5:48]
for(i in 1:44){
  if(!is.na(cgene[i])){
    if(cgene[i]==-1){
      cgene[i]=NA
    }
  }
}
#cgene[cgene==-1] = rep(NA, sum(cgene==-1))
xx = data.frame(tissue,as.numeric(results_meta[results_meta[,2]=='APOE',][5:48]),cgene)
names(xx) = c('tissue', 'metaxcan', 'v3')
write.table(xx, "/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/metaxcan_output/APOE.txt", quote=F, row.names = F, col.names = T, sep = '\t')
tmp = tmp[!is.na(tmp$pvalue),]
tmp[tmp$pvalue<1.277782e-07,2]
tmp[tmp[,2]=='APOE',]$pvalue

load("/Users/huyiming/Desktop/Msh/results3/va_pgc.RData")
pdf("/Users/huyiming/Desktop/Msh/results3/pgc_va_1e-4.pdf",8,8)
par(mfrow=c(2,2))
hist(bip2_tmp$pval, main = 'PGC_bip, pval', xlab='pval', probability=T)
hist(bip1_tmp$PVALUE, main = 'VA_bip, pval', xlab='pval', probability=T)
hist(scz2_tmp$p, main = 'PGC_scz, pval', xlab='pval', probability=T)
hist(scz1_tmp$PVALUE, main = 'VA_scz, pval', xlab='pval', probability=T)
dev.off()

pdf("/Users/huyiming/Desktop/Msh/results3/pgc_va_or_1e-4.pdf",8,8)
par(mfrow=c(2,2))
hist(bip2_tmp$or, main = 'PGC_bip, or (<1e-4)', xlab='or', probability=T)
hist(exp(bip1_tmp$ALT_EFFSIZE), main = 'VA_bip, or (<1e-4)', xlab='or', probability=T)
hist(scz2_tmp$or, main = 'PGC_scz, or (<1e-4)', xlab='pval', probability=T)
hist(exp(scz1_tmp$ALT_EFFSIZE), main = 'VA_scz, or (<1e-4)', xlab='or', probability=T)
dev.off()

pdf("/Users/huyiming/Desktop/Msh/results3/pgc_va_or_all.pdf",8,8)
par(mfrow=c(2,2))
hist(bip_or2, main = 'PGC_bip, or (all)', xlab='or', probability=T)
hist(bip_or1, main = 'VA_bip, or (all)', xlab='or', probability=T)
hist(scz_or2, main = 'PGC_scz, or (all)', xlab='or', probability=T)
hist(scz_or1, main = 'VA_scz, or (all)', xlab='or', probability=T)
dev.off()
