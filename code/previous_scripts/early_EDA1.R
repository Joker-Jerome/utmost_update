library(ggplot2)
library(gridExtra)
library(scales)
library(gplots)
library(grid)
options(stringsAsFactors = F)
#save(g_sig2, file="/Users/huyiming/Documents/workspace/MultiGExpr/Codes/glasso2_sig/sig_gene_list.RData")
load("/Users/huyiming/Documents/workspace/MultiGExpr/Codes/glasso2_sig/sig_gene_list.RData")
#setwd("/Users/huyiming/Box Sync/joint test/v3")
tissue <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
ssize = c(298,185,126,197,118,285,72,100, 89,103, 96, 92, 81, 81, 93, 82,183,114,272,124,169,127,241,218,159,190, 97,278,361,256, 85,149, 87, 87,196,302, 77, 89,170,157,278, 70, 79,338)
N <- length(tissue)
trait = c("AD","SCZ","BC")
res_path = "/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/metaxcan_output/"
t = 2
i=0
pt = c('AD Joint', 'SCZ Joint')
plots=list()
for(t in 1:2)
{
  results_meta <- read.table(paste("/Users/huyiming/Box Sync/joint test/metaxcan/", trait[t], ".txt", sep=""))
  results_v3 <- read.table(paste("/Users/huyiming/Box Sync/joint test/v3/", trait[t], ".txt", sep=""))
  results_meta = results_meta[!is.na(results_meta[,i+4]),]
  results_v3 = results_v3[!is.na(results_v3[,i+4]),]
  
  obs_meta <- -log10(sort(results_meta[,i+4],decreasing=F))
  exp_meta <- -log10(1:length(obs_meta)/length(obs_meta))
  data_meta <- cbind(obs_meta, exp_meta)
  sig_meta = sum(obs_meta>-log10(0.05/length(obs_meta)))
  writeLines(results_meta[results_meta[,i+4]<0.05/length(obs_meta),2], paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/",trait[t],'_joint_meta.txt'))
  obs_v3 <- -log10(sort(results_v3[,i+4],decreasing=F))
  exp_v3 <- -log10(1:length(obs_v3)/length(obs_v3))
  data_v3 <- cbind(obs_v3, exp_v3)
  sig_v3 = sum(obs_v3>-log10(0.05/length(obs_v3)))
  writeLines(results_v3[results_v3[,i+4]<0.05/length(obs_v3),2], paste0("/Users/huyiming/Documents/workspace/MultiGExpr/Downstream/",trait[t],'_joint_v3.txt'))
  
  data <- rbind(data_meta, data_v3)
  indi <- c(rep("meta", length(obs_meta)), rep("v3", length(obs_v3)))
  data <- data.frame(data, indi)
  colnames(data) <- c("observed", "expected", "indi")
  
  nam <- pt[t]
  nam <- ggplot(data=data, aes(expected, observed)) + geom_point(aes(colour = factor(indi))) + geom_abline(intercept = 0, slope = 1) + ggtitle(nam) +
    scale_colour_discrete(name  ="Methods",breaks=c("meta", "v3"), labels=c(paste0("meta:#sig=", sig_meta), paste0("v3:#sig=", sig_v3)))
  plots[[t]] = nam
}
png('/Users/huyiming/Box Sync/joint test/v3/AD_SCZ_joint.png', height = 300, width = 600)
grid.arrange(plots[[1]],plots[[2]], nrow = 1, ncol=2)
dev.off()      

trait = c(0.1,0.5,0.8)
for(t in 1:length(trait))
{
  plots=list()
  for(i in 1:45){
    results_meta <- read.table(paste("/Users/huyiming/Box Sync/joint test/simulation/simu_meta_", trait[t], ".txt", sep=""))
    results_v3 <- read.table(paste("/Users/huyiming/Box Sync/joint test/simulation/simu_v3_", trait[t], ".txt", sep=""))
    obs_meta <- -log10(sort(results_meta[,i+3],decreasing=F))
    exp_meta <- -log10(1:length(obs_meta)/length(obs_meta))
    data_meta <- cbind(obs_meta, exp_meta)
    sig_meta = sum(obs_meta>-log10(2.4e-7))
    
    obs_v3 <- -log10(sort(results_v3[,i+3],decreasing=F))
    exp_v3 <- -log10(1:length(obs_v3)/length(obs_v3))
    data_v3 <- cbind(obs_v3, exp_v3)
    sig_v3 = sum(obs_v3>-log10(1.28e-7))
    
    data <- rbind(data_meta, data_v3)
    indi <- c(rep("meta", length(obs_meta)), rep("v3", length(obs_v3)))
    data <- data.frame(data, indi)
    colnames(data) <- c("observed", "expected", "indi")
    
    if(i == 1){
      nam <- "Joint"
    }else{
      nam <- tissue[i-1]
    }
    nam <- ggplot(data=data, aes(expected, observed)) + geom_point(aes(colour = factor(indi))) + geom_abline(intercept = 0, slope = 1) + ggtitle(nam) +
      scale_colour_discrete(name  ="Methods",breaks=c("meta", "v3"), labels=c(paste0("meta:#sig=", sig_meta), paste0("v3:#sig=", sig_v3)))
    plots[[i]] = nam
  }
  png(paste0("/Users/huyiming/Box Sync/joint test/simulation/simu_", trait[t], '_1.png'), height = 900, width = 1200)
  grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],plots[[9]], nrow = 3, ncol=3)
  dev.off()      
  
  png(paste0("/Users/huyiming/Box Sync/joint test/simulation/simu_", trait[t], '_2.png'), height = 900, width = 1200)
  grid.arrange(plots[[10]],plots[[11]],plots[[12]],plots[[13]],plots[[14]],plots[[15]],plots[[16]],plots[[17]],plots[[18]], nrow = 3, ncol=3)
  dev.off()
  
  png(paste0("/Users/huyiming/Box Sync/joint test/simulation/simu_", trait[t], '_3.png'), height = 900, width = 1200)
  grid.arrange(plots[[19]],plots[[20]],plots[[21]],plots[[22]],plots[[23]],plots[[24]],plots[[25]],plots[[26]],plots[[27]], nrow = 3, ncol=3)
  dev.off()
  
  png(paste0("/Users/huyiming/Box Sync/joint test/simulation/simu_", trait[t], '_4.png'), height = 900, width = 1200)
  grid.arrange(plots[[28]],plots[[29]],plots[[30]],plots[[31]],plots[[32]],plots[[33]],plots[[34]],plots[[35]],plots[[36]], nrow = 3, ncol=3)
  dev.off()
  
  png(paste0("/Users/huyiming/Box Sync/joint test/simulation/simu_", trait[t], '_5.png'), height = 900, width = 1200)
  grid.arrange(plots[[37]],plots[[38]],plots[[39]],plots[[40]],plots[[41]],plots[[42]],plots[[43]],plots[[44]],plots[[45]], nrow = 3, ncol=3)
  dev.off()
  
}
