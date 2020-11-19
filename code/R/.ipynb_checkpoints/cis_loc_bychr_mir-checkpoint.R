#!/usr/bin/R

args <- commandArgs(trailingOnly=TRUE) 
chr <- args[1]
#genes<- read.table("gene_names_noVersion_locs_v8.txt",header=T,stringsAsFactors=F,sep="\t")
#genes<- read.table("/ysm-gpfs/project/wl382/GTEx_v8/genotype/gene_names_noVersion_locs_cisSNPs_v8.txt",header=T,stringsAsFactors=F,sep="\t")
genes <- read.table("gene_names_noVersion_locs_v8.txt",header=T,stringsAsFactors=F,sep="\t")

genes_1 <- genes[genes[,6]==chr,-2]
genes_1 <- unique(genes_1)

#genes_gtex <- read.table("gene_names_v8.txt",header=F,stringsAsFactors=F,sep="\t")
i <- chr

print("INFO: loading .bim file")
paste0("/ysm-gpfs/scratch60/wl382/GTEx/v8/genotype/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_maf0.01_grch38p7_chr",i,".bim") -> bim_file
#table <- read.table(file=paste0("locb38Tolochg19_chr",chr,".txt"),header=T,stringsAsFactors=F,sep="\t")
bim <- read.table(bim_file,header=F,stringsAsFactors=F,sep="\t")
bim <- na.omit(bim)
gene_names <- list.files(paste0("/gpfs/loomis/scratch60/zhao/zy92/GTEX/expr_gtex_mirna_updated/chr",chr,"/"))
remain <- c()

print("INFO: checking the cis snp locations")

for (i in 1:length(gene_names)){
#for(i in 1:dim(genes_1)[1]){
  #gene_in_gtex <- genes_gtex$V1[grep(tmp[1,1],genes_gtex$V1)]

  gene_name_gtex <- gene_names[i]
  gene_name <- gsub("\\.[0-9]+","",gene_names[i])  
  #gene_name <- gsub("\\.[0-9]+","",gene_names[i])
  tmp <- genes_1[genes_1[,1]==gene_name,]
  #gene_in_gtex <- genes_gtex$V1[grep(gene_name,gsub("\\.[0-9]+","",genes_gtex$V1))]
#  tmp_file <- paste0("/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/chr",chr,"/",gene_in_gtex,"/",gene_in_gtex,".txt")
#  tmp_df <- try(read.table(tmp_file,header=T,stringsAsFactors=F))
#  if(!inherits(tmp_df,"try-error")){
#    next()
#  }
  
  if((length(tmp[,1]))==0){
    remain <- c(remain,gene_name)
    next()
  }
  pos_1 <- tmp[,3] - 1000000
  pos_2 <- tmp[,4] + 1000000
#  cat(pos_1,"\t",pos_2,"\n")
  ind_1 <- which(bim[,4]>=pos_1)
  ind_2 <- which(bim[,4]<=pos_2)
  inds <- intersect(ind_1,ind_2)
#   cat(length(inds),"\n")
  if(length(inds)==0){
    next()
  }
  loc <- bim[inds,4]
   
#  gene_in_gtex <- genes_gtex$V1[grep(tmp[1,1],genes_gtex$V1)] 
  tmp_df <- data.frame(chr=rep(chr,length(loc)),loc=loc)
  #loc_new <- table[match(tmp_df[,2],table$loc_hg19),]$loc_b38
  #tmp_df[!is.na(loc_new),2]<- loc_new[!is.na(loc_new)]
  #tmp_df<- na.omit(tmp_df) 
  dir.create(paste0("/gpfs/loomis/project/zhao/zy92/GTEX/genotype/cis_loc_mirna/chr",chr,"/", gene_name_gtex,"/"), recursive = T, showWarnings = F)
  tmp_file <- paste0("/gpfs/loomis/project/zhao/zy92/GTEX/genotype/cis_loc_mirna/chr",chr,"/", gene_name_gtex,"/", gene_name_gtex, ".txt")
  write.table(tmp_df,file=tmp_file,col.names=F,row.names=F,sep="\t",quote=F)
 
}

write.table(remain,file=paste0("/gpfs/loomis/project/zhao/zy92/GTEX/genotype/cis_loc_mirna/chr",chr,"_remainGenes.txt"),col.names=F,row.names=F,sep="\t",quote=F)
print("INFO: complete!")

