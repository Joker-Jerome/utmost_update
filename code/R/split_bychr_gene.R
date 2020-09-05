#!/usr/bin/R
library(R.utils)
library(data.table)
args <- commandArgs()
input_dir <- args[6]
netid <- args[7]
#gunzip(input_dir)
#input_dir_0 <- gsub("\\.gz","",input_dir)
input <- as.data.frame(fread(input_dir, stringsAsFactors=F))
#gzip(input_dir_0)

#input <- read.table(gzfile(input_dir,"r"),header=F,stringsAsFactors=F)
tissue <- gsub("\\.v8\\.normalized_expression\\.bed","",basename(input_dir))
names_id <- colnames(input)[-c(1:4)]
names_id <- gsub("\\.","-",names_id)
for(i in 1:22){
  chr <- paste0("chr",i)
  tmp <- input[which(input[,1]==chr),]
  genes <- tmp[,4]
  for(gene in genes){
    result_dir <- paste0("/ysm-gpfs/pi/zhao-data/",netid,"/GTEx_V8/expr_normalized/", chr, "/", gene,"/")
    dir.create(result_dir,recursive=T)
    expr_tmp <- tmp[which(tmp[,4]==gene),-c(1:4)]
    tmp_df <- data.frame(id=names_id,expr=as.numeric(expr_tmp))
    print(paste0("INFO:", result_dir,tissue,".txt"))
    write.table(tmp_df,file=paste0(result_dir,tissue,".txt"),col.names=F,row.names=F,sep="\t",quote=F)
  }
}
