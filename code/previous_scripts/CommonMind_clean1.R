EXPR="/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/rna_tmp/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/"
"EXCLUDE_ANCESTRY_NoSVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv"

/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/rna_tmp/DorsolateralPrefrontalCortex/CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv

## 1. generate dosage matrix for each gene ##
cd /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed
tar -xzvf CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_Imputation1000genomes_chr21.tar.gz
for i in {1..20}
do
tar -xzvf CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_Imputation1000genomes_chr$i.tar.gz
done

## 1.1 generate merge commands for gtool ##
cmds = c()
for(chr in 1:18){
	chrdir = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr",chr)
	setwd(chrdir)
	lsd = dir(chrdir)
	info_index = sample_index = merged_index = c()
	for(i in 1:length(lsd)){
		if(length(grep(".gen_info",lsd[i]))){
			info_index = c(info_index, i)
		}
		if(length(grep(".sample",lsd[i]))){
			sample_index = c(sample_index, i)
		}
		if(length(grep("merged",lsd[i]))){
			merged_index = c(merged_index, i)
		}
	}
	infofile = lsd[info_index]
	samplefile = rep(lsd[sample_index], length(infofile))
	genfile = lsd[-c(info_index,sample_index,merged_index)]
	gen_cat = paste(genfile, collapse=" ")
	sample_cat = paste(samplefile, collapse=" ")
	#dir.create(paste0(chrdir, "/merged"))
	ini_gen = paste0("source ~/.bashrc; cd ", chrdir, "; gtool -M --g ", gen_cat," --s ", sample_cat, " --og merged/chr", chr, ".gen --os merged/chr", chr, ".sample")
	cmds = c(cmds, ini_gen)
	#ini_info = readLines(infofile[1])
	#for(i in 2:length(infofile)){
	#	tmp = readLines(infofile[i])[-1]
	#	ini_info = c(ini_info, tmp)
	#}
	#writeLines(ini_info, paste0(chrdir, "/merged/chr", chr, ".gen_info"))
}
dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/gtool_sq/"))
for(i in 1:18){
	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/gtool_sq/",i))
	writeLines(cmds[i], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/gtool_sq/",i,'/task',i,'.sh'))
}
for i in {1..18}
do
cd /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/gtool_sq/$i
sqCreateScript -q general -N task$i -n 1 -m 40000 -w 24:00:00 task$i.sh > task$i.pbs
sbatch task$i.pbs
done

## 1.2 convert merged .gen to .vcf.gz with qctool ##
"source ~/.bashrc; cd /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr",1:21,"/merged; qctool -g chr",1:21,".gen -og chr",1:21,".vcf; gzip chr",1:21,".vcf"
for i in {1..20}
do
mkdir /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/qctool_sq/$i
cd /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/qctool_sq/$i
echo "source ~/.bashrc; cd /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr$i/merged; qctool -g chr$i.gen -og chr$i.vcf; gzip chr$i.vcf" > task$i.sh
sqCreateScript -q general -N task$i -n 1 -m 40000 -w 24:00:00  task$i.sh > task$i.pbs
sbatch task$i.pbs
done
#awk -F' ' '{print NF; exit}' chr22.gen 

## 1.3 extract by position ##
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/cmds_prefix1.RData")
cmds = sapply(cmds_prefix, xx)
names(cmds) = NULL
tsk_idx = 1:17886
output_cmds = paste(cmds[tsk_idx], output_path)
dir.create(output_path, showWarnings = FALSE)
dir.create(sqfile_path, showWarnings = FALSE)
setwd(sqfile_path)
for(i in 1:length(tsk_idx)){
	dir.create(paste0(i), showWarnings = FALSE)
	writeLines(output_cmds[i], paste0(i, "/task", i, ".sh"))
}
cat("Simple Queue files written to ", sqfile_path, '\n')
cat("# of tasks: ", length(output_cmds), '\n')

xx = function(x){
	tmp = unlist(strsplit(x,' '))
	tmp[7:9]
}

chrV = gidV = gsymbolV = c()
for(i in 1:length(cmds_prefix)){
	tmp = xx(cmds_prefix[i])
	chrV = c(chrV, as.numeric(tmp[1]))
	gidV = c(gidV, as.numeric(tmp[2]))
	gsymbolV = c(gsymbolV, tmp[3])
}
save(chrV, gidV, gsymbolV, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/extract_vcf.RData")

load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/extract_vcf.RData") ## contains three variables: chrV, gidV, gsymbolV
ngenes = length(chrV)
## CHANGE this one ##
output_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/test/"
dir.create(output_path, showWarnings = FALSE)
cmds = c()
for(i in 1:ngenes){
	chr = chrV[i]
	k = gidV[i]
	gene_id = gsymbolV[i]
	glist = dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", chr));
	gid = glist[k];
	dir.create(paste0(output_path, "chr", chr), showWarnings = FALSE)
	dir.create(paste0(output_path, "chr", chr, "/", gene_id), showWarnings = FALSE)
	#gid = dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", chr))
	chr_tmp = paste0("source ~/.bashrc; vcftools --gzvcf /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr",chr,"/merged/chr", chr, ".vcf.gz --positions /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", chr, "/", gid, "/intersect --recode -c | gzip -c > ", output_path, "chr", chr, "/", gene_id, "/", gene_id, ".vcf.gz")
	cmds = c(cmds, chr_tmp)
}
## will generate ~17,000 commands, you can later parallalize them with simple queue ##

## DosageConvertor commands ##
DosageConvertor --vcfDose SCYL3.vcf.gz --prefix SCYL3 --type mach --format DS

#source ~/.bashrc; vcftools --gzvcf /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr1/merged/chr1.vcf.gz --positions /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr1/ENSG00000000457.9/intersect --recode -c | gzip -c > /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/test/chr1/SCYL3/SCYL3.vcf.gz
#source ~/.bashrc; vcftools --vcf /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr22/merged/chr22.vcf --chr 22 --out v2_chr22
#awk '{if($0 !~ /^#/) print "22"$0; else print $0}' chr22.vcf > with_chr.vcf
#awk '{if($0 !~ /^#/) gsub("NA","22",$0)}1' chr22.vcf > with_chr.vcf

## add chr number to vcf files (should have spent more time reading qctool instruction)
for i in {1..20}
do
mkdir /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr$i/merged/backup
cd /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr$i/merged/
gunzip chr$i.vcf.gz
awk '{if($0 !~ /^#/) gsub("NA", f ,$0)}1' f="$i" chr$i.vcf > with_chr.vcf
mv chr$i.vcf backup/
mv with_chr.vcf chr$i.vcf
gzip chr$i.vcf
done

## gunzip all dose files
dose_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/vcf_out/"
#setwd(dose_path)
for(chrom in 1:22){
	gene_path = paste0(dose_path, "chr", chrom)
	gene_chrom = dir(gene_path)
	for(g in gene_chrom){
		gpath = paste0(gene_path, '/', g, '/')
		system(paste0("cd ", gpath, "; gunzip ", g, ".mach.dose.gz"))
	}
}

## 2. regenerate weight matrix with rs number for both single and cross-tissue model ##
## single-tissue weights path: /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/single_rerun/results
## multi-tissue weights path: # yh367
## rerun single tissue elastic net
module load SimpleQueue/3.1
start=1
end=17886
output_path="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/single_rerun/results" ## change to your own dir
sqfile_path="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/single_rerun/sq" ## change to your own dir
for i in {1..17886}
do
cd $sqfile_path/$i
sqCreateScript -q general -N task$i -n 1 -m 20000 -w 24:00:00 task$i.sh > task$i.pbs
sbatch task$i.pbs
done
## see map_rs_to_weights.py
python map_rs_to_weights.py

## 3. prediction and evaluation ##

## 3.1 merge several id files
## from dose file: /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/vcf_out/chr1/A3GALT2/A3GALT2.mach.dose
## from vcf file: /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr1/merged/chr1.sample
## from expr file: /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv
options(stringsAsFactors=F)
dose = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/vcf_out/chr1/A3GALT2/A3GALT2.mach.dose", header=F)
dose_ID = paste0("sample_",1:621)
vcf = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/chr1/merged/chr1.sample", header=T)
vcf_ID = vcf[-1,1]
ID_map = data.frame(dose_ID, vcf_ID)
expr = read.csv("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv")
expr_ID = expr[,c(1,2)]
clinical = read.csv("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Clinical/CMC_MSSM-Penn-Pitt_Clinical.csv")
ctrl = clinical[clinical$Dx=="Control",1]
expr_ID = expr_ID[expr_ID[,1]%in%ctrl,]
shit <- function(xx){
	tmp0 = unlist(strsplit(xx, '_'))
	center = tmp0[2]
	num = as.numeric(tmp0[3])
	paste(center,num,sep='_')
}
tmp0 = sapply(expr_ID[,1], shit)
names(tmp0) = NULL
expr_ID = data.frame(tmp0, expr_ID[,2])
names(expr_ID) = c("vcf_ID", "DLPFC_RNA_Sequencing_Sample_ID")
ID_map_merged = merge(ID_map, expr_ID, by = "vcf_ID")
#ID_map[!ID_map[,2]%in%ID_map_merged[,1],2]
#expr_ID[!expr_ID[,1]%in%vcf_ID,1]
dose_ID2 = matrix(unlist(strsplit(ID_map_merged[,2],'_')), ncol=2, byrow=T)[,2]
ID_map_merged2 = cbind(ID_map_merged,dose_ID2)
write.table(ID_map_merged2[,c(1,2,4,3)], "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/genotype_expr_ID.map", quote=F, row.names=F, col.names=c("vcf_ID","dose_ID1","dose_ID2","DLPFC_RNA_Sequencing_Sample_ID"))

read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/genotype_expr_ID.map")

## 3.2 get expr files for each gene
## expr files:
options(stringsAsFactors=F)
expr_paths = c("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/EXCLUDE_ANCESTRY_NoSVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv",
	"/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/EXCLUDE_ANCESTRY_SVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv",
	"/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/REQUIRE_ANCESTRY_NoSVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-dataNormalization-includeAncestry-adjustedLogCPM.tsv",
	"/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/REQUIRE_ANCESTRY_SVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-includeAncestry-adjustedLogCPM.tsv")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/protein_coding.RData")
#pc_gene, pc_num, org_num
ENSG_SYBL = data.frame()
for(i in 1:22){
	ENSG_SYBL = rbind(ENSG_SYBL, pc_gene[[i]])
}
output_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/expr_by_gene/"
dir.create(output_path, showWarnings = FALSE)
for(i in 1:4){
	expr = read.table(expr_paths[i], header=T)
	names(expr)[1] = "ensembl_gene_id"
	expr_merged = merge(ENSG_SYBL, expr, by="ensembl_gene_id")
	for(chrom in 1:22){
		dir.create(paste0(output_path,"chr",chrom), showWarnings = FALSE)
		expr_chrom = expr_merged[expr_merged["chromosome_name"]==chrom,]
		for(g in 1:dim(expr_chrom)[1]){
			goutpath = paste0(output_path, "chr", chrom, '/', expr_chrom["hgnc_symbol"][g,], '/')
			dir.create(goutpath, showWarnings = FALSE)
			exprV = expr_chrom[g,-c(1,2,3)]
			sids = names(exprV)
			exprV = as.numeric(exprV)
			write.table(data.frame(sids,exprV), paste0(goutpath, "expr",i,".txt"), quote=F, row.names=F, col.names=c("DLPFC_RNA_Sequencing_Sample_ID", "Expr"))
		}
	}
}


## 3.3 predict gene expr with single/multi weights
## several paths needed: (dorsolateral prefrontal cortex / Brain_Frontal_Cortex_BA9)
## genotypes: /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/vcf_out/
## multi-weights: /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi/
## single-weights: /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/single/

## multi-weight prediction
options(stringsAsFactors=F)
weight_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi/"
gtype_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/vcf_out/"
with_pred1 = list()
for(chrom in 1:22){
	with_pred1[[chrom]] = c("")
}
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/with_pred1.RData")
#chrom=22
for(chrom in 1:22){
	glist = dir(paste0(weight_path, "chrom", chrom))
	for(g in glist){
		weights = read.table(paste0(weight_path, "chrom", chrom, '/', g, '/', g, '.est'), header=T)
		cond = ("Brain_Frontal_Cortex_BA9.adj_expr"%in%names(weights))&(sum(weights$Brain_Frontal_Cortex_BA9.adj_expr!=0))
		if(cond){
			idx = which(names(weights)=="Brain_Frontal_Cortex_BA9.adj_expr")
			genotypes = read.table(paste0(gtype_path, "chr", chrom, '/', g, '/', g, '.mach.dose'), header=F)
			info = read.table(paste0(gtype_path, "chr", chrom, '/', g, '/', g, '.mach.info'), header=T, sep='\t')
			info_trimmed = cbind(info[,1:3], t(genotypes[,-c(1,2)]))
			names(info_trimmed) = c("RS", "REF", "ALT", paste0("sample_", 1:621))
			wB = weights[,c(1,2,3,ncol(weights),idx)]
			coefX = merge(wB, info_trimmed, by="RS")
			qc1 = (coefX$REF.0. == coefX$ALT)&(coefX$ALT.1. == coefX$REF)
			cat(sum(qc1), "SNPs flipped due to reverse alleles\n")
			if(sum(qc1)){
				coefX[qc1, 8:ncol(coefX)] = 2 - coefX[qc1, 8:ncol(coefX)]
				coefX[qc1,]$REF = coefX[qc1,]$ALT.1.
				coefX[qc1,]$ALT = coefX[qc1,]$REF.0.	
			}		
			qc2 = (coefX$REF.0. == coefX$REF)&(coefX$ALT.1. == coefX$ALT)
			cat(nrow(coefX) - sum(qc2), "SNPs discarded due to mismatched alleles\n")
			if(nrow(coefX) - sum(qc2)){
				coefX = coefX[qc2,]
			}
			predicted = as.numeric(t(coefX[,8:ncol(coefX)])%*%coefX$Brain_Frontal_Cortex_BA9.adj_expr)
			write.table(data.frame(names(coefX)[8:ncol(coefX)],predicted), paste0(weight_path, "chrom", chrom, '/', g, '/', g, '.pred'), quote=F, row.names=F, col.names=c("dose_ID1", "predicted"))
			with_pred1[[chrom]] = c(with_pred1[[chrom]], g)
		}
	}
}
save(with_pred1, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/with_pred1.RData")


options(stringsAsFactors=F)
weight_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/single/"
gtype_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/Imputed/vcf_out/"
with_pred2 = list()
for(chrom in 1:22){
	with_pred2[[chrom]] = c("")
}

for(chrom in 1:22){
	glist = dir(paste0(weight_path, "chrom", chrom))
	for(g in glist){
		weights = read.table(paste0(weight_path, "chrom", chrom, '/', g, '/', g, '.est'), header=T)
		cond = ("Brain_Frontal_Cortex_BA9.adj_expr"%in%names(weights))&(sum(weights$Brain_Frontal_Cortex_BA9.adj_expr!=0))
		if(cond){
			idx = which(names(weights)=="Brain_Frontal_Cortex_BA9.adj_expr")
			genotypes = read.table(paste0(gtype_path, "chr", chrom, '/', g, '/', g, '.mach.dose'), header=F)
			info = read.table(paste0(gtype_path, "chr", chrom, '/', g, '/', g, '.mach.info'), header=T, sep='\t')
			info_trimmed = cbind(info[,1:3], t(genotypes[,-c(1,2)]))
			names(info_trimmed) = c("RS", "REF", "ALT", paste0("sample_", 1:621))
			wB = weights[,c(1,2,3,ncol(weights),idx)]
			coefX = merge(wB, info_trimmed, by="RS")
			qc1 = (coefX$REF.0. == coefX$ALT)&(coefX$ALT.1. == coefX$REF)
			cat(sum(qc1), "SNPs flipped due to reverse alleles\n")
			if(sum(qc1)){
				coefX[qc1, 8:ncol(coefX)] = 2 - coefX[qc1, 8:ncol(coefX)]
				coefX[qc1,]$REF = coefX[qc1,]$ALT.1.
				coefX[qc1,]$ALT = coefX[qc1,]$REF.0.	
			}		
			qc2 = (coefX$REF.0. == coefX$REF)&(coefX$ALT.1. == coefX$ALT)
			cat(nrow(coefX) - sum(qc2), "SNPs discarded due to mismatched alleles\n")
			if(nrow(coefX) - sum(qc2)){
				coefX = coefX[qc2,]
			}
			predicted = as.numeric(t(coefX[,8:ncol(coefX)])%*%coefX$Brain_Frontal_Cortex_BA9.adj_expr)
			write.table(data.frame(names(coefX)[8:ncol(coefX)],predicted), paste0(weight_path, "chrom", chrom, '/', g, '/', g, '.pred'), quote=F, row.names=F, col.names=c("dose_ID1", "predicted"))
			with_pred2[[chrom]] = c(with_pred2[[chrom]], g)
		}
	}
}
save(with_pred2, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/with_pred2.RData")

module restore Python_R_SQ; Rscript /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/sq0/single_pred.R
sqCreateScript -q general -N single_pred -n 1 -m 40000 -w 24:00:00 single.sh > single.sh.pbs
sbatch single.sh.pbs


## 3.4 merge predicted with actual gene expr
## multi-weights: /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi/
## single-weights: /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/single/
## actual: /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/expr_by_gene/
options(stringsAsFactors=F)
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/with_pred1.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/with_pred2.RData")
sum(unlist(lapply(with_pred1,length)))
sum(unlist(lapply(with_pred2,length)))
ovp_genes = list()
ovp_genes_df = data.frame()
ovp_genes_num = rep(0,22)
for(chrom in 1:22){
	ovp_genes[[chrom]] = intersect(intersect(with_pred1[[chrom]][-1],with_pred2[[chrom]][-1]),dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/expr_by_gene/chr", chrom)))
	ovp_genes_num[chrom] = length(ovp_genes[[chrom]])
	ovp_genes_df = rbind(ovp_genes_df,cbind(chrom, ovp_genes[[chrom]]))
}
write.table(ovp_genes_df, "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/merged_gene_list.txt", quote=F, row.names=F, col.names=c("chrom","gene"))
save(ovp_genes, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/useful_gene_list.RData")


### calculate PEER and PCs ##
## PEER (15) ##
options(stringsAsFactors=F)
library(peer)
expr_paths = c("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/EXCLUDE_ANCESTRY_NoSVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv",
	"/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/EXCLUDE_ANCESTRY_SVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv",
	"/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/REQUIRE_ANCESTRY_NoSVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-dataNormalization-includeAncestry-adjustedLogCPM.tsv",
	"/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/REQUIRE_ANCESTRY_SVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-includeAncestry-adjustedLogCPM.tsv")
expr = read.csv("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv")
clinical = read.csv("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Clinical/CMC_MSSM-Penn-Pitt_Clinical.csv")
IDs_map = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/genotype_expr_ID.map", header=T)
for(j in 1:4){
	tmp = read.table(expr_paths[j], header=T)
	xxx = names(tmp)[-1]
	tmp = t(tmp[,-1])
	tmp = data.frame(xxx,tmp)
	names(tmp)[1] = "DLPFC_RNA_Sequencing_Sample_ID"
	tmp = merge(tmp, IDs_map, by="DLPFC_RNA_Sequencing_Sample_ID")
	exprMat = as.matrix(tmp[,2:16424])
	exprMat = unname(exprMat)
	model = PEER()
	PEER_setPhenoMean(model,as.matrix(exprMat))
	dim(PEER_getPhenoMean(model))
	PEER_setNk(model,15)
	PEER_getNk(model)
	PEER_update(model)
	factors = PEER_getX(model)
	write.table(data.frame(tmp[,1], factors), paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/",j,".PEER"), quote=F, row.names=F, col.names=c("DLPFC_RNA_Sequencing_Sample_ID", paste0("PEER",1:15)))
}
Sample_ID = tmp[,1]
## PCA (3) ##
IDs_map1 = IDs_map[IDs_map[,4]%in%Sample_ID,]
fam = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/QCd/CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_QCed.fam", header=F)
extr_list = fam[fam[,2]%in%IDs_map1[,1],c(1,2)]
write.table(extr_list, "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/extr_list.txt", quote=F, row.names=F, col.names=F)
cd /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/QCd
plink --bfile CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_QCed --keep /ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/extr_list.txt --make-bed --out gp_filtered
bim = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/QCd/gp_filtered.bim", header=F)
keep_snp = bim[bim[,2]!='.',2]
keep_snp = keep_snp[-1]
writeLines(keep_snp, "/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/QCd/keep_snp.txt")
plink --bfile gp_filtered --extract keep_snp.txt --make-bed --out gp_filtered1
gcta64 --bfile gp_filtered1 --autosome --make-grm --out gp_filtered1_grm
gcta64 --grm gp_filtered1_grm --pca 3 --out gp_filtered1_pca3
## combine the covariate file ##
clinical = read.csv("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Clinical/CMC_MSSM-Penn-Pitt_Clinical.csv")
IDs_map = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/genotype_expr_ID.map", header=T)
tmp1 = merge(IDs_map, clinical, by="DLPFC_RNA_Sequencing_Sample_ID")
tmp1 = tmp1[,c(1:6,10)]
pca_file = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Genotypes/QCd/gp_filtered1_pca3.eigenvec", header=F)
names(pca_file) = c("FID","vcf_ID",paste0("PCA", 1:3))
tmp2 = merge(tmp1, pca_file, by="vcf_ID")
tmp2 = tmp2[,-8]

peer_files = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/",1:4,".PEER")
for(j in 2:4){
	peer_file = read.table(peer_files[j], header=T)
	tmp3 = merge(tmp2, peer_file, by="DLPFC_RNA_Sequencing_Sample_ID")
	tmp3 = tmp3[,-c(2:5)]
	write.table(tmp3, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/",j,".COV"), quote=F, row.names=F, col.names=T)	
}
chrom=5; i=67

results_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi_single_merged/"
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/useful_gene_list.RData")

Rsq = Rsq_raw = list()
for(i in 1:4){
	Rsq[[i]] = list()
	Rsq_raw[[i]] = list()
	for(chrom in 1:22){
		Rsq[[i]][[chrom]] = matrix(0,ovp_genes_num[chrom],2)
		Rsq_raw[[i]][[chrom]] = matrix(0,ovp_genes_num[chrom],2)
	}
}

k = 0
for(chrom in 1:22){
	glist = ovp_genes[[chrom]]
	for(i in 1:ovp_genes_num[chrom]){
		g = glist[i]
		dir.create(paste0(results_path, "chr", chrom, '/', g), showWarnings=F)
		multi_pred = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi/chrom",chrom,'/',g,'/',g,".pred"), header=T)
		single_pred = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/single/chrom",chrom,'/',g,'/',g,".pred"), header=T)
		names(multi_pred) = c("dose_ID1", "multi_pred")
		names(single_pred) = c("dose_ID1", "single_pred")
		act_expr_paths = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/expr_by_gene/chr", chrom, '/', g, "/expr", 1:4, ".txt")
		#tmp0 = merge(IDs_map, single_pred, by="dose_ID1")
		#tmp1 = merge(tmp0, multi_pred, by="dose_ID1")
		for(j in 1:4){
			act_expr = read.table(act_expr_paths[j], header=T)
			#cov_expr = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/",j,".COV"), header=T)
			#expr_with_cov = merge(act_expr, cov_expr, by="DLPFC_RNA_Sequencing_Sample_ID")
			#forreg = expr_with_cov[,-c(1)]
			#expr_with_cov["adjExpr"] = resid(lm(Expr ~ ., data=forreg))
			#tmp2 = merge(tmp1, expr_with_cov, by="DLPFC_RNA_Sequencing_Sample_ID")
			tmp2 = merge(tmp1, expr_with_cov, by="DLPFC_RNA_Sequencing_Sample_ID")
			Rsq[[j]][[chrom]][i,1] = cor(tmp2$single_pred,tmp2$adjExpr)^2
			Rsq[[j]][[chrom]][i,2] = cor(tmp2$multi_pred,tmp2$adjExpr)^2
			Rsq_raw[[j]][[chrom]][i,1] = cor(tmp2$single_pred,tmp2$Expr)^2
			Rsq_raw[[j]][[chrom]][i,2] = cor(tmp2$multi_pred,tmp2$Expr)^2
			write.table(tmp2, paste0(results_path, "chr", chrom, '/', g, "/expr", j, ".txt"))
		}
		k = k + 1
		print(k)
	}
}
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_list.RData")

table(expr_with_cov1$Gender)
table(expr_with_cov1$Ethnicity)
table(expr_with_cov1$Brain_Weight_g)
table(expr_with_cov1$PMI_hrs)
table(expr_with_cov1$pH)
table(expr_with_cov1$Year_of_Autopsy)
table(expr_with_cov1$Age_of_Death)
table(expr_with_cov1$Cause_of_Death)
table(expr_with_cov1$Height_in)
table(expr_with_cov1$Weight_lbs)
table(expr_with_cov1$Genotyping_Sample_ID)
table(expr_with_cov1$Notes)
table(expr_with_cov1$Individual_ID.y)
table(expr_with_cov1$DLPFC_Brodmann_Area)
table(expr_with_cov1$DLPFC_Hemisphere)
table(expr_with_cov1$DLPFC_Dissection_Operator)
table(expr_with_cov1$DLPFC_Date_Dissected)
table(expr_with_cov1$DLPFC_RNA_isolation_Batch)
table(expr_with_cov1$DLPFC_RNA_isolation_PrepDate)
table(expr_with_cov1$DLPFC_RNA_isolation_RIN)
table(expr_with_cov1$DLPFC_RNA_isolation_28S_18S_ratio)
table(expr_with_cov1$DLPFC_RNA_isolation_PrepOperator)
table(expr_with_cov1$DLPFC_RNA_isolation_TotalYield_ug)
table(expr_with_cov1$DLPFC_RNA_isolation_260_280_ratio)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Ribozero_Batch)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Library_Batch)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Flowcell_Batch)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Flowcell_Name)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Flowcell_Lane_A)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Flowcell_Lane_B)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Mapped_Reads)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Exonic_Rate)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Intronic_Rate)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Intergenic_Rate)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Genes_Detected)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Transcripts_Detected)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Expression_Profiling_Efficiency)
table(expr_with_cov1$DLPFC_RNA_Sequencing_rRNA_Rate)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Total_Reads)
table(expr_with_cov1$DLPFC_RNA_Sequencing_Percent_Aligned)


results_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi_single_merged/"
dir.create(results_path, showWarnings=F)

IDs_map = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/genotype_expr_ID.map", header=T)
Rsq = Rsq_raw = list()
for(i in 1:4){
	Rsq[[i]] = list()
	Rsq_raw[[i]] = list()
	for(chrom in 1:22){
		Rsq[[i]][[chrom]] = matrix(0,ovp_genes_num[chrom],2)
		Rsq_raw[[i]][[chrom]] = matrix(0,ovp_genes_num[chrom],2)
	}
}
k = 0
for(chrom in 1:22){
	dir.create(paste0(results_path, "chr", chrom), showWarnings=F)
	glist = ovp_genes[[chrom]]
	for(i in 1:ovp_genes_num[chrom]){
		g = glist[i]
		dir.create(paste0(results_path, "chr", chrom, '/', g), showWarnings=F)
		multi_pred = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi/chrom",chrom,'/',g,'/',g,".pred"), header=T)
		single_pred = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/single/chrom",chrom,'/',g,'/',g,".pred"), header=T)
		names(multi_pred) = c("dose_ID1", "multi_pred")
		names(single_pred) = c("dose_ID1", "single_pred")
		act_expr_paths = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/expr_by_gene/chr", chrom, '/', g, "/expr", 1:4, ".txt")
		tmp0 = merge(IDs_map, single_pred, by="dose_ID1")
		tmp1 = merge(tmp0, multi_pred, by="dose_ID1")
		for(j in 1:4){
			act_expr = read.table(act_expr_paths[j], header=T)
			tmp2 = merge(tmp1, act_expr, by="DLPFC_RNA_Sequencing_Sample_ID")
			Rsq[[j]][[chrom]][i,1] = cor(tmp2$single_pred,tmp2$Expr)^2
			Rsq[[j]][[chrom]][i,2] = cor(tmp2$multi_pred,tmp2$Expr)^2
			write.table(tmp2, paste0(results_path, "chr", chrom, '/', g, "/expr", j, ".txt"))
		}
		k = k + 1
		print(k)
	}
}
save(Rsq, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_list.RData")
Rsq_stack = list()
for(i in 1:4){
	Rsq_stack[[i]] = data.frame()
	for(chrom in 1:22){
		Rsq_stack[[i]] = rbind(Rsq_stack[[i]], Rsq[[i]][[chrom]])
	}
}
save(Rsq_stack, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack.RData")
Rsq_stack_not_na = list()
for(i in 1:4){
	not_na = (!is.na(Rsq_stack[[i]][,1]))&(!is.na(Rsq_stack[[i]][,2]))
#	Rsq_stack_not_na[[i]] = Rsq_stack[[i]][not_na,]
	cat("single.mean = ", median(Rsq_stack[[i]][not_na,1]), "; multi.mean = ", median(Rsq_stack[[i]][not_na,2]), "\n")
}
save(Rsq_stack_not_na, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack_not_na.RData")

#for(chrom in 1:22){
#	glist = ovp_genes[[chrom]]
#	for(i in 1:ovp_genes_num[chrom]){
#		for(j in 1:4){
#			if(is.na(Rsq[[j]][[chrom]][i,1])){
#				cat(chrom,", ", i, "\n")
#			}
#			if(is.na(Rsq[[j]][[chrom]][i,2])){
#				cat(chrom,", ", i, "\n")
#			}
#		}
#	}
#}
final = cbind(ovp_genes_df,Rsq_stack[[1]],Rsq_stack[[2]],Rsq_stack[[3]],Rsq_stack[[4]])
final1 = final[final[,2]%in%common_genes_id[[12]],]
shit = (!is.na(final1[,3]))&(!is.na(final1[,4]))&(!is.na(final1[,5]))&(!is.na(final1[,6]))&(!is.na(final1[,7]))&(!is.na(final1[,8]))&(!is.na(final1[,9]))&(!is.na(final1[,10]))
final2 = final1[shit,]

