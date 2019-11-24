#get all protein coding genes
library(biomaRt)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
## get ensembl ID
helper <- function(str){
	unlist(strsplit(str, '\\.'))[1]
}
get_gene_index <- function(chr){
	glist = dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/adjusted_expr1/chr", chr));
	tmp = sapply(glist, helper)
	unname(tmp)
}

pc_gene = list()
pc_num = c()
org_num = c()
for(chr in 1:22){
	esd = get_gene_index(chr)
	gene_info = getBM(filters= c("ensembl_gene_id","biotype"), attributes= c("ensembl_gene_id",'hgnc_symbol','chromosome_name'), values= list(esd,"protein_coding"), mart= grch37)
	pc_gene[[chr]] = gene_info
	pc_num = c(pc_num, nrow(gene_info))
	org_num = c(org_num, length(esd))
}
save(pc_gene, pc_num, org_num, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/protein_coding.RData")


load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/protein_coding.RData")
helper <- function(str){
	unlist(strsplit(str, '\\.'))[1]
}
get_esb_index <- function(chr, esbl){
	glist = dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/adjusted_expr1/chr", chr));
	tmp = sapply(glist, helper)
	which(tmp==esbl)
}

all_cmds = c()
for(chr in 1:22){
	for(j in 1:nrow(pc_gene[[chr]])){
		if(pc_gene[[chr]][j,2]!=""){
			gid = pc_gene[[chr]][j,2]
			esbl = pc_gene[[chr]][j,1]
			idx = get_esb_index(chr, esbl)
			all_cmds = c(all_cmds, paste0("module restore Python_R_SQ; Rscript --vanilla /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_prototype.R ", chr, ' ', idx, ' ', gid, ';'))
		}
	}
}


N = length(all_cmds)
N #17886
dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0/"))
for(i in 1:N){
	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0/",i))
	writeLines(all_cmds[i], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0/",i,"/cmd_",i,".sh"))
}
save.image("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/cmds.RData")
17886
general scavenge bigmem pi_zhao

#mkdir all_cmds
#for i in {1..17886}
#do
#mkdir all_cmds/i
#mv cmd_$i.sh all_cmds/i/
#done



cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0/

for i in {2..5000}
do
cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0/$i
sqCreateScript -q general -N name -n 1 -m 20000 -w 24:00:00 cmd_$i.sh > cmd_$i.pbs
sbatch cmd_$i.pbs
done

for i in {5001..10000}
do
cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0/$i
sqCreateScript -q general -N name -n 1 -m 20000 -w 24:00:00 cmd_$i.sh > cmd_$i.pbs
sbatch cmd_$i.pbs
done

for i in {10001..17886}
do
cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0/$i
sqCreateScript -q general -N task_$i -n 1 -m 20000 -w 24:00:00 cmd_$i.sh > cmd_$i.pbs
sbatch cmd_$i.pbs
done



cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds/
sqCreateScript -q general -N name -n 1 -m 20000 -w 24:00:00 cmd_1.sh > pbs_files/cmd_1.pbs
sbatch pbs_files/cmd_1.pbs


### unfinished ones ###
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/cmds.RData")
undone = c()
for(i in 1:length(all_cmds)){
	tmp = strsplit(all_cmds[i], ' ')
	chr = tmp[[1]][7]
	gene = strsplit(tmp[[1]][9], ';')[[1]]
	out = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso0/chr", chr, '/', gene)
	if(length(list.files(out))<2){
		undone = c(undone, i)
		print(i)
	}
}

dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0.2/"))
for(k in 1:length(undone)){
	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0.2/",k))
	writeLines(all_cmds[undone[k]], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0.2/",k,"/cmd_",k,".sh"))	
}

for i in {1..2036}
do
cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0.1/$i
sqCreateScript -q general -N task_$i -n 1 -m 20000 -w 24:00:00 cmd_$i.sh > cmd_$i.pbs
sbatch cmd_$i.pbs
done

for i in {1..267}
do
cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_proto_cmds0.2/$i
sqCreateScript -q general -N task_$i -n 1 -m 20000 -w 24:00:00 cmd_$i.sh > cmd_$i.pbs
sbatch cmd_$i.pbs
done


remaining1 = c()
for(i in 1:6000){
	if(file.exists(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1SQ/",i,"/task",i,".sh.REMAINING"))){
		lr = length(readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1SQ/",i,"/task",i,".sh.REMAINING")))
		if(lr){
			remaining1 = c(remaining1, i)
		}
	}else{
		remaining1 = c(remaining1, i)
	}
}
remaining2 = c()
for(i in 6001:12000){
	if(file.exists(paste0("/ysm-gpfs/pi/zhao/from_louise/ql68/ExpressionPrediction/glasso1SQ/",i,"/task",i,".sh.REMAINING"))){
		lr = length(readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/ql68/ExpressionPrediction/glasso1SQ/",i,"/task",i,".sh.REMAINING")))
		if(lr){
			remaining2 = c(remaining2, i)
		}
	}else{
		remaining2 = c(remaining2, i)
	}
}


remaining3 = c()
for(i in 12001:17886){
	if(file.exists(paste0("/ysm-gpfs/pi/zhao/from_louise/ml237/yiming/SQfile/",i,"/task",i,".sh.REMAINING"))){
		lr = length(readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/ml237/yiming/SQfile/",i,"/task",i,".sh.REMAINING")))
		if(lr){
			remaining3 = c(remaining3, i)
		}
	}else{
		remaining3 = c(remaining3, i)
	}
}

for(i in 1:length(cmds_prefix)){
	if(length(grep("HTRA3", cmds_prefix[i]))){
		print(i)
	}
}
remaining1=as.numeric(readLines("remaining1.txt"))
remaining2=as.numeric(readLines("remaining2.txt"))
remaining = c(remaining1, remaining2)
remaining = remaining3
output_path="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1" ## change to your own dir
sqfile_path="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1SQ_REM" ## change to your own dir
load("/ysm-gpfs/pi/zhao /from_louise/yh367/GTEX/codes/cmds_prefix.RData")
output_cmds = paste(cmds_prefix[remaining], output_path)
dir.create(output_path, showWarnings = FALSE)
dir.create(sqfile_path, showWarnings = FALSE)
setwd(sqfile_path)
for(i in 1:length(output_cmds)){
	dir.create(paste0(i+1373), showWarnings = FALSE)
	writeLines(output_cmds[i], paste0(i+1373, "/task", i+1373, ".sh"))
}
cat("Simple Queue files written to ", sqfile_path, '\n')

1373
sqfile_path="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1SQ_REM"
for i in {1373..2299}
do
cd $sqfile_path/$i
sqCreateScript -q general -N task$i -n 1 -m 25000 -w 24:00:00 task$i.sh > task$i.pbs
sbatch task$i.pbs
done












Rscript --vanilla /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso_prototype1.R 19 554 APOE 50 /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1

