for file in *.tgz; do tar -xvzf $file; done

vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --extract-FORMAT-info GT
vcftools --vcf chr1.info4.vcf --indv-freq-burden --out v4_indv
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --indv-freq-burden --out v2_indv

vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 1 --out v2_chr1
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 2 --out v2_chr2
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 3 --out v2_chr3
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 4 --out v2_chr4
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 5 --out v2_chr5
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 6 --out v2_chr6
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 7 --out v2_chr7
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 8 --out v2_chr8
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 9 --out v2_chr9
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 10 --out v2_chr10
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 11 --out v2_chr11
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 12 --out v2_chr12
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 13 --out v2_chr13
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 14 --out v2_chr14
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 15 --out v2_chr15
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 16 --out v2_chr16
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 17 --out v2_chr17
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 18 --out v2_chr18
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 19 --out v2_chr19
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 20 --out v2_chr20
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 21 --out v2_chr21
vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr 22 --out v2_chr22

for CHR in {1..21}; do vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr $CHR --out "v2_chr"$CHR; done
for CHR in {1..22}; do vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --chr $CHR --recode -c | gzip -c > "v2_chr"$CHR".vcf.gz"; done
cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations
cut -d$'\t' -f1 Adipose_Subcutaneous_Analysis.v6p.all_snpgene_pairs.txt > Adipose_Subcutaneous.gene_id

head -n2000 Adipose_Subcutaneous_Analysis.v6p.all_snpgene_pairs.txt

for file in *.txt; do cut -d$'\t' -f1 $file > "gene_id/"$file".gene"; done
for file in *.txt; do cut -d$'\t' -f2 $file > "variant_id/"$file".variant"; done
for file in *.txt; do cut -d$'\t' -f3 $file > "tss_distance/"$file".tss"; done


## File
filel = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/")
setwd("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/") 
nl = rep(0,44)
for(i in 1:44){
	com <- paste("wc -l ", filel[i], " | awk '{ print $1 }'", sep="")
	nl[i] <- system(command=com, intern=TRUE)
}


options(stringsAsFactors=F)
v2 = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/v2_indv", header=T)
v4 = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/v4_indv", header=T)

#gid1 = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/Testis_Analysis.v6p.all_snpgene_pairs.txt.gene")
#gid2 = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/Brain_Nucleus_accumbens_basal_ganglia_Analysis.v6p.all_snpgene_pairs.txt.gene")
#ugid1 = unique(gid1)[-1]
#ugid2 = unique(gid2)[-1]
#
#Ng = length(ugid2)
#Ns = length(gid2)
#coord = matrix(0, nrow = Ng, ncol = 2)
#ks = 2
#for(i in 1:Ng){
#	coord[i,1] = ks
#	while(gid2[ks]==ugid2[i]){
#		ks = ks + 1
#	}
#	coord[i,2] = ks - 1
#	print(i)
#}
#coord[i,2] = Ns

########### test Testis ##########
gid = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/Testis_Analysis.v6p.all_snpgene_pairs.txt.gene")
ugid = unique(gid)[-1]
Ng = length(ugid)
Ns = length(gid)
coord = matrix(0, nrow = Ng, ncol = 2)
ks = 2
for(i in 1:Ng){
	coord[i,1] = ks
	while(gid[ks]==ugid[i]){
		ks = ks + 1
	}
	coord[i,2] = ks - 1
	print(i)
}
coord[i,2] = Ns

vid = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/variant_id/Testis_Analysis.v6p.all_snpgene_pairs.txt.variant")
for(i in 1:Ng){
	writeLines(vid[coord[i,1]:coord[i,2]], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_snp_list/", ugid[i]))
	print(i)
}



########### all tissues ##########
#T_names = matrix(unlist(strsplit(filel[1:44], "_Analysis")), ncol=2, byrow=T)[,1]
#for(i in 1:44){
#	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_",T_names[i]))
#}
gid_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id")
vid_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/variant_id")
T_names = matrix(unlist(strsplit(gid_fl[1:44], "_Analysis")), ncol=2, byrow=T)[,1]
for(k in 1:11){
	gid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/", gid_fl[k]))
	ugid = unique(gid)[-1]
	Ng = length(ugid)
	Ns = length(gid)
	coord = matrix(0, nrow = Ng, ncol = 2)
	ks = 2
	for(i in 1:Ng){
		coord[i,1] = ks
		while(gid[ks]==ugid[i]){
			ks = ks + 1
			if(ks > Ns){
				break
			}
		}
		coord[i,2] = ks - 1
	}	
	vid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/variant_id/", vid_fl[k]))
	for(i in 1:Ng){
		writeLines(vid[coord[i,1]:coord[i,2]], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_", T_names[k], '/', ugid[i]))
	}
	print(k)
}

for(k in 12:22){
	gid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/", gid_fl[k]))
	ugid = unique(gid)[-1]
	Ng = length(ugid)
	Ns = length(gid)
	coord = matrix(0, nrow = Ng, ncol = 2)
	ks = 2
	for(i in 1:Ng){
		coord[i,1] = ks
		while(gid[ks]==ugid[i]){
			ks = ks + 1
			if(ks > Ns){
				break
			}
		}
		coord[i,2] = ks - 1
	}	
	vid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/variant_id/", vid_fl[k]))
	for(i in 1:Ng){
		writeLines(vid[coord[i,1]:coord[i,2]], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_", T_names[k], '/', ugid[i]))
	}
	print(k)
}

for(k in 23:33){
	gid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/", gid_fl[k]))
	ugid = unique(gid)[-1]
	Ng = length(ugid)
	Ns = length(gid)
	coord = matrix(0, nrow = Ng, ncol = 2)
	ks = 2
	for(i in 1:Ng){
		coord[i,1] = ks
		while(gid[ks]==ugid[i]){
			ks = ks + 1
			if(ks > Ns){
				break
			}
		}
		coord[i,2] = ks - 1
	}	
	vid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/variant_id/", vid_fl[k]))
	for(i in 1:Ng){
		writeLines(vid[coord[i,1]:coord[i,2]], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_", T_names[k], '/', ugid[i]))
	}
	print(k)
}

for(k in (34:44)[-8]){
	gid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/", gid_fl[k]))
	ugid = unique(gid)[-1]
	Ng = length(ugid)
	Ns = length(gid)
	coord = matrix(0, nrow = Ng, ncol = 2)
	ks = 2
	for(i in 1:Ng){
		coord[i,1] = ks
		while(gid[ks]==ugid[i]){
			ks = ks + 1
			if(ks > Ns){
				break
			}
		}
		coord[i,2] = ks - 1
	}	
	vid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/variant_id/", vid_fl[k]))
	for(i in 1:Ng){
		writeLines(vid[coord[i,1]:coord[i,2]], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_", T_names[k], '/', ugid[i]))
	}
	print(k)
}


############ extract tss ############
gid_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id")
tss_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/tss_distance")
T_names = matrix(unlist(strsplit(gid_fl[1:44], "_Analysis")), ncol=2, byrow=T)[,1]
for(i in 1:44){
	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/tss_distance/cis_", T_names[i]))
}
for(k in 1:44){
	gid = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id/", gid_fl[k]))
	ugid = unique(gid)[-1]
	Ng = length(ugid)
	Ns = length(gid)
	coord = matrix(0, nrow = Ng, ncol = 2)
	ks = 2
	for(i in 1:Ng){
		coord[i,1] = ks
		while(gid[ks]==ugid[i]){
			ks = ks + 1
			if(ks > Ns){
				break
			}
		}
		coord[i,2] = ks - 1
	}	
	tss = readLines(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/tss_distance/", tss_fl[k]))
	for(i in 1:Ng){
		writeLines(tss[coord[i,1]:coord[i,2]], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/tss_distance/cis_", T_names[k], '/', ugid[i]))
	}
	print(k)
}




############# test -- Testis ENSG00000235339.1 ############
cis_all = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_Testis/ENSG00000235339.1")
cis_org = matrix(unlist(strsplit(cis_all, "_")), ncol=5, byrow=T)
names(cis_org) = c('chr', 'bp', 'A1', 'A2', 'Assem')
cis_snp = cis_org[(nchar(cis_org[,3])==1)&(nchar(cis_org[,4])==1),]
cis_no_ambigu = cis_snp[apply(cis_snp, 1, function(x){(length(intersect(c(x[3],x[4]), c('A','T')))<2)&(length(intersect(c(x[3],x[4]), c('C','G')))<2)}),]
writeLines(paste(cis_no_ambigu[,1], cis_no_ambigu[,2], sep='\t'), paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_snps/Testis/chr", cis_no_ambigu[1,1], '_', "ENSG00000235339.1"))

#dir.create("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_snps")
#gid_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id")
#vid_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/variant_id")
#T_names = matrix(unlist(strsplit(gid_fl[1:44], "_Analysis")), ncol=2, byrow=T)[,1]
#for(i in 1:44){
#	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_snps/",T_names[i]))
#}

cis_snp = function(tissue, gene, input_dir){
	cis_all = readLines(paste0(input_dir, tissue, '/', gene))
	if(length(cis_all)==1){
		return(list())
	}
	cis_org = matrix(unlist(strsplit(cis_all, "_")), ncol=5, byrow=T)
	names(cis_org) = c('chr', 'bp', 'A1', 'A2', 'Assem')
	cis_snp = cis_org[(nchar(cis_org[,3])==1)&(nchar(cis_org[,4])==1),]
	cis_no_ambigu = cis_snp[apply(cis_snp, 1, function(x){(length(intersect(c(x[3],x[4]), c('A','T')))<2)&(length(intersect(c(x[3],x[4]), c('C','G')))<2)}),]
	if(!is.matrix(cis_no_ambigu)){
		cis_no_ambigu = matrix(cis_no_ambigu, nrow = 1)
	}
#	writeLines(paste(cis_no_ambigu[,1], cis_no_ambigu[,2], sep='\t'), paste0(output_dir, '/chr', cis_no_ambigu[1,1], '_', tissue))
	return(list(site=paste(cis_no_ambigu[,1], cis_no_ambigu[,2], sep='\t'), chr=cis_no_ambigu[1,1]))
}

cis_snp = function(tissue, gene, input_dir){
	cis_all = readLines(paste0(input_dir, tissue, '/', gene))
	cis_org = matrix(unlist(strsplit(cis_all, "_")), ncol=5, byrow=T)
	names(cis_org) = c('chr', 'bp', 'A1', 'A2', 'Assem')
	cis_snp = cis_org[(nchar(cis_org[,3])==1)&(nchar(cis_org[,4])==1),]
	cis_no_ambigu = cis_snp[apply(cis_snp, 1, function(x){(length(intersect(c(x[3],x[4]), c('A','T')))<2)&(length(intersect(c(x[3],x[4]), c('C','G')))<2)}),]
	return(list(site=paste(cis_no_ambigu[,1], cis_no_ambigu[,2], sep='\t'), chr=cis_no_ambigu[1,1]))
}

gid_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id")
T_names = matrix(unlist(strsplit(gid_fl[1:44], "_Analysis")), ncol=2, byrow=T)[,1]
Tg = list()
for(k in 1:44){
	Tg[[k]] = dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_", T_names[k]))
}

g_all = unique(unlist(Tg))
g_freq = table(unlist(Tg))

lower = seq(3, 39768, by = 4419)
upper = c(lower[-1]-1, 39768)

lower = seq(15280, 17678, by = 500)
upper = c(lower[-1]-1, 17678)
t = 3
for(k in lower[t]:upper[t]){
	tryCatch({
		u_tmp = c(); i_tmp = c(); i_ind = 0;
		for(i in 1:44){
			if(g_all[k] %in% Tg[[i]]){
				## save cis of each tissue
				tmp = cis_snp(T_names[i], g_all[k], "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_");
				if(length(tmp) != 0){
					dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr), showWarnings = FALSE);
					dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k]), showWarnings = FALSE);
					writeLines(tmp$site, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k], '/', T_names[i]));
					u_tmp = union(u_tmp, tmp$site);
					if(i_ind==0){
						i_tmp = tmp$site;
						i_ind = 1;
					}else{
						i_tmp = intersect(i_tmp, tmp$site);
					}
				}
			}
		}
		## save union of each tissue
		writeLines(u_tmp, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k], "/union"));
		## save intersection of each tissue
		writeLines(i_tmp, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k], "/intersect"));
		print(k);
	},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

for(k in lower[t]:upper[t]){
	u_tmp = c(); i_tmp = c(); i_ind = 0;
	for(i in 1:44){
		if(g_all[k] %in% Tg[[i]]){
			## save cis of each tissue
			tmp = cis_snp(T_names[i], g_all[k], "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/cis_")
			if(length(tmp) != 0){
				dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr), showWarnings = FALSE);
				dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k]), showWarnings = FALSE);
				writeLines(tmp$site, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k], '/', T_names[i]));
				u_tmp = union(u_tmp, tmp$site);
				if(i_ind==0){
					i_tmp = tmp$site
					i_ind = 1;
				}else{
					i_tmp = intersect(i_tmp, tmp$site)
				}
			}
		}
	}
	## save union of each tissue
	writeLines(u_tmp, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k], "/union"))
	## save intersection of each tissue
	writeLines(i_tmp, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", tmp$chr, '/', g_all[k], "/intersect"))
	print(k)
}


### vcf --positions ###
cmds = c()
chrom = c(1:22)
for(chr in chrom){
	gid = dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", chr))
	chr_tmp = paste0("source ~/.bashrc; cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", chr, "/", gid, "; vcftools --gzvcf /ysm-gpfs/pi/zhao/from_louise/yh367/dbGaP/dbGaP_data/6609/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/v2_chr", chr, ".vcf.gz --positions intersect --recode -c | gzip -c > ", gid, ".vcf.gz")
	cmds = c(cmds, chr_tmp)
}
DosageConvertor --vcfDose ENSG00000000457.9.vcf.gz --prefix ENSG00000000457.9 --type mach --format DS
dosegz = gzfile("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr1/ENSG00000000457.9/ENSG00000000457.9.mach.dose.gz", "rt")
dose = read.table(dosegz)

info = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr1/ENSG00000000457.9/ENSG00000000457.9.mach.info", header=T, sep='\t')
ites = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr1/ENSG00000000457.9/intersect", header=F)
xx = matrix(unlist(strsplit(info[,1], '_')), ncol=5, byrow=T)
####################################### gene expr ##############################################
options(stringsAsFactors=F)
v2 = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/v2_indv", header=T)[,1]
v2xx = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/v2_indv", header=T)
#grpkm = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/dbGaP/dbGaP_data/6609/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/ExpressionFiles/phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct")
grpkm = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct")
grpkm.hd = unlist(strsplit(grpkm[3], '\t'))[-c(1,2)]
#grpkm.eg = unlist(strsplit(grpkm[4], '\t'))
tt = unlist(lapply(strsplit(grpkm.hd, '-'), length))
grpkm.sample.id = unlist(lapply(strsplit(grpkm.hd, '-'), function(x){x[2]}))
v2.subject.id = matrix(unlist(strsplit(v2, '-')), ncol=5, byrow=T)[,2]
#trpkm = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/dbGaP/dbGaP_data/6609/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/ExpressionFiles/phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20150112_RNAseq_Flux1.6_transcript_rpkm.txt")
#gtrpkm = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct")
#gtrpkm.hd = unlist(strsplit(gtrpkm[3], '\t'))[-c(1,2)]

sample.info = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/dbGaP/6609/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/ExpressionFiles/phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_metrics.tsv", header=T, sep='\t')

#cd /ysm-gpfs/pi/zhao/from_louise/yh367/dbGaP/dbGaP_data/6609/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/ExpressionFiles/phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/
#cut -d$'\t' -f1 GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct > GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.id
#
#cut -d$'\t' -f1 /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct > /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.id


sid_match = rep('no', length(grpkm.sample.id))
for(i in 1:length(grpkm.sample.id)){
	if(grpkm.sample.id[i]%in%v2.subject.id){
		sid_match[i] = v2[which(v2.subject.id==grpkm.sample.id[i])]
	}
}
#sum(grpkm.sample.id==sample.info[,1])==length(grpkm.sample.id)

tis_match = rep('no', length(grpkm.hd))
for(i in 1:length(grpkm.hd)){
	if(grpkm.hd[i]%in%sample.info[,1]){
		tis_match[i] = sample.info[,2][which(sample.info[,1]==grpkm.hd[i])]
	}
}

gid_fl = list.files("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id")
T_names = matrix(unlist(strsplit(gid_fl[1:44], "_Analysis")), ncol=2, byrow=T)[,1]
tis_idx = list()
T_list = unique(sample.info[,2])
T_tmp0 = strsplit(T_list, "[^[:alnum:]]")
T_tmp1 = unlist(lapply(T_tmp0, paste, collapse=""))
#T_tmp2 = strsplit(T_tmp1, " ")
T_tmp2 = unlist(lapply(strsplit(T_names, '_'), paste, collapse=""))

TT_match = rep(0, 44)
for(i in 1:44){
	if(T_tmp2[i]%in%T_tmp1){
		TT_match[i] = which(T_tmp1==T_tmp2[i])
	}
}
TT_match[18] = 35
cbind(T_list[TT_match], T_names)
sum(T_tmp1[TT_match]==T_tmp2)
for(i in 1:44){
	tis_idx[[i]] = which(tis_match==T_list[TT_match][i]&sid_match!='no')
}
unlist(lapply(tis_idx, length))

#for(k in 4:length(grpkm)){
#	gexpr = unlist(strsplit(grpkm[k], '\t'))
#	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gene_expr/", gexpr[1]))
#	for(i in 1:44){
#		T_sub_id = sid_match[tis_idx[[i]]]
#		T_expr = gexpr[-c(1,2)][tis_idx[[i]]]
#		writeLines(paste(T_sub_id, T_expr, sep='\t')[T_expr!='0'], paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/gene_expr/", gexpr[1], '/', T_names[i]))
#	}
#}

#expr_id = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/dbGaP/dbGaP_data/6609/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/ExpressionFiles/phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.id")
expr_id = readLines("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.id")
expr_id = expr_id[-c(1,2,3)]

for(chr in 15:22){
	chr_list = dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", chr))
	dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/expr_gtex1/chr", chr))
#	tmp0 = which(expr_id %in% chr_list)
#	tmp1 = expr_id[tmp0]
	tmp0 = intersect(expr_id, chr_list)
	expr_id1 = expr_id[expr_id%in%tmp0]
#	if(sum(chr_list1==expr_id1)==length(tmp0)){
#		print("Good!")
#	}else{
#		O1 = order(expr_id1)
#		O2 = order(chr_list1)
#		O3 = order(O1)
#		chr_list1 = chr_list1[O2][O3]
#		if(sum(chr_list1==expr_id1)==length(tmp0)){
#			print("Shit happens but fixed!")
#		}else{
#			print("Super Shit!")
#			break
#		}
#	}
	expr1 = grpkm[3+which(expr_id%in%tmp0)]
	for(k in 1:length(expr_id1)){
		ge = expr_id1[k]
		tiss1 = setdiff(dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr", chr, '/', ge)), c('union', 'intersect'))
		dir.create(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/expr_gtex1/chr", chr, '/', ge))
		tis_idx1 = which(T_names%in%tiss1)
		gexpr1 = unlist(strsplit(expr1[k], '\t'))[-c(1,2)]
		for(tt in tis_idx1){
			T_sub_id = sid_match[tis_idx[[tt]]]
			T_expr = gexpr1[tis_idx[[tt]]]
			output = paste(T_sub_id, T_expr, sep='\t')[T_expr!='0']
			if(length(output) >= 70){
				writeLines(output, paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/expr_gtex1/chr", chr, '/', ge, '/', T_names[tt]))
			}
		}
	}
	print(chr)
}
