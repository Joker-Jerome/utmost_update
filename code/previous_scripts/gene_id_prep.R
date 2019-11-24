library(biomaRt)

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
#head(listFilters(grch37))
#listAttributes(grch37)

#gene_info <- getBM(filters= c("ensembl_gene_id","biotype"), attributes= c("ensembl_gene_id",'hgnc_symbol','chromosome_name'), values= list(c("ENSG00000130203.5"),"protein_coding"), mart= grch37)

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


