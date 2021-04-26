options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
chr = as.numeric(args[1])
info_list <- list()
genotype_dir <- "/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/"


count = 0

cat(paste0("INFO: chr ", chr))
info_list[[chr]] <- list()
gtex_dir <- "/gpfs/loomis/project/zhao/zy92/GTEX/" 
chr_str <- paste0("chr", chr, "/")
#gene_vec <- list.files(paste0(genotype_dir, chr_str))
gene_vec = dir(paste0(gtex_dir, "adjusted_expr/chr", chr))

for (gene_id in gene_vec) {
    count = count + 1
    if (count %% 50 == 0) {cat(paste0("INFO: chr", chr, " gene ", count))}
    dose_path <- paste0(genotype_dir, chr_str, gene_id, "/", gene_id)
    Yt <- dir(paste0(gtex_dir, "adjusted_expr/", chr_str, gene_id, "/"))
    P <- length(Yt)
    Y = list()

    for(t in 1:P){
        Y[[t]] = read.table(paste0(gtex_dir, "adjusted_expr/chr", chr, "/", gene_id, "/", Yt[t]), header=F)
    }

    ssize = unlist(lapply(Y, nrow))
    T_num = length(Yt)
    info_list[[chr]][[gene_id]] <- list(T_num, ssize, Y)
    }

save(info_list, 
     file = paste0("/gpfs/loomis/project/zhao/zy92/GTEX/gene_info/info_list_chr", chr, ".RData"))


