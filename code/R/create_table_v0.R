library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
chr = as.numeric(args[1])

ref_df = as.data.frame(fread(paste0("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/chr", chr, "_snp.txt")))

est_file_list = list.files(paste0("/gpfs/scratch60/zhao/zy92/GTEX/output_0715/chr", chr))
output_dir = paste0("/gpfs/project/zhao/zy92/GTEX/weight_sparse/chr", chr)
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}
   
weight_list = list()
weight_snps = list()
gene_vec = c()

for (i in 1:length(est_file_list)) {
    # read est
    if (i %% 20 == 0) {
        print(paste0("INFO: gene ", i))
    }
    gene = est_file_list[i]
    est_dir = paste0("/gpfs/project/zhao/zy92/GTEX/output/chr", chr, "/", est_file_list[i])
    est_file = list.files(est_dir, pattern = "*.est")
    if (length(est_file) > 0) {   
        est_file = paste0(est_dir, "/", est_file)
        est_df = as.data.frame(fread(est_file))
        joint_df = est_df %>%
            inner_join(ref_df, by = c("id" = "variant_id"))
            
        tissue_vec = colnames(est_df)[7:ncol(est_df)]
        tissue_vec = as.character(sapply(tissue_vec, function(x) unlist(strsplit(x, "\\."))[1]))
                                  
        # write est 

        for (tissue_idx in 1:length(tissue_vec)) {
            tissue = tissue_vec[tissue_idx]
            #IRdisplay::display_html(tissue)
            # weight table
            if (!tissue %in% names(weight_list)) {
                weight_list[[tissue]] = data.frame(rsid = joint_df$rs_id_dbSNP151_GRCh38p7,
                                                 gene = gene,
                                                 weight = joint_df[, tissue_idx + 6],
                                                 ref_allele = joint_df$ref.x,
                                                 eff_allele = joint_df$alt.x
                                                ) %>% 
                filter(weight != 0)
                    
                
            } else {
                tmp_df  = data.frame(rsid = joint_df$rs_id_dbSNP151_GRCh38p7,
                                                 gene = gene,
                                                 weight = joint_df[, tissue_idx + 6],
                                                 ref_allele = joint_df$ref.x,
                                                 eff_allele = joint_df$alt.x
                                     ) %>% 
                filter(weight != 0)
                
                weight_list[[tissue]] = rbind(weight_list[[tissue]], tmp_df)
            }
            
            # extra table
            if (!tissue %in% names(weight_snps)) {
                weight_snps[[tissue]] = data.frame(gene = gene,
                                                 genename = gene,
                                                 pred.perf.R2 = 0.5,
                                                 n.snps.in.model = nrow(joint_df),
                                                 pred.perf.pval = 1e-10,
                                                 pred.perf.qval = 1e-10
                                                )
                
            } else {
                tmp_df  = data.frame(gene = gene,
                                                 genename = gene,
                                                 pred.perf.R2 = 0.5,
                                                 n.snps.in.model = nrow(joint_df),
                                                 pred.perf.pval = 1e-10,
                                                 pred.perf.qval = 1e-10
                                                )
                                     
                weight_snps[[tissue]] = rbind(weight_snps[[tissue]], tmp_df)
            }
            
        }
        
    }
                   
                   
}


for (tissue in names(weight_list)) {
    weight_file = paste0(output_dir, "/", tissue, ".weight.txt")
    extra_file = paste0(output_dir, "/", tissue, ".extra.txt")
    write.table(weight_list[[tissue]], file = weight_file, quote = F, row.names = F)
    write.table(weight_snps[[tissue]], file = extra_file, quote = F, row.names = F)


}
