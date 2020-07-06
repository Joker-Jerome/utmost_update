library(data.table)
library(dplyr)

# arguments
args = commandArgs(trailingOnly=TRUE)
task_index = as.numeric(args[1]) ## a file contains task index
#task_index = 1 ## a file contains task index


output_path <- "/gpfs/project/zhao/zy92/GTEX/expr_normalized/" ## path for saving outputs of tasks
## e.g. Rscript --vanilla extract_normalized_expression.R chr_idx 

# exp files 
file_vec <- dir("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_expression_matrices",
               pattern = "*.bed$")
tissue_vec <- as.character(sapply(file_vec, function(x) unlist(strsplit(x, "\\."))[1]))
exp_list <- list()
for (i in 1:length(tissue_vec)) {
    exp_list[[tissue_vec[i]]] <- fread(paste0("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_expression_matrices/", 
                                            tissue_vec[i], 
                                           ".v8.normalized_expression.bed"))
}                                  

# extract the expression level and save in the directory
# loop over tissue
start_time <- Sys.time()
                               
for (i in 1) {   
    
    # create the chr path
    for (chr in c(1:22, "X", "Y")) {
        output_path_chr <- paste0(output_path, "chr", i, "/")
        if (!dir.exists(output_path_chr)) {
            dir.create(output_path_chr)
        }
    }
    
    # 
    tissue <- tissue_vec[i]
    tissue_exp <- exp_list[[tissue]]
    n_gene <- nrow(tissue_exp)
    id_vec <- colnames(tissue_exp)[5:ncol(tissue_exp)]
    colnames(tissue_exp)[1] <- "chr"
    for (j in 1:n_gene) {
        # info
        if (j %% 500 == 0) { 
            print(paste0("INFO: ", round(j*100/n_gene, 3), " % processed.")) 
        }
        chr <- tissue_exp[j, 1]
        gene <- tissue_exp[j, 4]
        tmp_df <- data.frame(id = id_vec, exp = as.numeric(tissue_exp[j, 5:ncol(tissue_exp)]))
    
        gene_dir <- paste0(output_path, "/", chr, "/", gene)
        tissue_file_name <- paste0(output_path, "/", chr, "/", gene, "/", tissue, ".txt")
        if (!dir.exists(gene_dir)) {
            dir.create(gene_dir)
        }
            
        write.table(tmp_df, file = tissue_file_name, row.names = F, col.names = F, quote = F)
    }


}
end_time <- Sys.time()

end_time - start_time