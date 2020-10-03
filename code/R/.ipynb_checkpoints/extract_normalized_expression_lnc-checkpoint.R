# libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

# arguments
args = commandArgs(trailingOnly=TRUE)
task_index = as.numeric(args[1]) ## a file contains task index
output_path <- "/gpfs/scratch60/zhao/zy92/GTEX/expr_gtex_lnc/" ## path for saving outputs of tasks
## e.g. Rscript --vanilla extract_raw_expression.R chr_idx 

# extract the expression level and save in the directory
# loop over tissue
start_time <- Sys.time()
                          

# load data
load("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/processed_data/exp_df_various_biotypes_various_tissues.RData")
load("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/processed_data/annotation_complete_chr.RData")

tissue_vec_target <- target_tissue_vec <- names(df_list_lnc)[1:49]
# extract the expression level and save in the directory
# loop over chr
#for (i in 1:22) {
for (i in task_index) {   
    
    # create the chr path
    print(paste0("INFO: chr ", i ))
    output_path_chr <- paste0(output_path, "chr", i, "/")
    if (!dir.exists(output_path_chr)) {
        dir.create(output_path_chr, recursive = T)
    }
    
    # subset the specific gene set for the transcript type
    transcriptome_summary <- res_summary %>%
        filter(ensembl_id %in% df_list_lnc[[1]]$Name)
    expression_info_tmp <- transcriptome_summary %>% 
        filter(chromosome_name %in% UQ(i))
    head(expression_info_tmp)
    
    # create the gene dir
    for (gene in expression_info_tmp$ensembl_id) {
        gene_dir <- paste0(output_path_chr, gene)
        if (!dir.exists(gene_dir)) {
            dir.create(gene_dir, recursive = T)
        }
    }
    
    # tissue-based extraction
    for (tissue in tissue_vec_target) {
        tmp_exp <- df_list_lnc[[tissue]] %>%
            filter(Name %in% expression_info_tmp$ensembl_id) %>%
            as.data.frame()
        sample_id <- as.character(sapply(colnames(tmp_exp)[3:ncol(tmp_exp)], function(x) 
                paste0(unlist(strsplit(x, split = "-"))[1:2], collapse = "-"))
        )
        tmp_gene_vec <- tmp_exp$Name
        for (j in 1:length(tmp_gene_vec) ) {
            if (j %% 100 == 0) { 
                IRdisplay::display_html(paste0("INFO: ", round(j*100/length(tmp_gene_vec), 3), " % processed.")) 
            }
            gene <- tmp_gene_vec[j]
            tmp_file <- paste0(output_path_chr, gene, "/", tissue, ".txt")
            exp_vec <- as.numeric(tmp_exp[j, 3:ncol(tmp_exp)])
            
            # filter out the transcripts with too many zeros
            n_obs <- length(exp_vec)
            if (sum(exp_vec == 0) > round(0.95 * n_obs)) {
                next
            }
            tmp_df <- data.frame(sample_id = sample_id, exp = exp_vec)
            IRdisplay::display_html(paste0("INFO: ", gene))
            write.table(tmp_df, file = tmp_file, row.names = F, col.names = F, quote = F)
            
        }
       
    }
}
    
   
     

end_time <- Sys.time()

# running time
end_time - start_time
