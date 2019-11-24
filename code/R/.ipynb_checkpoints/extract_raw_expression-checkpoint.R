# libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

# arguments
args = commandArgs(trailingOnly=TRUE)
task_index = args[1] ## a file contains task index
output_path <- "/gpfs/loomis/project/radev/zy92/utmost_update/GTEX/expr_gtex1/" ## path for saving outputs of tasks
sqfile_path = args[3] ## path for saving simple queue files
## e.g. Rscript --vanilla extract_raw_expression.R chr_idx 

load("/gpfs/loomis/project/radev/zy92/utmost_update/GTEX/raw_expression.RData")
# extract the expression level and save in the directory
# loop over chr
#for (i in 1:22) {
for (i in task_index) {    
    # create the chr path
    print(paste0("INFO: chr ", i ))
    output_path_chr <- paste0(output_path, "chr", i, "/")
    if (!dir.exists(output_path_chr)) {
        dir.create(output_path_chr)
    }
    # subset the gene for the chr
    expression_info_tmp <- gene_info_target %>% 
        filter(chromosome_name %in% UQ(i))
    head(expression_info_tmp)
    gene_vec <- expression_info_tmp$ensembl_gene_id
    
    expression_list_target_chr <- list()
    n_sample <- list()
    #expression_list_target_chr[[i]] <- list()
    for (tissue in tissue_vec_target) {
        tmp_exp <- expression_info_tmp %>%
            left_join(expression_list_target[[tissue]], by = c("ensembl_gene_id" = "ensembl_gene_id"))
        expression_list_target_chr[[tissue]] <- tmp_exp
        n_sample[[tissue]] <- ncol(tmp_exp) - 9
        #print(n_sample[[tissue]])
    }
    
    # tissue-specific sample id 
    id_vec <- list()
    for (tissue in tissue_vec_target) {
        id_vec[[tissue]] <- colnames(expression_list_target_chr[[tissue]])[8:(7+as.numeric(n_sample[[tissue]]))]

    }
    # loop over gene
    #gene2idx <- list()
    #idx <- match(gene, )
    for (j in 1:length(gene_vec)) {
    #for (j in 1:1) {
        gene <- gene_vec[j]
        # watch the process
        if (j %% 100 == 0) { 
            print(paste0("INFO: ", round(j*100/length(gene_vec), 3), " % processed.")) 
        }
        output_path_chr_gene <- paste0(output_path_chr, gene)
        if (!dir.exists(output_path_chr_gene)) {
            dir.create(output_path_chr_gene)
        } 
        
        # loop over tissue
        for (tissue in tissue_vec_target) {
            tissue_file_name <- paste0(output_path_chr_gene, "/", tissue, ".txt")
            #exp_vec <- expression_list_target_chr[[tissue]] %>% 
                #filter(ensembl_gene_id == gene)
            exp_vec <- expression_list_target_chr[[tissue]][j,]
            #print(dim(exp_vec))
            #print(n_sample[[tissue]])
            #print(exp_vec[, 1:20])
            exp_vec <- as.numeric(exp_vec[,8:(7+as.numeric(n_sample[[tissue]]))])
            #print(head(exp_vec))
            #id_vec <- colnames(expression_list_target_chr[[tissue]])[8:(7+as.numeric(n_sample[[tissue]]))]
            output_df <- data.frame(id = id_vec[[tissue]], exp = exp_vec)
            write.table(output_df, file = tissue_file_name, row.names = F, col.names = F, quote = F)
        }
        
    }
    
}


