args <- commandArgs(trailingOnly=TRUE)
idx <- as.numeric(args[1])

# library
suppressPackageStartupMessages(library(formattable))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))

predictDB_dir <- "/ysm-gpfs/pi/zhao-data/zy92/predictDB/"
db_vec <- list.files(predictDB_dir, pattern = "*.db")

## convenient query function
sqlite <- dbDriver("SQLite")
query <- function(...) dbGetQuery(db, ...)

setwd(predictDB_dir)
gene_list <- list()
for(i in 1:length(db_vec))
  {
    IRdisplay::display_html(paste0("INFO: reading ", db_vec[i]))
    dbname <- db_vec[i]
    db <- dbConnect(sqlite, dbname)
    #print(query('select count(*) from weights'))
    gene_list[[dbname]] <- query(paste('select distinct(gene) from weights'))
    
}

# gene set 
gene_set <- unique(unlist(gene_list))
length(gene_set)


# check function
check_weights <- function(gene)  { 
    snp_total <- c()
    snp_gene <- list()
    weights_gene <- list()
    #IRdisplay::display_html(paste0("INFO: reading snp info"))
    for(i in 1:length(db_vec))
      { 
        dbname <- db_vec[i]
        db <- dbConnect(sqlite,dbname)
        query <- function(...) dbGetQuery(db, ...)
        #print(query('select count(*) from weights'))
        snp_cur <- query(paste0('select * from weights where gene = "', gene, '"', sep = ""))$rsid
        snp_gene[[dbname]] <- snp_cur
        snp_total <- c(snp_total, snp_cur)
        weights_gene[[dbname]] <- query(paste('select * from weights where gene = "', gene, '"', sep = ""))$weight
      }
    
    snp_total <- unique(snp_total)
    weights <- matrix(0, length(snp_total), length(db_vec))
    #IRdisplay::display_html(paste0("INFO: finish loading snp info"))
    for(i in 1:length(db_vec))
      { 
        total_snp <- c()
        #IRdisplay::display_html(paste0("INFO: reading ", db_vec[i]))
        dbname <- db_vec[i]
        #db <- dbConnect(sqlite,dbname)
        #print(query('select count(*) from weights'))
        snp_cur <- snp_gene[[dbname]]
        index <- match(snp_cur, snp_total)
        indi <- is.na(index) == FALSE
        index <- match(snp_cur, snp_total)
        indi <- !is.na(index)
        weights[index[indi],i] <- weights_gene[[dbname]][indi]
    }
    #return(list(output1 = weights, output2 = weights_gene))  
    return(weights)
}

summary_list <- list()
# function 
analyze_weights <- function(gene, summary_list) {
    test_weights <- check_weights(gene)
    n_snp <- nrow(test_weights)
    n_tissue <- ncol(test_weights)
    
    snp_vec <- c()
    for (i in 1:n_snp) {
        vec_sign <- sign(test_weights[i,])
        res_positive <- length(which(vec_sign > 0))
        res_negative <- length(which(vec_sign < 0))
        # if the vector contains both postive values and negative values
        if (res_positive > 0 & res_negative > 0) {
             snp_vec <- c(snp_vec, i)
        }
    }
    output_list <- list(snp_vec = snp_vec, num = length(snp_vec)) 
    return(output_list)    
}

start_idx <- 3000*(idx-1) + 1
end_idx <- min(3000*idx, length(gene_set))

for (i in start_idx:end_idx) {
    if (i %% 500 == 0 ) { IRdisplay::display_html(paste0("INFO: gene", i)) }
    gene <- gene_set[i]
    summary_list[[gene]] <- analyze_weights(gene)
}

file_nam <- paste0("/ysm-gpfs/pi/zhao/zy92/projects/utmost_update/utmost_update/data/snp_direction_summary_", idx, ".RData")
save(summary_list, file = file_nam)
