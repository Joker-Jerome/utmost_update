##### adjust for covariates first #####
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(fastmatch))
## For each tissue, normalized gene expression data was adjusted for covariates such as 
# gender, 
# sequencing platform, 
# pcr, 
# the top 3 principal components from genotype data 
# and top PEER Factors. 
# the number of PEER Factors used was determined by sample size: 15 for n < 150, 30 for n between 150 and 250, and 35 for n > 250 Covariate data was provided by GTEx.
# e. g. Rscript --vanilla adjust_expression_gtex8.R chr_idx

### load expr ###
options(stringsAsFactors=F)
# arguments
args = commandArgs(trailingOnly=TRUE)
chr_index = as.numeric(args[1]) ## chr index
netid = args[2]
input_expression_dir = "/gpfs/loomis/scratch60/zhao/zy92/GTEX/expr_gtex_lnc_updated/"
cov_dir = "/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/"
output_prefix = "/GTEX/adjusted_lnc_expr/" # /gpfs/loomis/scratch60/zhao/netid/prefix                    

###### decide which tissue to be predicted ######

### load covariates ###
helper <- function(x){
    unlist(strsplit(x, '-'))[2]
}               
cov_fl = list.files(cov_dir)
idf = matrix(unlist(strsplit(cov_fl, "\\.")), ncol = 4, byrow = T)[,1]
cov_tissue = as.character(sapply(cov_fl, function(x) unlist(strsplit(x, '\\.'))[1]))

                                 
# read in the covariates  
covar = list()
for (i in 1:length(cov_fl)){
	covar[[cov_tissue[i]]] = read.table(paste0(cov_dir, cov_fl[i]), header=F)
}
                                 
# read in the related transcripts                       
# Xt_extra_list_tissue                               
load("/gpfs/loomis/project/zhao/zy92/GTEX/expr_related_genes.RData") 

                                 
                                 
#tissues <- gsub("\\.v8\\.covariates\\.txt","",cov_fl)
### load expression and adjusting the covariates ###
for (chr in chr_index) {
    glist = dir(paste0(input_expression_dir, "chr", chr))
	bgt = Sys.time()
    #subjob_idx <- (length(glist)//10 * (task_idx-1)+1) : min((length(glist)//10 )* (task_idx, length(glist)))
	for (k in 1:length(glist)) {
        # monitoring process
        if (k %% 50 == 0 ) {
            print(paste0("INFO: ", round(k/length(glist) * 100, 2), "% processed" ))
        }
        # try catch block
		tryCatch({
			g = glist[k]
            print(paste0("INFO: gene ", g))
			Yt = dir(paste0(input_expression_dir, "chr", chr, "/", g, "/"))
			#
			dir.create(file.path(paste0("/gpfs/loomis/project/zhao/", netid, output_prefix, "chr", chr, "/", g)), recursive = TRUE)
			setwd(paste0("/gpfs/loomis/project/zhao/", netid, output_prefix, "chr", chr, "/", g))
			#print("finish reading the exp")
			## expr files ##
			Y = list()
           
			for(t in 1:length(Yt)){
				Y[[t]] = as.data.frame(fread(paste0(input_expression_dir, "chr", chr, "/", g, "/", Yt[t]), header=F))#, sep = "\t"))
				print(paste0("chr", chr, "/", g, "/", Yt[t]))
			}

			print("section B")			
			ssize1 = unlist(lapply(Y, nrow))
			T_num = length(Yt)
			
            Yt_tissue = as.character(sapply(Yt, function(x) unlist(strsplit(x, '\\.'))[1]))
			for(t in 1:length(Yt)){
				itmp = which(idf==Yt_tissue[t])
                # idx for y and for cov
                #print(head(Y[[t]]))
				tmp1 = sapply(Y[[t]][,1], helper)
				tmp2 = sapply(t(covar[[itmp]][1,-1]), helper)                
                ss_tmp = nrow(Y[[t]])
				cs_tmp = ncol(covar[[itmp]])-1
                idx1 = fmatch(tmp2, tmp1)
                tofactor = c()
                  
                # check the consistency 
				if(sum(tmp2 == tmp1[idx1])==cs_tmp){
                    # extract numeric variables for 4 classes with different number of covariates
					if(cs_tmp < 150){
						X = t(covar[[itmp]][2:20,-1])
						tofactor = c(which(covar[[itmp]][,1] == "sex"), which(covar[[itmp]][,1] == "pcr"), which(covar[[itmp]][,1] == "platform"))
						X = cbind(X,t(covar[[itmp]][tofactor,-1]))
					}
					if(cs_tmp >= 150 & cs_tmp < 250){
						X = t(covar[[itmp]][2:35,-1])
						tofactor = c(which(covar[[itmp]][,1] == "sex"), which(covar[[itmp]][,1] == "pcr"), which(covar[[itmp]][,1] == "platform"))
						X = cbind(X,t(covar[[itmp]][tofactor,-1]))
					}
					if(cs_tmp >= 250  & cs_tmp < 350){
						X = t(covar[[itmp]][2:50,-1])
						tofactor = c(which(covar[[itmp]][,1] == "sex"), which(covar[[itmp]][,1] == "pcr"), which(covar[[itmp]][,1] == "platform"))
						X = cbind(X,t(covar[[itmp]][tofactor,-1]))
					}
                    if(cs_tmp >= 350){
						X = t(covar[[itmp]][2:65,-1])
						tofactor = c(which(covar[[itmp]][,1] == "sex"), which(covar[[itmp]][,1] == "pcr"), which(covar[[itmp]][,1] == "platform"))
						X = cbind(X,t(covar[[itmp]][tofactor,-1]))
					}
				}else{
					cat("Error in sample numbers", chr, "/", g, "/", Yt[t], "\n")
					break
				}
                
                # convert the data type (numeric and factor type for regression analysis) 
				X = apply(X, 2, as.numeric)
				X = as.data.frame(X)
				if(length(tofactor)==1){
					ds_tmp = ncol(X)
					rmv = c()
					if(var(X[,ds_tmp])==0){
						rmv = c(rmv, ds_tmp)
					}else{
						X[,ds_tmp] = as.factor(X[,ds_tmp])
					}
				}else{
					ds_tmp = ncol(X)
					rmv = c()
					if(var(X[,ds_tmp])==0){
						rmv = c(rmv, ds_tmp)
					}else{
						X[,ds_tmp] = as.factor(X[,ds_tmp])
					}
					if(var(X[,ds_tmp-1])==0){
						rmv = c(rmv, ds_tmp-1)
					}else{
						X[,ds_tmp-1] = as.factor(X[,ds_tmp-1])
					}
                    if(var(X[,ds_tmp-2])==0){
						rmv = c(rmv, ds_tmp-2)
					}else{
						X[,ds_tmp-2] = as.factor(X[,ds_tmp-2])
					}
					if(length(rmv)){
						X = X[,-rmv]
					}
				}
                
                # add extra columns to adjust the expression
                Y[[t]] <- Y[[t]][idx1, ]
                X_extra <- Y[[t]] %>% 
                    left_join(Xt_extra_list_tissue[[Yt[t]]], by = c("V1" = "sample_id")) %>%
                    as.data.frame()
                
                # final adjustment model y ~ cov + extra_x 
				X = data.frame(y=Y[[t]][,2], X, X_extra[, 3:6]) 
				adj_expr = resid(lm(y~., data=X))
				print("INFO: before writing")
				write.table(cbind(Y[[t]][,1], adj_expr), paste0(Yt_tissue[t], ".adj_expr"), quote=F, row.names=F, col.names=F)
				print(paste0(Yt_tissue[t], ".adj_expr"))
			}
			ssize2 = unlist(lapply(covar, ncol)) - 1
			#print(k)
		}, error=function(e){
            cat("INFO: ERROR from TryCatch Block:",conditionMessage(e), "\n")
        }
        )
	}
	print(Sys.time()-bgt)
}
