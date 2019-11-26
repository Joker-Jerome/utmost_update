##### adjust for covariates first #####
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
## For each tissue, normalized gene expression data was adjusted for covariates such as 
# gender, 
# sequencing platform, 
# the top 3 principal components from genotype data 
# and top PEER Factors. 
# the number of PEER Factors used was determined by sample size: 15 for n < 150, 30 for n between 150 and 250, and 35 for n > 250 Covariate data was provided by GTEx.
# e. g. Rscript --vanilla adjust_expression_gtex8.R chr_idx

### load expr ###
options(stringsAsFactors=F)
# arguments
args = commandArgs(trailingOnly=TRUE)
task_index = as.numeric(args[1]) ## a file contains task index


###### decide which tissue to be predicted ######

### load covariates ###
helper <- function(x){
    unlist(strsplit(x, '-'))[2]
}               
cov_fl = list.files("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/")
idf = matrix(unlist(strsplit(cov_fl, "\\.")), ncol = 4, byrow = T)[,1]
print(paste0("INFO: idf ", idf))
cov_tissue = as.character(sapply(cov_fl, function(x) unlist(strsplit(x, '\\.'))[1]))

# read in the covariates  
covar = list()
for (i in 1:length(cov_fl)){
	covar[[cov_tissue[i]]] = read.table(paste0("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/", cov_fl[i]), header=F)
}

### load expression and adjusting the covariates ###
for (chr in task_index) {
	glist = dir(paste0("/gpfs/loomis/scratch60/radev/zy92/GTEX/expr_gtex1/chr", chr))
	bgt = Sys.time()
	for(k in 1:length(glist)) {
        # print(paste0("INFO: gene", k ))
        # monitoring process
        if (k %% 500 == 0 ) {
            print(paste0("INFO: ", round(k/length(glist)) * 100), "% processed" )
        }
        # try catch block
		tryCatch({
			g = glist[k]
            #print(g)
			Yt = dir(paste0("/gpfs/loomis/scratch60/radev/zy92/GTEX/expr_gtex1/chr", chr, "/", g, "/"))
			dir.create(paste0("/gpfs/loomis/scratch60/radev/zy92/GTEX/adjusted_expr1/chr", chr), showWarnings = FALSE)
			dir.create(paste0("/gpfs/loomis/scratch60/radev/zy92/GTEX/adjusted_expr1/chr", chr, "/", g), showWarnings = FALSE)
			setwd(paste0("/gpfs/loomis/scratch60/radev/zy92/GTEX/adjusted_expr1/chr", chr, "/", g))
			## expr files ##
			Y = list()
			for(t in 1:length(Yt)){
                #print(paste0("INFO: Yt", t ))
				Y[[t]] = read.table(paste0("/gpfs/loomis/scratch60/radev/zy92/GTEX/expr_gtex1/chr", chr, "/", g, "/", Yt[t]), header=F, sep = " ")
			}
			ssize1 = unlist(lapply(Y, nrow))
			T_num = length(Yt)
			
            #print(Yt)
            Yt_tissue = as.character(sapply(Yt, function(x) unlist(strsplit(x, '\\.'))[1]))
            #print(Yt_tissue)                                    
			for(t in 1:length(Yt)){
				itmp = cov_fl[which(idf==Yt_tissue[t])]
                #print(paste0("INFO: itmp ", itmp))
                # idx for y and for cov
				tmp1 = sapply(Y[[t]][,1], helper)
				tmp2 = sapply(covar[[t]][1,-1], helper)
				ss_tmp = nrow(Y[[t]])
				cs_tmp = ncol(covar[[t]])-1
				#idx1 = rep(0, ss_tmp)
                #print(tmp1[8:18])
                #print(tmp2[8:18])
                #print(length(tmp1))
                #print(length(tmp2))
                #print(length(idx1))
				#for(i in 1:ss_tmp){
                #    print(i)
				#	idx1[i] = which(tmp2 == tmp1[i])
				#}
                # matching the samples
                #idx1 = na.omit(match(tmp1, tmp2))
                idx1 = match(tmp2, tmp1)
                
                #print(sum(is.na(idx1)))
				tofactor = c()
                
                #print(paste0("INFO: number in common: ", sum(tmp2 == tmp1[idx1])))
                #print(paste0("INFO: cs_tmp ", cs_tmp))
                # check the 
				if(sum(tmp2 == tmp1[idx1])==cs_tmp){
                    print(paste0("INFO: cs_tmp ", cs_tmp))
                    
                    # extract numeric variables 4 classes with different number of covariates
					if(cs_tmp < 150){
						X = t(covar[[t]][2:20,-1])
						tofactor = c(which(covar[[t]][,1] == "sex"), which(covar[[t]][,1] == "pcr"), which(covar[[t]][,1] == "platform"))
						X = cbind(X,t(covar[[t]][tofactor,-1]))
					}
					if(cs_tmp >= 150 & cs_tmp < 250){
						X = t(covar[[t]][2:35,-1])
						tofactor = c(which(covar[[t]][,1] == "sex"), which(covar[[t]][,1] == "pcr"), which(covar[[t]][,1] == "platform"))
						X = cbind(X,t(covar[[t]][tofactor,-1]))
					}
					if(cs_tmp >= 250  & cs_tmp < 350){
						X = t(covar[[t]][2:50,-1])
						tofactor = c(which(covar[[t]][,1] == "sex"), which(covar[[t]][,1] == "pcr"), which(covar[[t]][,1] == "platform"))
						X = cbind(X,t(covar[[t]][tofactor,-1]))
					}
                    if(cs_tmp >= 350){
						X = t(covar[[t]][2:65,-1])
						tofactor = c(which(covar[[t]][,1] == "sex"), which(covar[[t]][,1] == "pcr"), which(covar[[t]][,1] == "platform"))
						X = cbind(X,t(covar[[t]][tofactor,-1]))
					}
				}else{
					cat("Error in sample numbers", chr, "/", g, "/", Yt[t], "\n")
					break
				}
                print("INFO: X matrix")
                
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
                Y[[t]] <- Y[[t]][idx1, ]
                print(paste0("dim y :", nrow(Y[[t]])))
                print(paste0("dim x :", nrow(X)))
				X = data.frame(y=Y[[t]][,2],X)
				adj_expr = resid(lm(y~., data=X))
                #print(adj_expr)
				write.table(cbind(Y[[t]][,1], adj_expr), paste0(Yt_tissue[t], ".adj_expr"), quote=F, row.names=F, col.names=F)
			}
			ssize2 = unlist(lapply(covar, ncol)) - 1
			print(k)
		}, error=function(e)
        {cat("ERROR from TryCatch Block:",conditionMessage(e), "\n")}
        )
	}
	print(Sys.time()-bgt)
}
