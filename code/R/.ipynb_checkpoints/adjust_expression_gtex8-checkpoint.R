##### adjust for covariates first #####
library(dplyr)
library(data.table)
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
for(chr in task_index) {
#for(chr in 16:22){
	glist = dir(paste0("/home/fas/radev/zy92/scratch60/GTEX/expr_gtex1/chr", chr));
	bgt = Sys.time()
	for(k in 1:length(glist)){
        print(paste0("INFO: gene", k ))
		tryCatch({
			g = glist[k];
			Yt = dir(paste0("/gpfs/loomis/project/radev/zy92/utmost_update/GTEX/expr_gtex1/chr", chr, "/", g, "/"))
			
			dir.create(paste0("/gpfs/loomis/project/radev/zy92/utmost_update/GTEX/adjusted_expr1/chr", chr), showWarnings = FALSE)
			dir.create(paste0("/gpfs/loomis/project/radev/zy92/utmost_update/GTEX/adjusted_expr1/chr", chr, "/", g), showWarnings = FALSE)
			setwd(paste0("/gpfs/loomis/project/radev/zy92/utmost_update/GTEX/adjusted_expr1/chr", chr, "/", g))
			## expr files ##
			Y = list()
			for(t in 1:length(Yt)){
				Y[[t]] = read.table(paste0("/gpfs/loomis/project/radev/zy92/utmost_update/GTEX/adjusted_expr1/chr", chr, "/", g, "/", Yt[t]), header=F)
			}
			ssize1 = unlist(lapply(Y, nrow))
			T_num = length(Yt)
			### load covariates ###
			helper <- function(x){
				unlist(strsplit(x, '-'))[2]
			}
			cov_fl = list.files("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/")
			idf = matrix(unlist(strsplit(cov_fl, "\\."))[1], ncol = 4, byrow = T)[,1]
			covar = list()
			for(t in 1:length(Yt)){
				itmp = cov_fl[which(idf==Yt[t])]
				covar[[t]] = read.table(paste0("/ysm-gpfs/pi/zhao-data/zy92/GTEx_V8/GTEx_Analysis_v8_eQTL_covariates/", itmp), header=F)
				tmp1 = sapply(Y[[t]][,1], helper)
				tmp2 = sapply(covar[[t]][1,-1], helper)
				ss_tmp = nrow(Y[[t]])
				cs_tmp = nrow(covar[[t]])
				idx1 = rep(0, ss_tmp)
				for(i in 1:ss_tmp){
					idx1[i] = which(tmp2 == tmp1[i])
				}
				tofactor = c()
				if(sum(tmp2[idx1] == tmp1)==ss_tmp){
					if(ss_tmp < 150){
						X = t(covar[[t]][2:19,-1])[idx1,]
						tofactor = c(which(covar[[t]][,1] == "gender"), which(covar[[t]][,1] == "Platform"))
						X = cbind(X,t(covar[[t]][tofactor,-1])[idx1,])
					}
					if(ss_tmp >= 150 & ss_tmp < 250){
						X = t(covar[[t]][2:34,-1])[idx1,]
						tofactor = c(which(covar[[t]][,1] == "gender"), which(covar[[t]][,1] == "Platform"))
						X = cbind(X,t(covar[[t]][tofactor,-1])[idx1,])
					}
					if(ss_tmp >= 250){
						X = t(covar[[t]][2:39,-1])[idx1,]
						tofactor = c(which(covar[[t]][,1] == "gender"), which(covar[[t]][,1] == "Platform"))
						X = cbind(X,t(covar[[t]][tofactor,-1])[idx1,])
					}
				}else{
					cat("shit happens at ", chr, "/", g, "/", Yt[t], "\n")
					break
				}
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
					if(length(rmv)){
						X = X[,-rmv]
					}
				}
				X = data.frame(y=Y[[t]][,2],X)
				adj_expr = resid(lm(y~., data=X))
				write.table(cbind(Y[[t]][,1], adj_expr), paste0(Yt[t], ".adj_expr"), quote=F, row.names=F, col.names=F)
			}
			ssize2 = unlist(lapply(covar, ncol)) - 1
			print(k)
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	print(Sys.time()-bgt)
}
