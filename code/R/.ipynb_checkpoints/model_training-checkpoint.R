#-------------------------- arguments ---------------------------#
options(stringsAsFactors=F)
library(glmnet)
library(foreach)
library(genio)
library(Rcpp)



bgt_origin = Sys.time()

#-------------------------- functions ---------------------------#  

grad_prep <- function(X, Y){
	## pre-calculate some metrics for gradient
	ll = length(Y)
	P = ncol(X[[1]])
	XY = matrix(0,P,ll)
	for(i in 1:ll){
		XY[,i] = t(X[[i]])%*%Y[[i]]/nrow(X[[i]])
	}
	XY
}

cv_helper <- function(N, fold){
	valid_num = floor(N/fold)
	set.seed(123)
	perm = sample(1:N, size = N)
	idx1 = seq(1,N,valid_num)
	idx2 = c(idx1[-1]-1,N)
	list(perm=perm, idx=cbind(idx1,idx2))
}

minmax_lambda <- function(lst){
	max_lam = max(unlist(lapply(lst, function(x){max(x$lambda)})))
	min_lam = min(unlist(lapply(lst, function(x){min(x$lambda)})))
	c(min_lam, max_lam)
}

elastic_net_mse <- function(lst, X_tune, Y_tune, X_test, Y_test){
	P = length(lst)
	M = ncol(X_tune[[1]])
	lam_V = rep(0, P)
	test_res = list()
	test_beta = matrix(0, M, P)
	for(t in 1:P){
		ncv = length(lst[[t]]$lambda)
		tmp_mse = rep(0, ncv)
		for(k in 1:ncv){
			tmp_mse[k] = mean((Y_tune[[t]] - X_tune[[t]]%*%lst[[t]]$glmnet.fit$beta[,k])^2)
		}
		ss = which.min(tmp_mse)
		test_beta[,t] = lst[[t]]$glmnet.fit$beta[,ss]
		lam_V[t] = lst[[t]]$lambda[ss]
		predicted = X_test[[t]]%*%lst[[t]]$glmnet.fit$beta[,ss]
		test_res[[t]] = cbind(Y_test[[t]], predicted)
	}
	list(lam = lam_V, mse = test_res, est = test_beta)
}

multi_mse <- function(theta_est, X_test, Y_test){
	answer = list()
	P = ncol(theta_est)
	for(t in 1:P){
		predicted = X_test[[t]]%*%theta_est[,t]
		answer[[t]] = cbind(Y_test[[t]], predicted)
	}
	answer
}

avg_perm <- function(mse_lst){
	fd = length(mse_lst)
	P = length(mse_lst[[1]])
	rsq = mse = adj_mse = matrix(0, fd, P)
	for(f in 1:fd){
		for(t in 1:P){
			rsq[f,t] = (cor(mse_lst[[f]][[t]])[1,2])^2
			mse[f,t] = mean((mse_lst[[f]][[t]][,1]-mse_lst[[f]][[t]][,2])^2)
			adj_mse[f,t] = mse[f,t]/var(mse_lst[[f]][[t]][,1])
		}
	}
	cbind(apply(rsq, 2, mean), apply(mse, 2, mean), apply(adj_mse, 2, mean))

	#list(rsq = apply(rsq, 2, mean), mse = apply(mse, 2, mean), adj_mse = apply(adj_mse, 2, mean))
}

## group lasso on data with missing covariates, use validation data for stopping criterion
glasso <- function(X, Y, X1, Y1, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3, verbose = FALSE){
	Bgt = Sys.time()
	M = nrow(XY)
	P = length(X)
	NN = unlist(lapply(X, nrow))
	old_objV1 = rep(0,P)
	for(t in 1:P){
		old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
	}
  if (verbose) {
    cat("Training error: ", old_objV1, '\n')	
  }
	old_objV2 = rep(0,P)
	for(t in 1:P){
		old_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
	}
	if (verbose) {
    cat("Testing error: ", old_objV2, '\n')
  }
	beta_j_lasso = rep(0, P)
	tmp_XYj = 0
	if(!is.loaded("wrapper")){
		dyn.load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso.so") # change this to the abs path to optim.so
	}
	for(i in 1:maxiter){
		bgt = Sys.time()
		res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
		edt = Sys.time()
		
		new_objV1 = new_objV2 = rep(0,P)
		for(t in 1:P){
			new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
		}
		if (verbose) {print(paste0("Training error: ", new_objV1, '\n'))}
		for(t in 1:P){
			new_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
		}
		if (verbose) {print(paste0("Testing error: ", new_objV2, '\n'))}
		if(mean(new_objV2) > mean(old_objV2)|mean(new_objV1) > mean(old_objV1)){
			break
		}else{
			old_objV2 = new_objV2
		}
		if(max(abs(new_objV1-old_objV1)) < eps){
			break
		}else{
			old_objV1 = new_objV1
		}
	}
	Edt = Sys.time()
	cat("total training time: ", Edt-Bgt, "\n")
	list(est = theta, avg_tune_err = mean(new_objV2), tune_err=new_objV2)
}

## simpler version of glasso, train model until converges 
glasso_no_early_stopping <- function(X, Y, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3, verbose = FALSE){
  cat("running glasso_no_early_stopping\n")
	M = nrow(XY)
	P = length(X)
	NN = unlist(lapply(X, nrow))
	old_objV1 = rep(0,P)
	for(t in 1:P){
		old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
	}
	if (verbose) {cat("Training error: ", mean(old_objV1), '\n')}
	beta_j_lasso = rep(0, P)
	tmp_XYj = 0
	if(!is.loaded("wrapper")){
		dyn.load("optim.so")
	}
	for(i in 1:maxiter){
		res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
		new_objV1 = rep(0,P)
		for(t in 1:P){
			new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
		}
		if (verbose) {cat("Training error: ", mean(new_objV1), '\n')}
		if(max(abs(new_objV1-old_objV1)) < eps|mean(new_objV1) > mean(old_objV1)){
			break
		}else{
			old_objV1 = new_objV1
		}
	}
	list(est = theta, avg_train_err = mean(new_objV1), train_err = new_objV1)
}

# glasso final 
glasso_all <- function(X, Y, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3, verbose = FALSE){
	Bgt = Sys.time()
	M = nrow(XY)
	P = length(X)
	NN = unlist(lapply(X, nrow))
	old_objV1 = rep(0,P)
	for(t in 1:P){
		old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
	}
  if (verbose) {
    cat("Training error: ", old_objV1, '\n')	
  }
	
	beta_j_lasso = rep(0, P)
	tmp_XYj = 0
	if(!is.loaded("wrapper")){
		dyn.load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/glasso.so") # change this to the abs path to optim.so
	}
	for(i in 1:maxiter){
		bgt = Sys.time()
		res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
		edt = Sys.time()
		
		new_objV1 = rep(0,P)
		for(t in 1:P){
			new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
		}
		if (verbose) {print(paste0("Training error: ", new_objV1, '\n'))}
		
        # mean objective
		if(mean(new_objV1) > mean(old_objV1)){
			break
		}
        
        # update criteria
		if (max(abs(new_objV1-old_objV1)) < eps) {
			break
		} else {
			old_objV1 = new_objV1
		}
	}
	Edt = Sys.time()
	cat("total training time: ", Edt-Bgt, "\n")
	list(est = theta, avg_train_err = mean(new_objV1), train_err=new_objV1)
}

# extract genotype
extract_genotype <- function(gene_idx, chr_str, genotype_dir, gene_vec) {
    gene_dir <- paste0(genotype_dir, chr_str, gene_vec[gene_idx], "/", gene_vec[gene_idx])
    genotype_info <- read_plink(gene_dir)
    
    print(paste0("INFO: genotype matrix dimension:", dim(genotype_info$X)[1], " * ", dim(genotype_info$X)[2]))
    return(genotype_info)
}


#-------------------------- main ---------------------------#  

# testing
if (F) {
    args <- c("21", "1")
} 
args = commandArgs(trailingOnly = TRUE)
chr = as.numeric(args[1])
gene_idx = as.numeric(args[2])

#-------------------------- data import ---------------------------#  


k = "V1"
genotype_dir = "/ysm-gpfs/project/wl382/GTEx_v8/genotype/cis_loc/"
gtex_dir = "/gpfs/loomis/project/zhao/zy92/GTEX/" 
glist = dir(paste0(gtex_dir, "adjusted_expr/chr", chr))


# idx 
chr_str = paste0("chr", chr, "/")
gene_vec = list.files(paste0(genotype_dir, chr_str))
dose_idx = which(startsWith(gene_vec, glist[gene_idx]))
gene_id = gene_vec[dose_idx]
g = glist[gene_idx]
dose_path = paste0(genotype_dir, chr_str, gene_id, "/", gene_id)
Yt = dir(paste0(gtex_dir, "adjusted_expr/", chr_str, g, "/"))
P = length(Yt)
#output_dir = "/gpfs/loomis/project/zhao/zy92/GTEX/output/" 
output_dir = "/gpfs/scratch/zhao/zy92/GTEX/output_0715/"
ntune = 5
if_verbose = F
# loading fast matrix operations 
#sourceCpp("/gpfs/project/zhao/zy92/utmost_update/CTIMP/MatrixMtp.cpp") # call the C++ file and we have three functions as armaMatMult，eigenMatMult，eigenMapMatMult.
sourceCpp("/gpfs/loomis/project/zhao/zy92/utmost_update/CTIMP/MatrixMtp.cpp") # call the C++ file and we have three functions as armaMatMult，eigenMatMult，eigenMapMatMult.



# create dirs
dir.create(paste0(output_dir, "chr", chr), showWarnings = FALSE)
dir.create(paste0(output_dir, "chr", chr, "/", gene_id), showWarnings = FALSE)
setwd(paste0(output_dir, "chr", chr, "/", gene_id))

# test existing files
est_file = paste0(k, '.', gene_id, ".est")
if (file.exists(est_file)) {
    stop("INFO: Est file exists")
    return(0)
}

# expr files 
print("INFO: loading expression files...")
Y = list()
T_idx = c()
for(t in 1:P){
    tmp_exp = read.table(paste0(gtex_dir, "adjusted_expr/chr", chr, "/", g, "/", Yt[t]), header=F)
    # check if the y is constant
    if (!all(tmp_exp[,2] == tmp_exp[1,2])) {
        T_idx = c(T_idx, t)
        Y[[length(T_idx)]] = tmp_exp 
    }
     
}
Yt = Yt[T_idx]
ssize = unlist(lapply(Y, nrow))
T_num = length(Yt)
P = T_num

# genotype files 
print("INFO: loading genotype files...")	
genotype_info <- read_plink(dose_path)
dose = genotype_info$X #[1:1000, ]
dose_std = dose

# center the dose
for(j in 1:nrow(dose)){
    dose_std[j, is.na(dose[j, ])] <- mean(dose[j, ], na.rm = T)
    dose_std[j, ] <- dose_std[j, ] - mean(dose[j, ], na.rm = T)
}

# covariance matrix 
n_sample = ncol(dose_std)
tmp = t(as.matrix(dose_std))
XX = eigenMatMult(t(tmp), tmp)/n_sample 
Xnorm = diag(XX)
remove(tmp)
remove(XX)
sub_id = colnames(dose_std)

# number of snps
M = nrow(dose_std)
sub_id_map = list()
sub_id_map_exp = list()

# sample matching 
for(t in 1:T_num){
    #tmp = rep(0, nrow(Y[[t]]))
    exp_id <- as.character(sapply(Y[[t]][,1], function(x) substr(x, 1, 10)))
    # index of ids that have matched genotypes 
    # exp_id based on order of exp_id
    match_id <- match(exp_id, sub_id)
    sub_id_map[[t]] <- na.omit(match_id)
    sub_id_map_exp[[t]] <- !is.na(match_id)  
}

# cv                                   
# feature selection based on glmnet
X_all = list()
Y_all = list()
glmnet_list = list()
for(t in 1:T_num){
    X_all[[t]] = apply(as.matrix(dose_std[, sub_id_map[[t]]]),1,as.numeric)
    Y_all[[t]] = Y[[t]][which(sub_id_map_exp[[t]] == T), 2]
}
# initial values 
                                  
#alpha_vec = seq(0, 1, 0.1)
alpha_vec = c(1)
single_initial_est = matrix(0, ncol(X_all[[1]]), T_num)
for(t in 1:T_num) {
    print(paste0("INFO: Tissue ", t))
    cvm_min = c()
    for (alpha in alpha_vec) {
        set.seed(1)
        tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = alpha, nfolds = 5)
        cvm_min = c(cvm_min, min(tt$cvm))
        glmnet_list[[t]] = tt
        # update the beta if the current cvm is the min
        if (min(tt$cvm) == min(cvm_min)) {
            single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
        }
        
    }
}

# feature selection results
initial_numeric_df = data.frame(matrix(single_initial_est, nrow= 49))
zero_summary = apply(initial_numeric_df, 2, function(x) sum(x == 0))
print(paste0("INFO: Number of SNPs ", sum(zero_summary != T_num)))

# subset index
index_snp = which(zero_summary != T_num)
M = length(index_snp)

# subsetting the genotype 
dose = genotype_info$X[index_snp, ] #[1:1000, ]
dose_std = dose

# center the dose
for(j in 1:nrow(dose)){
    dose_std[j, is.na(dose[j, ])] <- mean(dose[j, ], na.rm = T)
    dose_std[j, ] <- dose_std[j, ] - mean(dose[j, ], na.rm = T)
}

# covariance matrix 
n_sample = ncol(dose_std)
tmp = t(as.matrix(dose_std))
XX = eigenMatMult(t(tmp), tmp)/n_sample 
Xnorm = diag(XX)
remove(tmp)
remove(XX)
sub_id = colnames(dose_std)

# tuning based on CV
fold = 5
#IR                                
print("INFO: CV preparation")
cv_config = cv_helper(n_sample, fold)
cv_perm = cv_config$perm
cv_idx = cv_config$idx

single_res_test = list()
single_lam = matrix(0, fold, T_num)
single_theta_est = list()
	
multi_res_test = list()
multi_lam = matrix(0, fold, 2)
multi_theta_est = list()

multi_res_test2 = list()
multi_lam2 = array(0, dim=c(fold, T_num, 2))
multi_theta_est2 = list()
                                  

res_tune = list()
rec_lamv = matrix(0, fold, ntune)                                 
                                  
for(f in 1:fold){
#for(f in 1:1){
    bgt = Sys.time()
    print(paste0("INFO: fold ", f))
    test_index = cv_perm[cv_idx[f,1]:cv_idx[f,2]]
    test_id = sub_id[test_index]

    # move the tuning idx to another idx block
    tuning_index = cv_perm[cv_idx[f%%fold+1,1]:cv_idx[f%%fold+1,2]]
    tuning_id = sub_id[tuning_index]

    # idx list
    X_test = list()
    Y_test = list()
    X_tune = list()
    Y_tune = list()
    X_train = list()
    Y_train = list()
    for(t in 1:T_num){
        # idx matching  
        id_filter = !(sub_id_map[[t]]%in%c(test_index,tuning_index))
        X_train_tmp = sub_id_map[[t]][id_filter]
        Y_train_tmp = which(sub_id_map_exp[[t]] == T)[id_filter]
        tuning_id_filter = sub_id_map[[t]]%in%tuning_index
        X_tuning_tmp = sub_id_map[[t]][tuning_id_filter]
        Y_tuning_tmp = which(sub_id_map_exp[[t]] == T)[tuning_id_filter]
        testing_id_filter = sub_id_map[[t]]%in%test_index
        X_test_tmp = sub_id_map[[t]][testing_id_filter]
        Y_test_tmp = which(sub_id_map_exp[[t]] == T)[testing_id_filter]
        # training data
        X_train[[t]] = apply(as.matrix(dose_std[,X_train_tmp]),1,as.numeric)
        Y_train[[t]] = Y[[t]][Y_train_tmp, 2]
        X_tune[[t]] = apply(as.matrix(dose_std[,X_tuning_tmp]),1,as.numeric)
        Y_tune[[t]] = Y[[t]][Y_tuning_tmp, 2]
        X_test[[t]] = apply(as.matrix(dose_std[,X_test_tmp]),1,as.numeric)
        Y_test[[t]] = Y[[t]][Y_test_tmp, 2]
    }

    ## model training ##	
    ## train elastic net and used average lambda as tuning parameters ##
    single_initial_est = matrix(0, ncol(X_train[[1]]), T_num)
    single_summary = list()
    for(t in 1:T_num) {
        #print(Y_train[[t]])
        #print(sum(is.na(X_train[[t]])))
        if (t %% 5 == 0) {
            print(paste0("INFO: glmnet cv tissue", t))
        }
        
        
        tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5)
        single_summary[[t]] = tt
        single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
    }
    ## performance of Elastic net on tuning and testing data with various tuning parameters
    els_output = elastic_net_mse(single_summary, X_tune, Y_tune, X_test, Y_test)
    single_res_test[[f]] = els_output$mse
    single_lam[f,] = els_output$lam
    single_theta_est[[f]] = els_output$est
    remove(els_output)

    
    ## use elastic net ests row norm as weights ##
    lam_range = minmax_lambda(single_summary)
    sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
    sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
    sig_norm = sig_norm/sum(sig_norm)
    weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);

    tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
    tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
    tis_norm = tis_norm/sum(tis_norm)
    weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
    lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)
    
    # initial value
    initial_numeric = as.numeric(single_initial_est)
	remove(single_summary)
    remove(single_initial_est)


    #-------------------------- train - validate - test ---------------------------#
    print('INFO: starting train-validate-test\n')
    ## preparation
    XY = grad_prep(X_train, Y_train)
    #bgt = Sys.time()
    XX_train = lapply(X_train, function(x) { eigenMatMult(t(x), x)/nrow(x)})
    #edt = Sys.time()
    #print(edt-bgt)
    
    spsz = unlist(lapply(X_train,nrow))
    #res_tune = rep(0, ntune)
    res_tune[[f]] = array(-1, dim=c(ntune, ntune, T_num))
    #best.lam = 0
    rec_lamv[f,] = lam_V
    for(lam1 in 1:ntune){
        for(lam2 in 1:ntune){
            single_est = matrix(initial_numeric, M, T_num )
            ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est, verbose = if_verbose)
            if(sum(ans$est!=0)>0){
                res_tune[[f]][lam1,lam2, ] = ans$tune_err
                if (if_verbose) { cat("lambda1=",lam_V[lam1], "; lambda2=", lam_V[lam2], "; avg tune err=", ans$avg_tune_err, '\n') }
                remove(single_est); remove(ans);
            }else{
                remove(single_est); remove(ans);
                break
            }			
        }
    }
    
    
    #-------------------------- save results on test set for evaluation ---------------------------#
    print("INFO: saving cross-validation results for evaluation and analysis")
    avg_tune_res = apply(res_tune[[f]], c(1,2), mean)
    best.lam = which(avg_tune_res == min(avg_tune_res[avg_tune_res>=0]), arr.ind = TRUE)[1,]
    single_est = matrix(initial_numeric, M, P)
    ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_est, verbose = if_verbose)
    multi_res_test[[f]] = multi_mse(ans$est, X_test, Y_test)
    multi_lam[f,] = lam_V[best.lam]
    multi_theta_est[[f]] = ans$est
    remove(single_est)
    remove(ans)

    
    edt = Sys.time()
    print(edt-bgt)
} # end of cv loop
edt_origin = Sys.time() 
print("Final Time Cost: ")                                 
print(edt_origin-bgt_origin)
   
bgt = Sys.time()
X_all = list()
Y_all = list()
for(t in 1:T_num){
    X_all_tmp = sub_id_map[[t]]
    X_all[[t]] = apply(as.matrix(dose_std[,X_all_tmp]), 1, as.numeric)
    Y_all[[t]] = Y[[t]][which(sub_id_map_exp[[t]] == T), 2]
}
# initial values 
single_initial_est = matrix(0, ncol(X_train[[1]]), T_num)
for(t in 1:T_num){
    print(paste0("INFO: Tissue ", t))
    tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5)
    single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
}

sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
sig_norm = sig_norm/sum(sig_norm)
weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);

tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
tis_norm = tis_norm/sum(tis_norm)
weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);

spsz = unlist(lapply(X_all,nrow))
initial_numeric = as.numeric(single_initial_est)
#remove(single_initial_est)
XY = grad_prep(X_all, Y_all)

XX_all = lapply(X_all, function(x){eigenMatMult(t(x), x)/nrow(x)})
tmp_res = rep(0, fold)
for(f in 1:fold){


    ans = glasso_all(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=multi_lam[f,1]/spsz, lambda2=multi_lam[f,2], theta=matrix(initial_numeric,M,P), verbose = if_verbose)
    tmp_res[f] = ans$avg_train_err
}
final.lam = multi_lam[which.min(tmp_res),]
ans = glasso_all(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=final.lam[1]/spsz, lambda2=final.lam[2], theta=matrix(initial_numeric,M,P), verbose = if_verbose)
#info = read.table(info_path, header=T, sep='\t')
#downstream_est = data.frame(info[,1:3], ans$est)
multi_all_res = multi_mse(ans$est, X_all, Y_all)
single_all_res = multi_mse(single_initial_est, X_all, Y_all)

info = genotype_info$bim
downstream_est = data.frame(info[index_snp, ], ans$est)
edt = Sys.time()
print(edt-bgt)

# save the results
save(single_res_test, single_lam, single_theta_est, multi_res_test, multi_lam, multi_theta_est, index_snp, file = paste0(k, 'chr', chr, '.', gene_id, ".RData"))
write.table(downstream_est, paste0(k, '.', gene_id, ".est"), quote=F, row.names=F, col.names=c(colnames(genotype_info$bim), Yt))
print('done!')