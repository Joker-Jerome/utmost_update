### chr:gene ###
options(stringsAsFactors=F)
library(pscl)
rb <- function(n, prob){
	as.numeric(runif(n) < prob)
}
gibbs <- function(initial, burnin, total, X, Y, a0, b0, alp0, bet0)
{
	library(MASS); library(pscl);
	TT = length(Y)
	XX = lapply(X, function(x){x%*%t(x)})
	s0_current = initial$s0_0; se_current = initial$se_0; sg_current = initial$sg_0; p_current = initial$p_0; Z_current = initial$Z_0; b_current = initial$b_0; B_current = initial$B_0;
	M = nrow(B_current)
	NN = unlist(lapply(Y, length))
	XY = list()
	for(t in 1:TT){
		XY[[t]] = as.numeric(t(X[[t]])%*%Y[[t]])
	}
	s0_all = 0; se_all = rep(0, TT); sg_all = rep(0, TT); p_all = 0; B_all = matrix(0, M, TT)
#	pb = txtProgressBar(0, total, style=3)
	for(k in 1:total){
#		setTxtProgressBar(pb, k)
		## beta ##
		for(t in 1:TT){
			s = se_current[t]*diag(NN[t]) + sg_current[t]*XX[[t]]
			s_inv = chol2inv((chol(s)))
			Sig = sg_current[t]*(diag(M) - sg_current[t]*t(X[[t]])%*%s_inv%*%X[[t]])
			V = XY[[t]]/se_current[t] + Z_current[,t]*b_current/sg_current[t]
			A = t(chol(Sig))
        	B_current[ ,t] = Sig%*%V + A%*%rnorm(M)
		}
		## Z ##
		for(t in 1:TT){
			tmp1 = exp(-(B_current[,t]-b_current)^2/(2*sg_current[t]))
			tmp2 = exp(-(B_current[,t])^2/(2*sg_current[t]))
			Z_current[, t] = rb(M, p_current*tmp1/(p_current*tmp1 + (1-p_current)*tmp2))
		}
		## b ##
		if(TT>1){
			V1 = 1/(colSums(apply(Z_current, 1, function(x, a){x/a}, a=sg_current)) + 1/s0_current)
			V2 = colSums(apply(Z_current*B_current, 1, function(x, a){x/a}, a=sg_current))
		}else{
			V1 = 1/(apply(Z_current, 1, function(x, a){x/a}, a=sg_current) + 1/s0_current)
			V2 = apply(Z_current*B_current, 1, function(x, a){x/a}, a=sg_current)
		}
		b_current = rnorm(M, mean = V1*V2, sd = sqrt(V1))
		## p ##
		SZ = sum(Z_current)
		p_current = rbeta(1, a0+SZ, b0+TT*M-SZ)
		## s0, se, sg ##
		s0_current = rigamma(1, alp0+M/2, bet0 + sum(b_current^2)/2)
		for(t in 1:TT){
			se_current[t] = rigamma(1, alp0 + NN[t]/2, bet0 + as.numeric(sum((Y[[t]]-X[[t]]%*%B_current[,t])^2)))
			sg_current[t] = rigamma(1, alp0 + M/2, bet0 + as.numeric(sum((B_current[,t]-Z_current[,t]*b_current)^2)/2))
		}
		## output difference to true values
#		mean(abs(B_current - B)), mean(abs(s0_current - s0)), mean(abs(sg_current - sg)), mean(abs(se_current - se)), abs(s0_current - s0), abs(p_current - p)
#		if(TT>1){
#			cat("Mean abs err of ", k, "th iter -- B:", mean(abs(B_current - B)), "sg:", mean(abs(sg_current - sg)), "se: ", mean(abs(se_current - se)), "s0: ", abs(s0_current - s0), "p: ", abs(p_current - p), "\n")
#		}else{
#			cat("Mean abs err of ", k, "th iter -- B:", mean(abs(B_current - B[,pred_index])), "sg:", mean(abs(sg_current - sg[pred_index])), "se: ", mean(abs(se_current - se[pred_index])), "s0: ", abs(s0_current - s0), "p: ", abs(p_current - p), "\n")			
#		}
		if(k > burnin){
			s0_all = s0_all + s0_current
			se_all = se_all + se_current
			sg_all = sg_all + sg_current
			p_all = p_all + p_current
			B_all = B_all + B_current
		}
	}
	list(s0_avg = s0_all/(total-burnin), se_avg = se_all/(total-burnin), sg_avg = sg_all/(total-burnin), p_avg = p_all/(total-burnin), B_avg = B_all/(total-burnin))
}

###### decide which tissue to be predicted ######
gid_fl = list.files("/net/zhao/yh367/GTEX/GTEx_Analysis_v6p_all-associations/gene_id")
T_names = matrix(unlist(strsplit(gid_fl[1:44], "_Analysis")), ncol=2, byrow=T)[,1]
tis_pred = T_names[1]; ## the tissue to be predicted
#################################################
a0 = 1; b0 = 3; alp0 = 1; bet0 = 1; burnin = 100; total = 200; ## mcmc hyperparameters
s_ths = 70;
fold = 5;
chr = 1; glist = dir(paste0("/net/zhao/yh367/GTEX/cis_snp_by_gene/chr", chr));
k = 1; g = glist[k];
dose_path = paste0("/net/zhao/yh367/GTEX/cis_snp_by_gene/chr", chr, "/", g, "/", g, ".mach.dose")
Yt = dir(paste0("/net/zhao/yh367/GTEX/expr_gtex_portal/chr", chr, "/", g, "/"))

if(tis_pred%in%Yt){
	dir.create(paste0("/net/zhao/yh367/GTEX/Results/", tis_pred), showWarnings = FALSE)
	dir.create(paste0("/net/zhao/yh367/GTEX/Results/", tis_pred, "/chr", chr), showWarnings = FALSE)
	dir.create(paste0("/net/zhao/yh367/GTEX/Results/", tis_pred, "/chr", chr, "/", g), showWarnings = FALSE)
	setwd(paste0("/net/zhao/yh367/GTEX/Results/", tis_pred, "/chr", chr, "/", g))
	## expr files ##
	Y = list()
	for(t in 1:length(Yt)){
		Y[[t]] = read.table(paste0("/net/zhao/yh367/GTEX/expr_gtex_portal/chr", chr, "/", g, "/", Yt[t]), header=F)
	}
	ssize = unlist(lapply(Y, nrow))
	Yt = Yt[ssize >= s_ths]
	Y = Y[ssize >= s_ths]
	ssize = ssize[ssize >= s_ths]
	T_num = length(Yt)

	## genotype files ##
	dose = read.table(dose_path, header=F)
	sub_id = matrix(unlist(strsplit(dose[,1], "->")), ncol=2, byrow=T)[,1]
	M = ncol(dose) - 2
	pred_index = which(Yt==tis_pred)
	sub_id_map = list()
	for(t in 1:T_num){
		tmp = rep(0, nrow(Y[[t]]))
		for(j in 1:length(tmp)){
			tmp[j] = which(sub_id==Y[[t]][j,1])
		}
		sub_id_map[[t]] = tmp
	}

	pred_ids = sub_id_map[[pred_index]]

	valid_num = floor(length(pred_ids)/fold) # number of individuals in validation 
	s0_0 = rigamma(1, alpha = alp0, beta = bet0); se_0 = rigamma(T_num, alpha = alp0, beta = bet0); sg_0 = rigamma(T_num, alpha = alp0, beta = bet0); p_0 = rbeta(1, a0, b0); 
	Z_0 = matrix(rb(M*T_num, p_0), M, T_num); b_0 = rnorm(M, mean = 0, sd = sqrt(s0_0)); B_0 = matrix(0, M, T_num)
	for(t in 1:T_num){
		B_0[,t] = rnorm(M, mean=0, sd = sqrt(sg_0[t])) + Z_0[,t]*b_0
	}
	initial0 = list(s0_0 = s0_0, se_0 = se_0[pred_index], sg_0 = sg_0[pred_index], p_0 = p_0, Z_0 = as.matrix(Z_0[,pred_index]), b_0 = b_0, B_0 = as.matrix(B_0[,pred_index]))

	for(f in 1:fold){
		test_index = sample(1:length(pred_ids), size = valid_num)
		test_id = pred_ids[test_index]
		writeLines(sub_id[test_id], paste0("ids_cv", f))
		X_test = apply(as.matrix(dose[sub_id_map[[pred_index]][test_index],-c(1,2)]),2,as.numeric)
		Y_test = Y[[pred_index]][test_index, 2]
		X_train = list()
		Y_train = list()
		for(t in 1:T_num){
			X_include = sub_id_map[[t]][!(sub_id_map[[t]]%in%test_id)]
			Y_include = !(sub_id_map[[t]]%in%test_id)
			X_train[[t]] = apply(as.matrix(dose[X_include,-c(1,2)]),2,as.numeric)
			Y_train[[t]] = Y[[t]][Y_include, 2]
		}
		cat("Training sample sizes: ", unlist(lapply(Y_train, length)), "\n")
		cv_res = gibbs(initial = initial0, burnin = burnin, total = total, X = list(X_train[[pred_index]]), Y = list(Y_train[[pred_index]]), a0=a0, b0=b0, alp0=alp0, bet0=bet0)
		save(cv_res, X_test, Y_test, file = paste0("cv_res_single", f, ".RData"))
		## calculate Rsq ##
		beta_est = cv_res$B_avg#[, pred_index]
		Rsq = 1 - sum((Y_test - X_test%*%beta_est)^2)/sum((Y_test - mean(Y_test))^2)
		cat("R square = ", Rsq, "\n")
	}
}
#unlist(lapply(sub_id_map, length))
#unlist(lapply(sub_id_map, function(x){length(unique(x))}))
#unlist(lapply(Y, function(x){length(unique(x[,1]))}))
#unlist(lapply(Y, function(x){nrow(x)}))
else{
	print(paste0("Gene ", g, " not expressed in ", tis_pred))
}
### cross-validation ###
## training ##

