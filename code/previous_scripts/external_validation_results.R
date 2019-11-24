options(stringsAsFactors=F)
library(scales)
library(gplots)
paths = paste0("/ysm-gpfs/pi/zhao/from_louise/jw2372/results_geuv/chr", 1:22, "/regress.txt")

res = c()
for(i in 1:22){
	tmp = read.table(paths[i], header=F, skip=1)
	res = rbind(res, tmp)
}
res = res[complete.cases(res), ]

which.max(res[,2])
which.max(res[,3])
## qq_plot ##
sigP = sort(res[,2])
mulP = sort(res[,3])
N = 379
expP = rep(0, length(sigP))
set.seed(100)
for(i in 1:length(sigP)){
	tmp1 = rnorm(N)
	tmp2 = rnorm(N)
	expP[i] = cor(tmp1,tmp2)^2
}
expP = sort(expP)
difP = ks.test(sigP, mulP, alternative = c("greater"))
saveRDS(res, '/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/revision/simulation/causal_tissue/geuv.rds')
png("/ysm-gpfs/pi/zhao/from_louise/yh367/GEUV/external_results/geuv_original_expr.png", width = 800, height = 800)
par(mar=c(5,5,4,2))
plot(expP, sigP, pch=16, col=alpha("#E64A45",0.8),
     xlab="Expected R2", ylab="Observed predictive R2", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
points(expP, mulP, pch=16, col=alpha("#3665AD",0.8), cex=1.5)
abline(0,1)
text(0.005,0.6,"P=3.21e-5", cex=1.5)
legend("topleft", legend=c("Single-tissue elastic net","Multi-tissue joint prediction"), col = c("#E64A45","#3665AD"), cex=1.5, pch=c(16,16))
dev.off()


## CommonMind
options(stringsAsFactors=F)
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/with_pred1.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/with_pred2.RData")
sum(unlist(lapply(with_pred1,length)))
sum(unlist(lapply(with_pred2,length)))
ovp_genes = list()
ovp_genes_df = data.frame()
ovp_genes_num = rep(0,22)
for(chrom in 1:22){
	ovp_genes[[chrom]] = intersect(intersect(with_pred1[[chrom]][-1],with_pred2[[chrom]][-1]),dir(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/expr_by_gene/chr", chrom)))
	ovp_genes_num[chrom] = length(ovp_genes[[chrom]])
	ovp_genes_df = rbind(ovp_genes_df,cbind(chrom, ovp_genes[[chrom]]))
}
write.table(ovp_genes_df, "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/merged_gene_list.txt", quote=F, row.names=F, col.names=c("chrom","gene"))
save(ovp_genes, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/useful_gene_list.RData")

results_path = "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi_single_merged/"
dir.create(results_path, showWarnings=F)
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/useful_gene_list.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_list.RData")
IDs_map = read.table("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/genotype_expr_ID.map", header=T)

Rsq = Rsq_raw = list()
for(i in 1:4){
	Rsq[[i]] = list()
	Rsq_raw[[i]] = list()
	for(chrom in 1:22){
		Rsq[[i]][[chrom]] = matrix(0,ovp_genes_num[chrom],2)
		Rsq_raw[[i]][[chrom]] = matrix(0,ovp_genes_num[chrom],2)
	}
}

k = 0
for(chrom in 1:22){
	glist = ovp_genes[[chrom]]
	for(i in 1:ovp_genes_num[chrom]){
		g = glist[i]
		dir.create(paste0(results_path, "chr", chrom, '/', g), showWarnings=F)
		multi_pred = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/multi/chrom",chrom,'/',g,'/',g,".pred"), header=T)
		single_pred = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/single/chrom",chrom,'/',g,'/',g,".pred"), header=T)
		names(multi_pred) = c("dose_ID1", "multi_pred")
		names(single_pred) = c("dose_ID1", "single_pred")
		act_expr_paths = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/expr_by_gene/chr", chrom, '/', g, "/expr", 1:4, ".txt")
		tmp0 = merge(IDs_map, single_pred, by="dose_ID1")
		tmp1 = merge(tmp0, multi_pred, by="dose_ID1")
		for(j in 1:4){
			act_expr = read.table(act_expr_paths[j], header=T)
			cov_expr = read.table(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/CommonMind/Expressions/DorsolateralPrefrontalCortex/NormalizedExpression/Gene/",j,".COV"), header=T)
			expr_with_cov = merge(act_expr, cov_expr, by="DLPFC_RNA_Sequencing_Sample_ID")
			forreg = expr_with_cov[,-c(1)]
			expr_with_cov["adjExpr"] = resid(lm(Expr ~ ., data=forreg))
			tmp2 = merge(tmp1, expr_with_cov, by="DLPFC_RNA_Sequencing_Sample_ID")
			Rsq[[j]][[chrom]][i,1] = cor(tmp2$single_pred,tmp2$adjExpr)^2
			Rsq[[j]][[chrom]][i,2] = cor(tmp2$multi_pred,tmp2$adjExpr)^2
			Rsq_raw[[j]][[chrom]][i,1] = cor(tmp2$single_pred,tmp2$Expr)^2
			Rsq_raw[[j]][[chrom]][i,2] = cor(tmp2$multi_pred,tmp2$Expr)^2
			write.table(tmp2, paste0(results_path, "chr", chrom, '/', g, "/expr", j, ".txt"))
		}
		k = k + 1
		print(k)
	}
}


save(Rsq, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_adj.RData")
save(Rsq_raw, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_raw.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_adj.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/useful_gene_list.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/sig_gene_list.RData")
t = 12
tissue_names[t]
g1 = g_sig1[[tissue_names[t]]]; g2 = g_sig2[[tissue_names[t]]];
ovp_genes[[chrom]]
Rsq_stack_adj  = list()
for(i in 1:4){
	Rsq_stack_adj[[i]] = data.frame()
	for(chrom in 1:22){
		tmp = data.frame(ovp_genes[[chrom]], Rsq[[i]][[chrom]])
		names(tmp) = c('gene', 'sig', 'multi')
		Rsq_stack_adj[[i]] = rbind(Rsq_stack_adj[[i]], tmp)
	}
}

for(i in 1:4){
	print(sum(is.na(Rsq_stack_adj[[i]][,2])))
	print(sum(is.na(Rsq_stack_adj[[i]][,3])))
	print(sum(Rsq_stack_adj[[i]][,1]%in%g1))
	print(sum(Rsq_stack_adj[[i]][,1]%in%g2))
	print('\n')
}

Rsq_stack_test = list()
for(i in 1:4){
	Rsq_stack_test[[i]] = list()
	Rsq_stack_test[[i]][[1]] = Rsq_stack_adj[[i]][!is.na(Rsq_stack_adj[[i]][,2])&Rsq_stack_adj[[i]][,1]%in%g1,c(1,2)]
	Rsq_stack_test[[i]][[2]] = Rsq_stack_adj[[i]][!is.na(Rsq_stack_adj[[i]][,3])&Rsq_stack_adj[[i]][,1]%in%intersect(g1,g2),c(1,3)]
}
for(i in 1:4){
	print(c(nrow(Rsq_stack_test[[i]][[1]]),nrow(Rsq_stack_test[[i]][[2]])))
	print(c(mean(Rsq_stack_test[[i]][[1]][,2]),mean(Rsq_stack_test[[i]][[2]][,2])))
	print(ks.test(Rsq_stack_test[[i]][[1]][,2],Rsq_stack_test[[i]][[2]][,2],alternative='greater')$p)
	print(wilcox.test(Rsq_stack_test[[i]][[1]][,2],Rsq_stack_test[[i]][[2]][,2],alternative='greater')$p.value)
}
save.image('/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/plot_0721.RData')
for(i in 1:4){
	sigP = sort(Rsq_stack_test[[i]][[1]][,2])
	mulP = sort(Rsq_stack_test[[i]][[2]][,2])
	N = 100
	expP = rep(0, length(sigP))
	set.seed(100)
	for(j in 1:length(sigP)){
		tmp1 = rnorm(N)
		tmp2 = rnorm(N)
		expP[j] = cor(tmp1,tmp2)^2
	}
	expP = sort(expP)
	difP = ks.test(sigP, mulP, alternative = c("greater"))
	print(difP$p.value)
	png(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/commonmind/commonmind_expr",i,".png"), width = 800, height = 800)
	par(mar=c(5,5,4,2))
	plot(expP, sigP, pch=16, col=alpha("#E64A45",0.8),
	     xlab="Expected R2", ylab="Observed predictive R2", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
	points(expP, mulP, pch=16, col=alpha("#3665AD",0.8), cex=1.5)
	abline(0,1)
	#text(0.2,0.6,"P=3.21e-5", cex=1.5)
	legend("topleft", legend=c("Single-tissue elastic net","Multi-tissue joint prediction"), col = c("#E64A45","#3665AD"), cex=1.5, pch=c(16,16))
	dev.off()
}






Rsq_stack_adj = Rsq_stack_raw = list()
for(i in 1:4){
	Rsq_stack_adj[[i]] = Rsq_stack_raw[[i]] = data.frame()
	for(chrom in 1:22){
		Rsq_stack_adj[[i]] = rbind(Rsq_stack_adj[[i]], Rsq[[i]][[chrom]])
		Rsq_stack_raw[[i]] = rbind(Rsq_stack_raw[[i]], Rsq_raw[[i]][[chrom]])
	}
}
save(Rsq_stack_adj, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack_adj.RData")
save(Rsq_stack_raw, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack_raw.RData")
Rsq_stack_adj_cmp = Rsq_stack_raw_cmp = list()
for(i in 1:4){
	Rsq_stack_adj_cmp[[i]] = Rsq_stack_raw_cmp[[i]] = data.frame()
	#not_na1 = (!is.na(Rsq_stack_adj[[i]][,1]))&(!is.na(Rsq_stack_adj[[i]][,2]))
	#not_na2 = (!is.na(Rsq_stack_raw[[i]][,1]))&(!is.na(Rsq_stack_raw[[i]][,2]))
	Rsq_stack_adj_cmp[[i]] = Rsq_stack_adj[[i]][complete.cases(Rsq_stack_adj[[i]]), ]
	Rsq_stack_raw_cmp[[i]] = Rsq_stack_raw[[i]][complete.cases(Rsq_stack_raw[[i]]), ]
	cat("adj.single.mean = ", mean(Rsq_stack_adj_cmp[[i]][,1]), "; adj.multi.mean = ", mean(Rsq_stack_adj_cmp[[i]][,2]), "\n")
	cat("raw.single.mean = ", mean(Rsq_stack_raw_cmp[[i]][,1]), "; raw.multi.mean = ", mean(Rsq_stack_raw_cmp[[i]][,2]), "\n")
}

save(Rsq_stack_adj_cmp, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack_adj_cmp.RData")
save(Rsq_stack_raw_cmp, file="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack_raw_cmp.RData")

load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack_adj_cmp.RData")
load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/WeightsMatrix/Rsq_stack_raw_cmp.RData")
for(i in 1:4){
	print(dim(Rsq_stack_adj_cmp[[i]]))
	print(c(mean(Rsq_stack_adj_cmp[[i]][,1]),mean(Rsq_stack_adj_cmp[[i]][,2])))
	print(ks.test(Rsq_stack_adj_cmp[[i]][,1],Rsq_stack_adj_cmp[[i]][,2],alternative='less'))
}
for(i in 1:4){
	sigP = sort(Rsq_stack_adj_cmp[[i]][,1])
	mulP = sort(Rsq_stack_adj_cmp[[i]][,2])
	N = 146
	expP = rep(0, length(sigP))
	set.seed(100)
	for(j in 1:length(sigP)){
		tmp1 = rnorm(N)
		tmp2 = rnorm(N)
		expP[j] = cor(tmp1,tmp2)^2
	}
	expP = sort(expP)
	difP = ks.test(sigP, mulP, alternative = c("greater"))
	print(difP$p.value)
	png(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GEUV/external_results/commonmind_adjusted_expr",i,".png"), width = 800, height = 800)
	par(mar=c(5,5,4,2))
	plot(expP, sigP, pch=16, col=alpha("#E64A45",0.8),
	     xlab="Expected R2", ylab="Observed predictive R2", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
	points(expP, mulP, pch=16, col=alpha("#3665AD",0.8), cex=1.5)
	abline(0,1)
	#text(0.2,0.6,"P=3.21e-5", cex=1.5)
	legend("topleft", legend=c("Single-tissue elastic net","Multi-tissue joint prediction"), col = c("#E64A45","#3665AD"), cex=1.5, pch=c(16,16))
	dev.off()
}


for(i in 1:4){
	sigP = sort(Rsq_stack_raw_cmp[[i]][,1])
	mulP = sort(Rsq_stack_raw_cmp[[i]][,2])
	N = 146
	expP = rep(0, length(sigP))
	set.seed(100)
	for(j in 1:length(sigP)){
		tmp1 = rnorm(N)
		tmp2 = rnorm(N)
		expP[j] = cor(tmp1,tmp2)^2
	}
	expP = sort(expP)
	difP = ks.test(sigP, mulP, alternative = c("greater"))
	print(difP$p.value)
	png(paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GEUV/external_results/commonmind_original_expr",i,".png"), width = 800, height = 800)
	par(mar=c(5,5,4,2))
	plot(expP, sigP, pch=16, col=alpha("#E64A45",0.8),
	     xlab="Expected R2", ylab="Observed predictive R2", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
	points(expP, mulP, pch=16, col=alpha("#3665AD",0.8), cex=1.5)
	abline(0,1)
	#text(0.2,0.6,"P=3.21e-5", cex=1.5)
	legend("topleft", legend=c("Single-tissue elastic net","Multi-tissue joint prediction"), col = c("#E64A45","#3665AD"), cex=1.5, pch=c(16,16))
	dev.off()
}
