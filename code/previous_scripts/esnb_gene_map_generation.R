load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/protein_coding.RData")
xx = data.frame()
for(i in 1:22){
	xx = rbind(xx, pc_gene[[i]])
}
write.table(xx, "/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/esnb_gene.mapping", sep='\t', quote=F, row.names=F, col.names=T)