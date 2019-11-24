args = commandArgs(trailingOnly=TRUE)
output_path = args[1] ## path for saving outputs of tasks
cmds_path = args[2]
data = args[3]
jobfile_path = args[4] ## path for saving job files
partition = args[5]

#"tar -zcvf cv_folder.tar.gz cv_folder && rm -rf cv_folder"
data='/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX'
#output_path='/SAY/archive/hz27-CC0937-Biostatistics-A/yh367/revision/BSLMM'
cmds_path='/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/revision/BSLMM/gemma_cmds_by_gene'
#library(devtools)
#devtools::install_github("gabraham/plink2R/plink2R")
load("/gpfs/loomis/project/fas/zhao/yh367/utmost_revision/gene_info.RData")
tsk_idx = 1:nrow(useful)
cv_path = paste0(output_path, '/chr', useful[,1], '/', useful[,3], '/cv_folder')
output_cmds = paste0("module load R; Rscript --vanilla /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/revision/BSLMM/bslmm_generate_cmds.R ", useful[,1], ' ', useful[,2], ' ', useful[,3], ' ', output_path, ' ', data, ' ', cmds_path)

writeLines(output_cmds, '/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/revision/BSLMM/get_shit_cmds.sh')


#dvd = matrix(unlist(strsplit(cmds_prefix, ' ')), nrow=17886, byrow=T)
#useful = dvd[,7:9]
#save(useful, file='/gpfs/loomis/project/fas/zhao/yh367/utmost_revision/gene_info.RData')


cd /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/revision/BSLMM/
module load SimpleQueue ## for farnam
sqCreateScript -q general -N shit -n 80 -m 5000 -w 24:00:00 get_shit_cmds.sh > get_shit_cmds.sh.pbs
sbatch get_shit_cmds.sh.pbs

writeLines(paste0(1:length(gene_list), '\t', gene_list), '/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/revision/BSLMM/gene_index.ref')


15000   TBC1D3


/gpfs/loomis/scratch60/fas/zhao/yh367/


.fam, 