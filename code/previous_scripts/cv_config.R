args = commandArgs(trailingOnly=TRUE)
output_path = args[1] ## path for saving outputs of tasks
cmds_path = args[2]
data = args[3]
jobfile_path = args[4] ## path for saving job files
partition = args[5]

#"tar -zcvf cv_folder.tar.gz cv_folder && rm -rf cv_folder"


#library(devtools)
#devtools::install_github("gabraham/plink2R/plink2R")
load("/gpfs/loomis/project/fas/zhao/yh367/utmost_revision/gene_info.RData")
tsk_idx = 1:nrow(useful)
cv_path = paste0(output_path, '/chr', useful[,1], '/', useful[,3], '/cv_folder')
output_cmds = paste0("module load R; Rscript --vanilla /gpfs/loomis/scratch60/fas/zhao/yh367/bslmm_generate_cmds.R ", useful[,1], ' ', useful[,2], ' ', useful[,3], ' ', output_path, ' ', data, ' ', cmds_path)
dir.create(output_path, showWarnings = FALSE)
dir.create(jobfile_path, showWarnings = FALSE)
setwd(jobfile_path)
for(i in tsk_idx){
	dir.create(paste0(i), showWarnings = FALSE)
	fname = paste0(i, "/job", i, ".sh")
	#file.remove(fname, showWarnings = FALSE)
	job_content = c()
	job_content = c(job_content, '#!/bin/bash')
	job_content = c(job_content, paste0('#SBATCH -p ', partition))
	job_content = c(job_content, paste0('#SBATCH -J ', i))
	job_content = c(job_content, '#SBATCH -t 24:00:00')
	job_content = c(job_content, '#SBATCH -N 1 --mem-per-cpu 10G')
	job_content = c(job_content, output_cmds[i])
	writeLines(job_content, fname)
}
cat("Job files written to ", jobfile_path, '\n', sep = "")
cat("# of tasks: ", length(tsk_idx), '\n', sep = "")



#dvd = matrix(unlist(strsplit(cmds_prefix, ' ')), nrow=17886, byrow=T)
#useful = dvd[,7:9]
#save(useful, file='/gpfs/loomis/project/fas/zhao/yh367/utmost_revision/gene_info.RData')