args = commandArgs(trailingOnly=TRUE)
task_index = args[1] ## a file contains task index
output_path = args[2] ## path for saving outputs of tasks
sqfile_path = args[3] ## path for saving simple queue files
## e.g. Rscript --vanilla cmds_gen.R xxx.txt /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1 /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/sq
## generates the commands for task index contained in xxx.txt and saves the output to /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1 and all the simple queue files will be saved in /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/sq

load("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/cmds_prefix.RData")
tsk_idx = as.numeric(readLines(task_index))
output_cmds = paste(cmds_prefix[start:end], output_path)
dir.create(output_path, showWarnings = FALSE)
dir.create(sqfile_path, showWarnings = FALSE)
setwd(sqfile_path)
for(i in 1:length(output_cmds)){
	dir.create(paste0(i), showWarnings = FALSE)
	writeLines(output_cmds[i], paste0(i, "/task", i, ".sh"))
}
cat("Simple Queue files written to ", sqfile_path, '\n')
cat("# of tasks: ", length(output_cmds), '\n')


