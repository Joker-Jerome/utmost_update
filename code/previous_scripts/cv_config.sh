
output_path=/gpfs/loomis/scratch60/fas/zhao/yh367/BSLMM
cmds_path=/gpfs/loomis/scratch60/fas/zhao/yh367//gemma_cmds_by_gene
jobfile_path=/gpfs/loomis/scratch60/fas/zhao/yh367/cv_jobs
bslmm_script=/gpfs/loomis/scratch60/fas/zhao/yh367/cv_config.R
data=/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX

partition=general
module load R
Rscript --vanilla $bslmm_script $output_path $cmds_path $data $jobfile_path $partition

for i in {10001..12000}
do
	cd $jobfile_path/$i
	sbatch job${i}.sh
done


netid=yh367
output_path=/gpfs/loomis/scratch60/fas/zhao/${netid}/BSLMM
cmds_path=/gpfs/loomis/scratch60/fas/zhao/${netid}/gemma_cmds_by_gene
jobfile_path=/gpfs/loomis/scratch60/fas/zhao/${netid}/cv_jobs_shit.sh
bslmm_script=/gpfs/loomis/scratch60/fas/zhao/${netid}/cv_config_onesh.R
data=/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX
idx1=17619  #12001  #10001
idx2=17886  #15000  #12000

module load R
Rscript --vanilla $bslmm_script $output_path $cmds_path $data $jobfile_path $idx1 $idx2
module load SimpleQueue ## for farnam
cd /gpfs/loomis/scratch60/fas/zhao/${netid}
sqCreateScript -q general -N shit -n 60 -m 10000 -w 24:00:00 cv_jobs_shit.sh > cv_jobs_shit.sh.pbs
sbatch cv_jobs_shit.sh.pbs

find . -type f -name '*log*' -exec rm -f {} \;
find . -type f -name '*log*' -print -delete
find . -type f -name '*keep*' -print -delete
find . -type f -name '*nosex*' -print -delete


chmod -R 777 yh367









options(stringsAsFactors=F)
finished = dir("/SAY/archive/hz27-CC0937-Biostatistics-A/yh367/revision/gemma_cmds_by_gene/")
ref = read.table("/SAY/archive/hz27-CC0937-Biostatistics-A/yh367/revision/gene_index.ref", header=F)
unfinished = ref[!ref[,2]%in%finished,1]
out = paste0('cd /SAY/archive/hz27-CC0937-Biostatistics-A/yh367/revision/cv_jobs/', unfinished, '; sbatch job', unfinished, '.sh;')
writeLines(out, '/SAY/archive/hz27-CC0937-Biostatistics-A/yh367/revision/unfinished_cv_config.sh')


./gemma.linux -no-check -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile /gpfs/loomis/scratch60/fas/zhao/yh367/BSLMM/chr1/SCYL3/cv_folder//cv_1_tissue_1 -bslmm 1 -o cv_1_tissue_1 -outdir /gpfs/loomis/scratch60/fas/zhao/yh367/BSLMM/chr1/SCYL3/cv_folder/