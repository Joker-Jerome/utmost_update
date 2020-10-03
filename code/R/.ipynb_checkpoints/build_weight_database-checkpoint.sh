# change dir 
weight_dir=/gpfs/scratch60/zhao/zy92/GTEX/weight_ENET/
cd $weight_dir

# scripts
bash /gpfs/project/zhao/zy92/GTEX/weight_normalized_pruned/merge_file.sh
cd merged
cp /gpfs/scratch60/zhao/zy92/GTEX/weight_ENET_0925/merged/build_database.sh ./
cp /gpfs/loomis/project/zhao/zy92/GTEX/weight_normalized_pruned/merged/make_sqlite_db.py ./
cp /gpfs/scratch60/zhao/zy92/GTEX/weight_ENET_0925/merged/build_database_task.txt ./

# dsq files 
dSQ --jobfile build_database_task.txt -J database --mem-per-cpu=16G -t 12:00:00 -p bigmem,general,pi_zhao --batch-file build_database_run.sh
sbatch build_database_run.sh

