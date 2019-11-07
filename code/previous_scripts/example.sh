start=1
end=6000
output_path="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1"
sqfile_path="/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/glasso1SQ"

Rscript --vanilla /ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/codes/cmds_gen.R $start $end $output_path $sqfile_path

for i in {1..6000}
do
cd $sqfile_path/$i
sqCreateScript -q general -N task$i -n 1 -m 20000 -w 24:00:00 task$i.sh > task$i.pbs
sbatch task$i.pbs
done
