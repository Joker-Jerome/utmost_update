#!/bin/bash
OUTPUT_FILE=$1
#START_IDX=$2
#END_IDX=$3


if [ -f "$OUTPUT_FILE" ]; then
	echo "$OUTPUT_FILE exists"
	rm $OUTPUT_FILE
fi

for i in {1..22}
	do echo "Rscript --vanilla /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/R/gene_info.R ${i}" >> $OUTPUT_FILE
done

