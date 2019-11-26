#!/bin/bash
OUTPUT_FILE=$1
#START_IDX=$2
#END_IDX=$3


if [ -f "$OUTPUT_FILE" ]; then
	echo "$OUTPUT_FILE exists"
	rm $OUTPUT_FILE
fi

for i in {1..11}
	do echo "Rscript --vanilla /ysm-gpfs/project/zhao/zy92/utmost_update/code/R/adjust_expression_gtex8.R ${i}" >> $OUTPUT_FILE
done

