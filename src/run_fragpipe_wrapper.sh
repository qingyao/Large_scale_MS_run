#!/bin/bash

year_month_folder=$1

while IFS=$'\t' read -r col1 col2; do
    echo $col1 $col2
    ./run_fragpipe.sh $col1 $col2 $year_month_folder || echo "Error processing: $col1 $col2" >&2
done < ../$year_month_folder/pxdid_taxid_to_run_fragpipe.txt

python3 write_done_file.py $year_month_folder