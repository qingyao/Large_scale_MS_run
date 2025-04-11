#!/bin/bash

PXDID=$1
TAXID=$2
YEAR_MONTH=$3
TAXID_SHORT=$(echo $TAXID|cut -d'-' -f1)
workdir=../$YEAR_MONTH/fragpipe_processing_output/$PXDID

WORKFLOW_FILE=${workdir}/${PXDID}-${TAXID_SHORT}.workflow;
rm -r $workdir/${TAXID_SHORT}; ## if previously has processed and needs reprocessing.
mkdir -p $workdir/${TAXID_SHORT};
../fragpipe/bin/fragpipe --headless \
--workflow $WORKFLOW_FILE \
--ram 128 \
--manifest $workdir/${PXDID}-${TAXID}.manifest \
--workdir $workdir/${TAXID_SHORT} \
2>&1> $workdir/${TAXID_SHORT}/out.log & pid=$(echo $!);
python3 ./monitor_pid.py $pid fragpipe $workdir/${TAXID_SHORT}/ram_cpu_usage_log.txt ## this pid doesn't account for all threads, thus use keywords