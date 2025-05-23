#!/bin/bash

TAG=$1
NFILE=$2
EVT=$3
if [ $# -lt 3 ]; then
    echo "Need arguments (in order): tag (string), nfile (int), EVT (int)"
    exit 1
fi

BASENAME="condor_${TAG}_${NFILE}_${EVT}_${4}_${5}"
PREFIX="." #"/sphenix/user/hanpuj/test"
SUBNAME="${BASENAME}.sub"

#echo "executable = containerscripts/earlydata.sh" > $PREFIX/$SUBNAME
echo "executable = earlysim.sh" > $PREFIX/$SUBNAME
#echo "concurrency_limits=CONCURRENCY_LIMIT_DEFAULT:250" >> $PREFIX/$SUBNAME
echo "arguments = ${TAG} \$(Process) ${EVT} ${4} ${5}" >> $PREFIX/$SUBNAME
echo "priority = 1000" >> $SUBNAME
echo "output = /sphenix/user/jocl/projects/multiColStudy/output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $PREFIX/$SUBNAME
echo "request_memory = 3000MB" >> $PREFIX/$SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $PREFIX/$SUBNAME
echo "error = /sphenix/user/jocl/projects/multiColStudy/output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $PREFIX/$SUBNAME
echo "queue ${NFILE}" >> $PREFIX/$SUBNAME

condor_submit $PREFIX/$SUBNAME
