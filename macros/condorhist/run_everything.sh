#!/bin/bash

RN=$1
BASENAME="condor_histmaker_${RN}"
SUBNAME="${BASENAME}.sub"
PREFIX=.
echo "executable = earlysim.sh" > $PREFIX/$SUBNAME
echo "concurrency_limits=CONCURRENCY_LIMIT_DEFAULT:250" >> $PREFIX/$SUBNAME
echo "arguments = ${RN}" >> $PREFIX/$SUBNAME
echo "priority = 1000" >> $SUBNAME
echo "output = /sphenix/user/jocl/projects/multiColStudy/output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $PREFIX/$SUBNAME
echo "request_memory = 3000MB" >> $PREFIX/$SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $PREFIX/$SUBNAME
echo "error = /sphenix/user/jocl/projects/multiColStudy/output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $PREFIX/$SUBNAME
echo "queue 1" >> $PREFIX/$SUBNAME

condor_submit $PREFIX/$SUBNAME
