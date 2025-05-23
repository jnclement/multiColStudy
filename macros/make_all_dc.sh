#!/bin/bash

bit=$1
shift
TAG=$1
shift
RNS=($@)
declare -a NFILE
for (( i=0; i<$#; i++ )); do
    NFILE[i]=`wc -l < /sphenix/user/jocl/projects/multiColStudy/run/lists/dst_calofitting_run2pp-000${RNS[i]}.list`
done

root -l -b -q 'make_hists.C("'$TAG'",{'`(IFS=','; echo "${RNS[*]}")`'},{'`(IFS=','; echo "${NFILE[*]}")`'},'$bit')'
if [ $bit -eq 18 ]; then
    root -l -b -q 'make_tturn.C("'$TAG'",{'`(IFS=','; echo "${RNS[*]}")`'},{'`(IFS=','; echo "${NFILE[*]}")`'})'
fi
