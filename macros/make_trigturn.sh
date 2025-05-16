#!/bin/bash

RNS=($@)
declare -a NFILE
for (( i=0; i<$#; i++ )); do
    NFILE[i]=`wc -l < ../run/lists/dst_calofitting_run2pp-000${RNS[i]}.list`
done
root -l -b -q 'make_tturn.C("20250513",{'`(IFS=','; echo "${RNS[*]}")`'},{'`(IFS=','; echo "${NFILE[*]}")`'})'
