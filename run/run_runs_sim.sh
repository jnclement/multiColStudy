#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Need tag argument (string), evtnum (int)"
    exit 1
fi
evt=$2
nfile=$3
issim=1
rn=0
#for i in {0..25000}; do
#    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/$i
#    mkdir -p /sphenix/tg/tg01/jets/jocl/err/$i
#    mkdir -p /sphenix/tg/tg01/jets/jocl/out/$i
#    if [ $(($i % 100)) -eq 0 ]; then
#	echo $i
#    fi
#done
if [ $evt -gt 1000 ]; then
    evt=0
fi
bash run_everything.sh $1 $nfile $evt $issim $rn
