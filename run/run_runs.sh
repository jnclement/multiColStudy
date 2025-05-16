#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Need tag argument (string), evtnum (int)"
    exit 1
fi

nmax=20000
evt=$2
filecounter=0
if [ $evt -gt 100000 ]; then
    evt=0
fi
echo $evt
for rn in `ls  lists/dst_calofitting_run2pp*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}' | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l lists/dst_calofitting_run2pp-000${rn}.list | awk '{print $1}'`
    njob=$(( $nfile + 4 ))
    njob=$(( $njob / 5 ))
    filecounter=$(( $filecounter + $njob ))
    if [ $filecounter -gt $nmax ]; then
	break
    fi
#    mkdir -p /sphenix/tg/tg01/jets/jocl/multiCol/$rn
#    mkdir -p /sphenix/tg/tg01/jets/jocl/multiCol/$rn
    echo $rn $filecounter
    bash run_everything.sh $1 $njob $evt 0 $rn
done


