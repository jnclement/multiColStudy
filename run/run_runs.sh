#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Need tag argument (string), evtnum (int)"
    exit 1
fi

nmax=20
evt=$2
filecounter=0
if [ $evt -gt 10000 ]; then
    evt=0
fi
echo $evt
for rn in `ls  lists/dst_calo_run2pp*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}' | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l lists/dst_calo_run2pp-000${rn}.list | awk '{print $1}'`
    njob=$(( $nfile + 9 ))
    njob=$(( $njob / 10 ))
    filecounter=$(( $filecounter + $njob ))
    if [ $filecounter -gt $nmax ]; then
	break
    fi
#    nfile=$(( ($nfile + 9) / 10 ))
    mkdir -p /sphenix/tg/tg01/jets/jocl/multiCol/$rn
    echo $rn $filecounter
    bash run_everything.sh $1 $njob $evt
done


