#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Need tag argument (string), evtnum (int), nlist (int)"
    exit 1
fi

nmax=20
evt=$2
filecounter=0
if [ $evt -gt 100000 ]; then
    evt=0
fi
echo $evt
for rn in `cat lists/rns${3}.list`; do #`ls  lists/dst_calofitting_run2pp*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}' | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    #if [ $rn -gt 47514 ]; then
#	break
 #   fi
  #  if [ $rn -ne 47514 ]; then
#	continue
 #   fi
    nfile=`wc -l lists/dst_calofitting_run2pp-000${rn}.list | awk '{print $1}'`
    njob=$(( $nfile + 3 ))
    njob=$(( $njob / 4 ))
    filecounter=$(( $filecounter + $njob ))
    if [ $filecounter -gt $nmax ]; then
	break
    fi
#    if [ $rn -gt 49349 ]; then
#	continue
#    fi
#    if [ $rn -gt 49111 ]; then
#	continue
#    fi
#    mkdir -p /sphenix/tg/tg01/jets/jocl/multiCol/$rn
#    mkdir -p /sphenix/tg/tg01/jets/jocl/multiCol/$rn
    echo $rn $filecounter
    bash run_everything.sh $1 $njob $evt 0 $rn
done


