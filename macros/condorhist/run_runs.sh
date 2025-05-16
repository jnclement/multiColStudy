#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Need region argument"
    exit 1
fi

REG=$1

echo $evt
for rn in `cat ../../run/lists/rns${REG}.list`; do
    echo $rn
    bash run_everything.sh $rn
done


