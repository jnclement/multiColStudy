#!/bin/bash

echo $evt
for rn in {0..1}; do
    echo $rn
    bash run_everything.sh $rn
done


