#!/bin/bash

echo $evt
for rn in {0..99}; do
    echo $rn
    bash run_everything.sh $rn
done


