#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.468
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
RN=$1
TAG=20250513
mkdir -p multicoltree
mkdir -p multicolhist
cp /sphenix/tg/tg01/jets/jocl/multiCol/$RN/*$TAG* multicoltree
mkdir -p /sphenix/user/jocl/projects/multiColStudy/output/hists
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/make_all_dc.sh .
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/make_hists.C .
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/make_tturn.C .

bash make_all_dc.sh 18 $RN
bash make_all_dc.sh 10 $RN

cp multicolhist/* /sphenix/user/jocl/projects/multiColStudy/output/hists
