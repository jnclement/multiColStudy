#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh new
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
echo "got into condor scratch dir"
RN=$1
TAG=20250702
mkdir -p multicoltree
mkdir -p multicolhist
echo "copy file from tg to here"
cp /sphenix/tg/tg01/jets/jocl/multiCol/$RN/*$TAG* multicoltree
mkdir -p /sphenix/user/jocl/projects/multiColStudy/output/hists
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/make_all_dc.sh .
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/make_hists.C .
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/make_tturn.C .
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/make_njl_only.C .
echo "got all files, run code"
bash make_all_dc.sh 22 $TAG $RN
#bash make_all_dc.sh 18 $TAG $RN
#bash make_all_dc.sh 26 $TAG $RN
echo "copy file back"
cp multicolhist/* /sphenix/user/jocl/projects/multiColStudy/output/hists
