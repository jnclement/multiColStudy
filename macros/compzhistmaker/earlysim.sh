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
TAG=20250612sim
#mkdir -p multicoltree
mkdir -p ./multicolhist
echo "copy file from tg to here"
#cp /sphenix/tg/tg01/jets/jocl/multiCol/$RN/*$TAG* multicoltree
mkdir -p /sphenix/user/jocl/projects/multiColStudy/output/hists
cp -r /sphenix/user/jocl/projects/multiColStudy/macros/comp_zvtx.C .
echo "got all files, run code"
root -b -q -l 'comp_zvtx.C("'${TAG}'",'$RN')'
echo "copy file back"
cp multicolhist/* /sphenix/user/jocl/projects/multiColStudy/output/hists
