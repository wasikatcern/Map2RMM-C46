#!/bin/bash

echo "Setup PyROOT and numpy"
echo "Setup ROOT, PyROOT tensorflow"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "views LCG_105 x86_64-el9-gcc13-opt"
