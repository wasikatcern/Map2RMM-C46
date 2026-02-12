#!/bin/bash


echo "Pythia8 setup"
# run it as: source setup.sh
echo "Set ROOT enviroment for Dijet+Lepton program"

HH=`hostname -A`
echo "HOST=$HH"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source $SCRIPT_DIR/lib/promc/setup.sh

echo "Setup ROOT, PyROOT tensorflow"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "views LCG_105 x86_64-el9-gcc13-opt"
echo "Setup Fastjet"
export FASTJET=/cvmfs/sft.cern.ch/lcg/latest/fastjet/3.4.1-5af57/x86_64-el9-gcc13-opt/
export PYTHIA8=/cvmfs/sft.cern.ch/lcg/latest/MCGenerators/pythia8/312-9ebb7/x86_64-el9-gcc13-opt
export PYTHIADIR=$PYTHIA8
export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
# LHAPDF6 configuration.
export LHAPDF6_USE=true
export LHAPDF6=/cvmfs/sft.cern.ch/lcg/latest/MCGenerators/lhapdf/6.5.4-f62b6/x86_64-el9-gcc13-opt
export LHAPDF6_BIN=$LHAPDF6/bin/
export PATH=$LHAPDF6_BIN:$PATH
export LHAPDF6_INCLUDE=$LHAPDF6/include/
export LHAPDF6_LIB=$LHAPDF6/lib
export LD_LIBRARY_PATH=$LHAPDF6_LIB/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH


