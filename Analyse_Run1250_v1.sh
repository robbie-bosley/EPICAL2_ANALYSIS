#!/bin/bash

echo "ILCSoft:"
source /cvmfs/clicdp.cern.ch/iLCSoft/builds/2020-02-07/x86_64-slc6-gcc7-opt/init_ilcsoft.sh

echo "====================================="
echo "starting job Analyse_mTower($1)"
echo "====================================="

root -l -b <<SHELL
.!pwd
.include /afs/cern.ch/user/r/rbosley/public/fastjet-install/include
gSystem->Load("/afs/cern.ch/user/r/rbosley/public/fastjet-install/lib/libfastjet.so");
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionToPush/classes/mTowerHit.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionToPush/classes/mTowerEvent.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionToPush/classes/mTowerClusterRobbie.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionToPush/classes/mTowerChipRobbie.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionToPush/processors/EventSelection_v1.cxx+
.L /afs/cern.ch/user/r/rbosley/public/EventSelectionToPush/Analyse_EPICAL2_v1.cxx+
Analyse_mTower(1250, 5.0)
.q
SHELL
