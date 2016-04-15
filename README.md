# IHeartNY_13TeV
13 TeV implementation of code for boosted ttbar differential cross section
Currently designed to run on B2G EDM ntuple files and create slimmed trees for later analysis
Description of files:
  PSet.py -- defines input files for standalone running
  iheartny_topxs_fwlite.py -- analysis code to make trees
  execute_iheartNY.sh -- commands for standalone running
  crabConfig_Boosted.py -- crab3 config file
To run as standalone:
source execute_iheartNY.sh
To run using Crab3:
crab submit -c crabConfig_Boosted.py
