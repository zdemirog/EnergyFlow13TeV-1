#!/bin/csh
##############################################################
# In order to send into the batch system the following files 
# bsub -o Output.txt -q 8nh bsub_script.csh
# to check the status of files use the command: bjobs
#############################################################

cd /afs/cern.ch/work/z/zdemirog/CFF_04June/CMSSW_7_5_0_pre1_ROOT5/src/EnergyFlow13TeV/cpp_macro/HFOR_EventSelectionwithJINSTbinning/Run255019_3p8Tesla
cmsenv
root -l EnergyFlow_DetLevel_HFOR_ZeroBias1_r255019_TreeProducer_with5GeV.C++

