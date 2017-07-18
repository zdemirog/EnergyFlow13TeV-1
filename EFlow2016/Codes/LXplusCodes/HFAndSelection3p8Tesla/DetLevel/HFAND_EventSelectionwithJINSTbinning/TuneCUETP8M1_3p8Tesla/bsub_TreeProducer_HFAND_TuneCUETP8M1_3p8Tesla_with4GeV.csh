#!/bin/csh
##############################################################
# In order to send into the batch system the following files 
# bsub -o Output.txt -q 8nh bsub_script.csh
# to check the status of files use the command: bjobs
#############################################################

cd /afs/cern.ch/work/z/zdemirog/CFF_04June/CMSSW_7_5_0_pre1_ROOT5/src/EnergyFlow13TeV/cpp_macro/HFAND_EventSelectionwithJINSTbinning/TuneCUETP8M1_3p8Tesla
cmsenv
root -l EnergyFlow_DetLevel_TreeProducer_HFAND_TuneCUETP8M1_3p8Tesla_with4GeV.C++

