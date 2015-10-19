#! /usr/bin/env python
import os
import getopt
import sys

path        = "/tmp/higgs313/"
prefix      = "/tmp/higgs313/"
destination = "/tmp/higgs313/Herwig_0T_2015A_ReConst"

DIR = os.listdir(path)

ss = "Herwig2015AMagnetOff"

for dd in DIR:
  print "Reading directory "+dd
  ROOTFILES = os.listdir(path+"/"+dd)
  command = "hadd -f "+destination+"/"+dd+"_"+ss+".root "
  for ll in ROOTFILES:
    command += prefix+dd+"/"+ll+" "
  print command
  os.system(command)
