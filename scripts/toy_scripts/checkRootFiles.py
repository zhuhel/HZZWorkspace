#!/usr/bin/env python
from ROOT import gSystem
from ROOT import TFile
import commands
import os
import glob

badfiles=[]
for filename in glob.glob("*root"):
    #print filename
    f1 = TFile.Open(filename,"read")
    if  not f1 or f1.IsZombie():
        print "bad file: "+filename
        badfiles.append(filename)
        #os.system("rm "+filename)
    else:
        f1.Close()

print "total bad files: "+str(len(badfiles))
answer=str(raw_input("remove the bad files[y/n]?"))
if answer == 'y' or answer == 'Y':
    for f1 in badfiles:
        os.system("rm "+f1)

