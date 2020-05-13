#!/usr/bin/env python
"""
print nuisance parameters in the workspace
"""

import ROOT
import sys
from optparse import OptionParser

usage = "%prog [options] file_name"
parser = OptionParser(usage=usage, description="Print nuisance parameters")
parser.add_option("-w", "--wsname", dest='ws_name', default='combined')
parser.add_option("-m", '--mcname', dest='mc_name', default='ModelConfig')

options,args = parser.parse_args()
#print options
#print args

if len(args) < 1:
    print parser.print_help()
    exit(1)


fin = ROOT.TFile.Open(args[0])
if not fin:
    print args[0],"does not exist"
    exit(1)

ws = fin.Get(options.ws_name)
mc = ws.obj(options.mc_name)
np_sets = mc.GetNuisanceParameters()
iter_np = ROOT.TIter(np_sets.createIterator())
obj = iter_np()
while obj:
    print obj.GetName()
    obj = iter_np()

fin.Close()
