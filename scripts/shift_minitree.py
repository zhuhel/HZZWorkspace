#!/usr/bin/env python
import ROOT
from array import array
ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c") 

def shift_m4l(in_name, out_name):
    tree_name = "tree_incl_all"
    infile = ROOT.loader(in_name, tree_name)
    nentries = infile.GetEntries()

    outfile = ROOT.TFile.Open(out_name, "recreate")
    outtree = ROOT.TTree(tree_name, tree_name)
    m4l = array('f', [0.])
    m4l_cst = array('f', [0.])
    weight = array('f', [0.])
    event_type = array('f', [0.])
    prod_type = array('f', [0.])
    met = array('f', [0.])
    w_lumi = array('f', [0.])
    run = array('i', [0])
    event = array('i', [0])
    outtree.Branch("m4l_constrained", m4l_cst, "m4l_constrained/F")
    outtree.Branch("m4l_unconstrained", m4l, "m4l_unconstrained/F")
    outtree.Branch("weight", weight, "weight/F")
    outtree.Branch("w_lumi", w_lumi, "w_lumi/F")
    outtree.Branch("event_type", event_type, "event_type/F")
    outtree.Branch("prod_type", prod_type, "prod_type/F")
    outtree.Branch("met_et", met, "met_et/F")
    outtree.Branch("run", run, "run/I")
    outtree.Branch("event", event, "event/I")
    for ientry in range(nentries):
        if infile.LoadTree(ientry) < 0:
            print ientry," has problem"
            break
        infile.GetEntry(ientry)
        m4l[0] = infile.m4l_unconstrained + 0.09
        m4l_cst[0] = infile.m4l_constrained + 0.09
        weight[0] = infile.weight
        event_type[0] = infile.event_type
        prod_type[0] = infile.prod_type
        met[0] = infile.met_et
        w_lumi[0] = infile.w_lumi
        run[0] = infile.run
        event[0] = infile.event
        outtree.Fill()

    outtree.Write()
    outfile.Close()

if __name__ == "__main__":
    #shift_m4l("signal_mH125.root", "signal_mH125_09.root")
    #shift_m4l("signal_mH125_notau.root", "signal_mH125_09_notau.root")
    shift_m4l("signal.list", "signal_mH125_09_notau.root")
