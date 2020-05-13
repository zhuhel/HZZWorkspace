import math
import os

import ROOT
import AtlasStyle
ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c") 

import monoH
import common
#import PyCintex

from ROOT import *


tree = "tree_incl_all"
def get_met_change(fileName1, fileName2, fileName3, outname):
    print 'File1: ', fileName1
    print 'File2: ', fileName2
    print 'File3: ', fileName3
    # open the file
    file1 = TFile.Open(fileName1)
    file2 = TFile.Open(fileName2)
    file3 = TFile.Open(fileName3)
    tree_f1 = file1.Get(tree)
    tree_f2 = file2.Get(tree)
    tree_f3 = file3.Get(tree)

    tree_f2.BuildIndex("run", "event")
    tree_f3.BuildIndex("run", "event")
    n_f1 = tree_f1.GetEntries()
    n_f2 = tree_f2.GetEntries()
    if n_f1 != n_f2:
        print "different number of events"
    met_nom = TH1F("met_nom", "met_nom", 30, 0, 300)
    met_up = TH1F("met_up", "met_up", 30, 0, 300)
    met_down = TH1F("met_down", "met_down", 30, 0, 300)

    for iEvt in range(n_f1):
        tree_f1.GetEntry(iEvt)
        if (tree_f1.met_et < 30):
            met_nom.Fill(tree_f1.met_et, tree_f1.weight)
            found1 = tree_f2.GetEntryWithIndex(tree_f1.run, tree_f1.event)
            if found1 < 0:
                print "Cannot find",tree_f1.event," in file2"
                continue
            else:
                met_up.Fill(tree_f2.met_et, tree_f2.weight)

            found2 = tree_f3.GetEntryWithIndex(tree_f1.run, tree_f1.event)
            if found2 < 0:
                print "Cannot find",tree_f1.event," in file3"
                continue
            else:
                met_down.Fill(tree_f3.met_et, tree_f3.weight)

    fileout = TFile.Open(outname, "recreate")
    met_nom.Write()
    met_up.Write()
    met_down.Write()
    fileout.Close()

    file1.Close()
    file2.Close()
    file3.Close()

def plot_met_cmp(file_name, npname):
    ROOT.gROOT.SetBatch(True)
    fin = TFile.Open(file_name, "read")
    met_nom = fin.Get("met_nom")
    met_up = fin.Get("met_up")
    met_down = fin.Get("met_down")
    met_nom.GetYaxis().SetRangeUser(1e-4, 10)
    met_nom.SetXTitle("E_{T}^{miss} [GeV]")
    canvas = TCanvas("canvas", "canvas", 600, 600)
    canvas.SetLogy()
    met_nom.Draw()
    met_up.SetLineColor(2)
    met_down.SetLineColor(4)
    met_up.Draw('same')
    met_down.Draw('same')
    legend = myLegend(0.4, 0.7, 0.85, 0.9)
    legend.AddEntry(met_nom, "events E_{T}^{miss} < 30 GeV", "l")
    legend.AddEntry(met_up, "same events: "+npname+" up", "l")
    legend.AddEntry(met_down, "same events: "+npname+" down", "l")
    legend.Draw()
    canvas.SaveAs(file_name.replace("root","pdf"))

    ibin = met_nom.FindBin(80)
    total_bin = met_nom.GetNbinsX()
    nabove1 = met_up.Integral(ibin, total_bin)
    nabove2 = met_down.Integral(ibin, total_bin)

    print nabove1/met_up.Integral(), nabove2/met_down.Integral()
    fin.Close()

def compare_met(fileName1, fileName2, fileName3, npname, sample_name):
    ROOT.gROOT.SetBatch(True)
    ch1 = loader(fileName1, tree)
    ch2 = loader(fileName2, tree)
    ch3 = loader(fileName3, tree)
    max_x = 300
    met_f1 = TH1F("met_f1", "f1", int(max_x/10), 0, max_x)
    met_f2 = TH1F("met_f2", "f2", int(max_x/10), 0, max_x)
    met_f3 = TH1F("met_f3", "f3", int(max_x/10), 0, max_x)
    canvas = TCanvas("canvas", "canvas", 600, 600)
    canvas.SetLogy()
    weight_name=""
    ch1.Draw("met_et >> met_f1", weight_name)
    ch2.Draw("met_et >> met_f2", weight_name)
    ch3.Draw("met_et >> met_f3", weight_name)
    #met_f1.SetLineColor(2)
    met_f2.SetLineColor(2)
    met_f3.SetLineColor(4)
    met_f1.SetMarkerSize(0.5)
    met_f2.SetMarkerSize(0.5)
    met_f3.SetMarkerSize(0.5)
    met_f2.SetMarkerColor(2)
    met_f3.SetMarkerColor(4)
    met_f1.SetXTitle("E_{T}^{miss} [GeV]")
    h_bkg_list = ROOT.TList()
    h_bkg_list.Add(met_f2)
    h_bkg_list.Add(met_f3)
    h_bkg_list.SetOwner()
    
    pad1 = ROOT.add_ratio_pad(met_f1, h_bkg_list)
    pad1.cd()
    pad1.SetLogy()
    met_f1.GetYaxis().SetRangeUser(3, 10e5)
    #met_f1.GetYaxis().SetRangeUser(3e-4, 10e3)
    met_f1.Draw("EP")
    met_f2.Draw("same EP")
    met_f3.Draw("same EP")
    legend = myLegend(0.4, 0.7, 0.85, 0.9)
    legend.AddEntry(met_f1, "Nominal", "l")
    legend.AddEntry(met_f2, npname+" Up", "l")
    legend.AddEntry(met_f3, npname+" Down", "l")
    legend.Draw()
    canvas.SaveAs("pdf/"+sample_name+"_met_"+npname+".pdf")
    h_bkg_list.Add(met_f1)
    h_bkg_list.Clear()

def compare_met_sys(file_name, sample_name):
    #shape_list = common.shape_sys_list
    shape_list = [
        ("MET_SoftTrk_Scale","ATLAS_MET_SoftTrk_Scale"),
        ("MET_SoftTrk_Reso","ATLAS_MET_SoftTrk_Reso"),
    ]
    for shape_sys,shape_name in shape_list:
        file_down, file_up = common.get_shape_downup(file_name, shape_sys)
        if not os.path.isfile(file_up): continue
        if not os.path.isfile(file_down): file_down = file_up 
        # print file_name, file_down, file_up
        compare_met(file_name, file_down, file_up, shape_sys, sample_name)

def compare_monoH():
    #sample_name = "ZHlvlv"
    #compare_met_sys(monoH.samples_bkg[sample_name], sample_name)
    #for sample in monoH.bkg_samples:
        #compare_met_sys(monoH.samples_bkg[sample], sample)

    sample_name = "shxx"
    compare_met_sys(monoH.samples_sig["200"]["shxx"], sample_name+"_200")

    #sample_name = "zphxx"
    #compare_met_sys(monoH.samples_sig["200"]["zphxx"], sample_name+"_200")

if __name__ == '__main__':
    #main()
    compare_monoH()
