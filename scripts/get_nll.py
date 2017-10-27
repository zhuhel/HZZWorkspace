#!/usr/bin/env python
import math

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c") 

def get_nll(file_name, poi_name):
    f1 = ROOT.TFile.Open(file_name, "read")
    ws = f1.Get("combined")
    mc = ws.obj("ModelConfig")
    obs_data = ws.data("obsData")
    poi = ws.var(poi_name)
    nuis = mc.GetNuisanceParameters()
    global_obs = mc.GetGlobalObservables()
    ws.saveSnapshot("nominalGlobs", global_obs)
    ws.saveSnapshot("nominalNuis", nuis)
    nll = mc.GetPdf().createNLL(
        obs_data, 
        ROOT.RooFit.Constrain(nuis),
        ROOT.RooFit.GlobalObservables(global_obs)
    )
    np = ws.var("alpha_syst1")
    n_x = 100
    n_y = 1000
    h1 = ROOT.TH2F("h1","h1", n_x, 1., 2., n_y, 1E-4, 1E-1)
    for i in range(n_x):
        for j in range(n_y):
            x_val = h1.GetXaxis().GetBinCenter(i+1)
            y_val = h1.GetYaxis().GetBinCenter(j+1)
            poi.setVal(x_val)
            np.setVal(y_val)
            h1.SetBinContent(i+1, j+1, math.exp(-1*nll.getVal()))

    fout = ROOT.TFile.Open(file_name.replace("test","test_nll"), "recreate")
    h1.Write()
    fout.Close()
    f1.Close()

def plot_nll(file_name):
    ROOT.SetAtlasOpt()
    f1 = ROOT.TFile.Open(file_name, "read")
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    ROOT.SetAtlasStyleCanvas(canvas, True)
    h1 = f1.Get("h1")
    h1.SetXTitle("#mu")
    h1.SetYTitle("#theta")
    h1.SetZTitle("Likelihood")
    h1.SetNdivisions(8, "Z")
    ROOT.SetAtlasStyleHist(h1)
    h1.Draw("colz")
    canvas.SaveAs(file_name.replace("root","pdf"))
    f1.Close()

def get():
    get_nll("test_100.root", "SigXsecOverSM")
    get_nll("test_200.root", "SigXsecOverSM")
    get_nll("test_300.root", "SigXsecOverSM")

def plot():
    plot_nll("test_nll_5.root")
    plot_nll("test_nll_300.root")
if __name__ == "__main__":
    plot()
