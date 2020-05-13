#!/usr/bin/env python
import ROOT
import low_mass
import high_mass
import AtlasStyle
import common
import os
ROOT.gROOT.LoadMacro("/afs/cern.ch/user/x/xju/tool/loader.c") 

out_dir = "pdf"
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
config = high_mass
def compare(file_name1, file_name2, hist_name, cut, postfix):
    if not os.path.isfile(file_name1) or not os.path.isfile(file_name2):
        return
    ROOT.gROOT.SetBatch()
    chain = ROOT.loader(file_name1, "tree_incl_all")
    f2 = ROOT.TFile.Open(file_name2, "read")
    h1 = f2.Get(hist_name)
    if not h1:
        print hist_name, " does not exist"
        return
    h2 = h1.Clone("h2")
    chain.Draw("m4l_constrained>>h2", cut)
    c2 = ROOT.TCanvas("c2", "c2", 600, 600)
    h1.SetLineColor(2)
    h2.SetLineColor(4)
    if h2.Integral() == 0: 
        print cut,h1.GetName()
        return
    h2.Scale(h1.Integral()/h2.Integral())
    h2.Draw()
    h1.Draw("same")
    if h1.GetMaximumBin() > h1.GetNbinsX()/2.:
        legend = ROOT.myLegend(0.2, 0.8, 0.4, 0.9)
    else:
        legend = ROOT.myLegend(0.65, 0.8, 0.85, 0.9)
    legend.AddEntry(h1, "smoothed shape", "l")
    legend.AddEntry(h2, "raw shape", "l")
    legend.UseCurrentStyle()
    legend.Draw()
    c2.SaveAs(out_dir+"/"+postfix+"_"+hist_name+".pdf")

def compare_sample(sample, mass):
    if float(mass) > 1000: return
    global config
    is_bkg = True
    try:
        mass_value = float(mass)
        if mass_value > 10: is_bkg = False
    except ValueError:
        pass
    for sp_key, sp_value in sample.iteritems():
        if sp_value is None: continue
        if is_bkg:
            f2 = "test_"+sp_key+".root"
        else:
            f2 = "test_"+sp_key+"_"+mass+".root"
        print sp_value, f2
        for key, value in config.categories.iteritems():
            if True:
                compare(sp_value, f2, "m4l_"+key, "weight*"+value, sp_key+mass)
            else:
                compare(sp_value, f2, "m4l_"+key, "weight/w_xs/w_br*"+value, sp_key+mass)

def do_cmp():
    global config
    config = low_mass 
    sample = {}
    sample["all"] = "signal.root"
    compare_sample(sample, "125.09")

def cmp_low():
    global config
    config = low_mass 
    #config = high_mass
    #mass = "125.09"
    #compare_sample(config.samples_sig[mass], mass)
    for im in range(config.mass_points):
        mass = str(config.get_mass(im))
        #compare_sample(config.samples_sig[mass], mass)
    compare_sample(config.samples_bkg, "-1")

if __name__ == '__main__':
    #do_cmp()
    cmp_low()
