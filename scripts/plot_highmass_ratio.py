#!/usr/bin/env python
import ROOT
import AtlasStyle
import math
import re
import os
from array import array
from ROOT import gStyle,gDirectory,TCanvas,TPad,TFile,TH1D,TF1,TLegend,TPaveText,TGaxis,TLine
if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("$HZZWSCODEDIR/scripts/loader.c")

def create_TGraphAsymmErrors(x_, nominal_, up_, down_):
    zero_ = array('f', [0]*len(x_))
    up_var_ = [x-y for x,y in zip(up_, nominal_)]
    down_var_ = [x-y for x,y in zip(nominal_, down_)]
    gr_error = ROOT.TGraphAsymmErrors(
        len(x_), array('f',x_), array('f',nominal_),
        zero_, zero_,
        array('f', down_var_), array('f', up_var_))
    return gr_error

def get_graph_from_list(mass_list, obs_list, exp_list, 
                        up_1sig_list, up_2sig_list,
                        down_1sig_list, down_2sig_list):

    zero_list = [0]*len(mass_list)
    print mass_list
    print obs_list
    gr_obs = ROOT.TGraph(len(mass_list), array('f', mass_list), array('f',obs_list))
    gr_exp = ROOT.TGraph(len(mass_list), array('f', mass_list), array('f', exp_list))
    gr_1sig = create_TGraphAsymmErrors(mass_list, exp_list, up_1sig_list, down_1sig_list)
    gr_2sig = create_TGraphAsymmErrors(mass_list, exp_list, up_2sig_list, down_2sig_list)
    gr_obs.SetLineWidth(2)
    gr_obs.SetLineStyle(1)
    gr_obs.SetMarkerStyle(20)

    gr_exp.SetLineWidth(2)
    gr_exp.SetLineStyle(2)

    gr_1sig.SetFillStyle(1001)
    gr_1sig.SetFillColor(ROOT.kGreen)
    gr_2sig.SetFillStyle(1001)
    gr_2sig.SetFillColor(ROOT.kYellow)
    
    return (gr_obs, gr_exp, gr_1sig, gr_2sig)


def make_limit_graph(file1):
    ROOT.gROOT.SetBatch()
    mass_list = []
    obs_list = []
    exp_list = []
    up_2sig_list = []
    up_1sig_list = []
    down_1sig_list = []
    down_2sig_list = []
    
    with open(file1, 'r') as f:
        for line in f:
            items = line.split()

            if (len(items)<9): continue
            if (math.isnan(float(items[3]))): continue
            if (math.isnan(float(items[4]))): continue
            if (math.isnan(float(items[5]))): continue
            if (math.isnan(float(items[6]))): continue
            if (math.isnan(float(items[7]))): continue
            if (math.isnan(float(items[8]))): continue
            
            mass = items[0]
            mass = re.sub(",gamma.*,","",mass)
            mass = re.sub(",gamma.*","",mass)
            mass = re.sub(".*mH:","",mass)
            if (("gamma" in items[0]) and (float(mass)<400)): continue

            mass_list.append(float(mass))
            down_2sig_list.append(float(items[3]))
            down_1sig_list.append(float(items[4]))
            exp_list.append(float(items[5]))
            up_1sig_list.append(float(items[6]))
            up_2sig_list.append(float(items[7]))
            obs_list.append(float(items[8]))

  #sort according to mass
    mass_list , down_2sig_list, down_1sig_list, exp_list, up_1sig_list, up_2sig_list, obs_list = (list(t) for t in zip(*sorted(zip(mass_list , down_2sig_list, down_1sig_list, exp_list, up_1sig_list, up_2sig_list, obs_list))))
    
    return get_graph_from_list(mass_list, obs_list, exp_list,
                               up_1sig_list, up_2sig_list,
                               down_1sig_list, down_2sig_list)

def plot_limit(path, file_name, xslabel, typelabel, minx=200, maxx=2000):

    low_y = 0.001
    hi_y = 5
    #low_y = 0.5
    #hi_y = 50
    sigma = "#sigma_{"+xslabel+"}"
    unit = "95% CL limits on "+sigma+" #times BR(S#rightarrow ZZ) [pb]"
    #unit = "95% CL limits on "+sigma+" #times BR(S#rightarrow ZZ #rightarrow 4l) [fb]"
    x_axis_title = "m_{S} [GeV]"

    dummy=ROOT.TH2F("dummy",";"+";"+unit,
                    160, float(minx),float(maxx),3000,low_y,hi_y);
    #dummy.GetXaxis().SetNdivisions(8);

    hist_obs,hist_exp,hist_1s,hist_2s = make_limit_graph(path+"DNNnop_v6_ggF.txt")
    hist_obs2,hist_exp2,hist_1s2,hist_2s2 = make_limit_graph(path+"DNN2020_v6_ggF.txt")
    for i in range(hist_exp.GetN()): 
	hist_exp.GetY()[i] *= 1./4.52
	hist_exp2.GetY()[i] *= 1./4.52
	hist_obs.GetY()[i] *= 1./4.52
	hist_obs2.GetY()[i] *= 1./4.52

    canvas = ROOT.TCanvas("canvas2", " ", 600, 600)
    canvas.SetLogy()

    ##========================= add ratio pad ===========================
    ratio=0.2
    if ratio>0:
      fraction=ratio+0.2
      Pad1 = TPad("p1","p1",0,fraction*1.0/(fraction+1),1,1,0,0) # x1,y1,x2,y2
      Pad1.SetMargin(0.15,0.10,0.03,0.05)
      Pad1.SetLogy()
      Pad2 = TPad("p2","p2",0,0,1,fraction*1.0/(fraction+1),0,0)
      Pad2.SetMargin(0.15,0.10,0.15/fraction,0.04)
      Pad2.SetGrid()
      Pad1.Draw()
      Pad2.Draw()

    hist_exp.SetLineWidth(3)
    hist_exp.SetMarkerStyle(20)
    hist_exp.SetMarkerSize(0.5)
    hist_exp.SetLineColor(1)
    hist_exp.SetLineStyle(1)
    hist_exp2.SetLineWidth(3)
    hist_exp2.SetMarkerStyle(20)
    hist_exp2.SetMarkerSize(0.5)
    hist_exp2.SetLineColor(2)

    if ratio>0: Pad1.cd()
    else: canvas.cd()
    dummy.Draw()
    hist_exp. Draw("CP")
    hist_exp2. Draw("CP")

    legend = ROOT.myLegend(0.56, 0.70, 0.83, 0.90)
    legend.AddEntry(hist_exp,  "Fix other POI", "l")
    legend.AddEntry(hist_exp2, "Float", "l")
    legend.Draw()

    lumi = 138.97
    x_off_title = 0.20
    ROOT.myText(x_off_title, 0.85, 1, "#bf{#it{ATLAS}} Internal")
    ROOT.myText(x_off_title, 0.80, 1, "13 TeV, {:.1f} fb^{{-1}}".format(lumi))
    ROOT.myText(x_off_title, 0.75, 1, typelabel)
    #dummy.Draw("AXIS SAME")

    if ratio>0: 
      x1=hist_exp.GetX()
      y1=hist_exp.GetY()
      y2=hist_exp2.GetY()
      hist_ratio2 =  ROOT.TGraph(len(x1))
      for i in range(len(x1)):
         hist_ratio2.SetPoint(i, x1[i], y2[i]/y1[i])
      Pad2.cd()
      mg = ROOT.TMultiGraph()
      mg.GetXaxis().SetLimits(minx,maxx)
      mg.SetMinimum(0.8)
      mg.SetMaximum(1.2)
      hist_ratio2.SetMarkerSize(0.5)
      hist_ratio2.SetLineColor(2)
      hist_ratio2.SetLineStyle(2)
      hist_ratio2.SetLineWidth(2)
      mg.Add(hist_ratio2)
      mg.Draw("ALP")
      mg.GetYaxis().SetTitle("Ratio")
      mg.GetXaxis().SetTitle(x_axis_title)

    canvas.SaveAs("Pdf_plots/"+file_name+".pdf")


if __name__ == "__main__":
    #path="/afs/cern.ch/work/h/hezhu/public/workplace/H4lAna/CMakeWS/Condor/Comb_DNN_ggF/"
    path="/afs/cern.ch/work/h/hezhu/public/workplace/H4lAna/CMakeWS/Condor/DNN_NWA_ggF/"
    plot_limit(path,"ratio_4lexp_profilepoi_ggF","ggF","NWA, ggF production")

