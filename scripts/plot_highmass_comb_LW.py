#!/usr/bin/env python
import ROOT
import AtlasStyle
import math
import re
import os
from array import array
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
    gr_obs.SetLineWidth(3)
    gr_obs.SetLineStyle(1)
    gr_obs.SetMarkerStyle(20)
    gr_obs.SetMarkerSize(0.6)

    gr_exp.SetLineWidth(2)
    gr_exp.SetLineStyle(2)
    gr_exp.SetMarkerSize(0)

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

def plot_limit(path, file_name, xslabel, typelabel, minx=400,maxx=2000):

    low_y = 0.001
    hi_y = 5
    sigma = "#sigma_{"+xslabel+"}"
    unit = "95% CL limits on "+sigma+" #times BR(S#rightarrow ZZ) [pb]"
    x_axis_title = "m_{S} [GeV]"

    dummy=ROOT.TH2F("dummy",";"+x_axis_title+";"+unit,
                    160, float(minx),float(maxx),3000,low_y,hi_y);
    dummy.GetXaxis().SetNdivisions(8);

    hist_obs,hist_exp,hist_1s,hist_2s = make_limit_graph(path+"4l_LWA_w15.txt")
    hist_obs2,hist_exp2,hist_1s2,hist_2s2 = make_limit_graph(path+"llvv_LWA_w15.txt")
    hist_obs3,hist_exp3,hist_1s3,hist_2s3 = make_limit_graph(path+"comb_LWA_w15.txt")
    #hist_exp.Scale(1./4.52)
    for i in range(hist_exp.GetN()): 
	hist_exp.GetY()[i] *= 1./4.52
	hist_obs.GetY()[i] *= 1./4.52

    canvas = ROOT.TCanvas("canvas2", " ", 600, 600)
    canvas.SetLogy()

    dummy.Draw()
    #hist_obs.SetLineWidth(3)
    hist_obs.SetLineColor(4)
    #hist_obs2.SetLineWidth(3)
    hist_obs2.SetLineColor(2)
    #hist_obs3.SetLineWidth(3)
    hist_obs3.SetLineColor(1)
    hist_exp.SetLineColor(4)
    hist_exp2.SetLineColor(2)
    hist_exp3.SetLineColor(1)

    hist_2s3.Draw("3")
    hist_1s3.Draw("3")
    hist_exp3. Draw("LP")
    hist_exp. Draw("LP")
    hist_exp2. Draw("LP")
    hist_obs3. Draw("CP")
    #hist_obs. Draw("CP")
    #hist_obs2. Draw("CP")

    dummy.Draw("AXIS SAME")

    #legend = ROOT.myLegend(0.56, 0.60, 0.83, 0.90)
    legend = ROOT.myLegend(0.56, 0.70, 0.83, 0.88)
    legend.AddEntry(hist_exp,  "llll", "l")
    legend.AddEntry(hist_exp2, "llvv", "l")
    legend.AddEntry(hist_exp3, "Combined", "l")
    legend.AddEntry(hist_1s, "Expected #pm 1 #sigma", "f")
    legend.AddEntry(hist_2s, "Expected #pm 2 #sigma", "f")
    legend.Draw()

    lumi = 139.0
    x_off_title = 0.20
    ROOT.myText(x_off_title, 0.85, 1, "#bf{#it{ATLAS}} Internal")
    ROOT.myText(x_off_title, 0.80, 1, "13 TeV, {:.1f} fb^{{-1}}".format(lumi))
    ROOT.myText(x_off_title, 0.75, 1, typelabel)

    #file_name = re.sub("txt","pdf",file_name)
    canvas.SaveAs("Pdf_plots/"+file_name+".pdf")


if __name__ == "__main__":
    path="/afs/cern.ch/work/h/hezhu/public/workplace/H4lAna/CMakeWS/Condor/Comb_LWA/"
    plot_limit(path,"LWA_comb_4l_llvv_w15","ggF","LWA, #Gamma_{S} = 0.15 #times m_{S}")

