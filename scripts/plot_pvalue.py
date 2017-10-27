#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch()

import AtlasStyle
from array import array

import sys
if not hasattr(ROOT, "myText"):
    ROOT.gROOT.LoadMacro("$HZZWSCODEDIR/scripts/loader.c")
    ROOT.gROOT.LoadMacro("$HZZWSCODEDIR/scripts/AtlasUtils.C")

def create_TGraph(x_, nominal_):
    gr = ROOT.TGraph(len(x_), array('f', x_), array('f', nominal_))
    return gr

def get_p0_graph(file_name):
    mass_list = []
    pvalue_list = []
    min_pvalue = 999
    min_mass = -999
    with open(file_name, 'r') as f:
        for line in f:
            items = line.split()
            #find the mass
            print items
            mass = float(items[0].split(',')[0].split(':')[1])
            pvalue = float(items[1])

            mass_list.append(mass)
            pvalue_list.append(pvalue)
            if pvalue < min_pvalue:
                min_mass = mass
                min_pvalue = pvalue

    print "max significance: ",ROOT.RooStats.PValueToSignificance(min_pvalue),min_mass,"GeV"
    return create_TGraph(mass_list, pvalue_list)

def select(file_name, gamma_frac):
    out_lines = " "
    with open(file_name, 'r') as f:
        for line in f:
            items = line.split()
            #find mass and gamma
            mass = float(items[0].split(',')[0].split(':')[1])
            gamma = float(items[0].split(',')[1].split(':')[1])
            #print mass, gamma, (mass*gamma_frac - gamma)
            if abs(mass*gamma_frac - gamma) < 1E-1:
                out_lines += line

    new_file_name = file_name.replace(".txt", "_"+str(gamma_frac)+".txt")
    with open(new_file_name, 'w') as f:
        f.write(out_lines)

    return new_file_name

def plot_p0( file_name, out_name, fraction, lumi):
    #file_name = select(file_name, fraction)

    cms = 13
    yoffset = 0.35
    gr_obs = get_p0_graph(file_name)

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    canvas.SetLogy()

    dummy = ROOT.TH1F("dummy",";m_{S} [GeV];Local p_{0}",
                      800, 200., 1000.)
    dummy.GetYaxis().SetRangeUser(1E-5, 1.2)
    dummy.Draw()

    ## setup sigma
    p7_sigma = 2.419637e-01
    one_sigma = 1.58655253931457074e-01
    two_sigma = (1.-0.954499736104)/2
    three_sigma = (1.-0.997300203937)/2
    four_sigma = 3.16712418331199785e-05
    ROOT.AddLine(dummy, p7_sigma, 11, 2)
    ROOT.AddLine(dummy, one_sigma, 11, 2)
    ROOT.AddLine(dummy, two_sigma, 11, 2)
    ROOT.AddLine(dummy, three_sigma, 11, 2)
    ROOT.AddLine(dummy, four_sigma, 11, 2)
    latex = ROOT.TLatex()
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.SetTextColor(ROOT.kGray);
    latex.SetTextAlign(12);
    latex.DrawLatex(1000., one_sigma, "1#sigma")
    latex.DrawLatex(1000., two_sigma, "2#sigma")
    latex.DrawLatex(1000., three_sigma, "3#sigma")
    latex.DrawLatex(1000., four_sigma, "4#sigma")
    latex.Draw()

    #plot observed p0
    gr_obs.Draw("L")
    gr_obs.SetLineWidth(2)
    legend = ROOT.myLegend(0.65, yoffset-0.05, 0.8, yoffset)
    legend.AddEntry(gr_obs, "observed", "L")
    legend.Draw("same")

    ROOT.myText(0.2, yoffset, 1, "#bf{#it{ATLAS}} Internal")
    ROOT.myText(0.2, yoffset-0.05, 1, str(cms)+" TeV, "+str(lumi)+" fb^{-1}")
    ROOT.myText(0.2, yoffset-0.05*2, 1, "LWA {:.2f}%".format(fraction*100))
    canvas.SaveAs("p0_"+out_name+".eps")
    canvas.SaveAs("p0_"+out_name+".png")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print sys.argv[0],"file_name out_name fraction [lumi]"
        sys.exit(1)

    file_name = sys.argv[1]
    out_name = sys.argv[2]
    fraction = float(sys.argv[3])

    if len(sys.argv) > 4:
        lumi = float(sys.argv[4])
    else:
        lumi = 14.8
        #lumi = 9.5

    plot_p0(file_name, out_name, fraction, lumi)
